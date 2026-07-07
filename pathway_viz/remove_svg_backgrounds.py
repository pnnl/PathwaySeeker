#!/usr/bin/env python3
"""
Remove backgrounds from SVG files and crop to content.

The SVGs exported from the pathway viewer have this structure:
  <svg>
    <g class="zoom-g" transform="translate(tx,ty) scale(s)">
      <g class="canvas-group">
        <rect id="canvas" width="W" height="H" .../>   ← background to remove
        ...
      </g>
      <g id="reactions"> ... </g>
      <g id="nodes"> ... </g>
    </g>
  </svg>

Strategy:
  1. Parse the outer zoom-g transform (tx, ty, scale).
  2. Walk all visible content elements (circles, paths, text, images, rects
     that are NOT the canvas/mouse-node/resize-rect) and collect their
     bounding boxes in zoom-g local coordinates.
  3. Apply the zoom-g transform to get SVG-space bounds.
  4. Set viewBox + width/height to those bounds (+ small padding).
  5. Make the canvas rect transparent (fill:none, no stroke).
"""

import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

SVG_NS = 'http://www.w3.org/2000/svg'
XLINK_NS = 'http://www.w3.org/1999/xlink'

# Tags we skip entirely when computing bounds
SKIP_TAGS = {'defs', 'style', 'title', 'desc', 'script'}

# IDs whose elements (and children) we skip
SKIP_IDS = {'mouse-node', 'canvas', 'brush-container',
            'beziers', 'text-labels', 'direction-arrow-container'}

# Classes that indicate invisible/UI-only elements
SKIP_CLASSES = {'resize-rect', 'direction-arrow', 'stoichiometry-label-rect'}


def _parse_zoom_g_transform(transform_str):
    """
    Parse 'translate(tx,ty) scale(s)' or 'translate(tx ty) scale(s)'.
    Returns (tx, ty, scale).
    """
    tx, ty, scale = 0.0, 0.0, 1.0
    t_match = re.search(r'translate\(\s*([+-]?[\d.eE+-]+)[,\s]+([+-]?[\d.eE+-]+)\s*\)',
                        transform_str)
    if t_match:
        tx = float(t_match.group(1))
        ty = float(t_match.group(2))
    s_match = re.search(r'scale\(\s*([+-]?[\d.eE+-]+)\s*\)', transform_str)
    if s_match:
        scale = float(s_match.group(1))
    return tx, ty, scale


def _parse_translate(transform_str):
    """Parse 'translate(x,y)' or 'translate(x y)' from a transform string."""
    m = re.search(r'translate\(\s*([+-]?[\d.eE+-]+)[,\s]+([+-]?[\d.eE+-]+)\s*\)',
                  transform_str)
    if m:
        return float(m.group(1)), float(m.group(2))
    return 0.0, 0.0


def _tag(elem):
    return elem.tag.replace('{' + SVG_NS + '}', '')


def _get_classes(elem):
    return set(elem.get('class', '').split())


def _should_skip(elem):
    tag = _tag(elem)
    if tag in SKIP_TAGS:
        return True
    eid = elem.get('id', '')
    if eid in SKIP_IDS:
        return True
    classes = _get_classes(elem)
    if classes & SKIP_CLASSES:
        return True
    style = elem.get('style', '')
    if 'display: none' in style or 'visibility: hidden' in style:
        return True
    return False


def _collect_bounds(elem, ox=0.0, oy=0.0, bounds_list=None):
    """
    Recursively walk elem, accumulating (x1,y1,x2,y2) tuples in bounds_list.
    ox, oy: accumulated translation offset from parent transforms.
    """
    if bounds_list is None:
        bounds_list = []

    if _should_skip(elem):
        return bounds_list

    # Accumulate this element's own translate
    transform = elem.get('transform', '')
    dx, dy = _parse_translate(transform)
    cx, cy = ox + dx, oy + dy

    tag = _tag(elem)

    if tag == 'circle':
        try:
            ex = float(elem.get('cx', 0)) + cx
            ey = float(elem.get('cy', 0)) + cy
            r  = float(elem.get('r', 0))
            bounds_list.append((ex - r, ey - r, ex + r, ey + r))
        except (ValueError, TypeError):
            pass

    elif tag == 'rect':
        # Skip canvas/mouse-node/resize rects (already handled by SKIP_IDS/SKIP_CLASSES)
        try:
            rx = float(elem.get('x', 0)) + cx
            ry = float(elem.get('y', 0)) + cy
            rw = float(elem.get('width', 0))
            rh = float(elem.get('height', 0))
            if rw > 0 and rh > 0:
                bounds_list.append((rx, ry, rx + rw, ry + rh))
        except (ValueError, TypeError):
            pass

    elif tag == 'text':
        try:
            tx = float(elem.get('x', 0)) + cx
            ty = float(elem.get('y', 0)) + cy
            # Rough text bounds — use font-size from style if available
            style = elem.get('style', '')
            fs_match = re.search(r'font-size:\s*([\d.]+)px', style)
            fs = float(fs_match.group(1)) if fs_match else 12.0
            text_content = elem.text or ''
            char_w = fs * 0.6
            text_w = len(text_content) * char_w
            bounds_list.append((tx - text_w / 2, ty - fs, tx + text_w / 2, ty + fs * 0.3))
        except (ValueError, TypeError):
            pass

    elif tag == 'image':
        try:
            ix = float(elem.get('x', 0)) + cx
            iy = float(elem.get('y', 0)) + cy
            iw = float(elem.get('width', 0))
            ih = float(elem.get('height', 0))
            if iw > 0 and ih > 0:
                bounds_list.append((ix, iy, ix + iw, iy + ih))
        except (ValueError, TypeError):
            pass

    elif tag == 'path':
        # Parse M/C/L commands to get rough bounds
        d = elem.get('d', '')
        coords = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', d)
        pts = []
        try:
            for i in range(0, len(coords) - 1, 2):
                pts.append((float(coords[i]) + cx, float(coords[i+1]) + cy))
        except (ValueError, IndexError):
            pass
        if pts:
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            bounds_list.append((min(xs), min(ys), max(xs), max(ys)))

    # Recurse into children (passing accumulated cx, cy)
    for child in elem:
        _collect_bounds(child, cx, cy, bounds_list)

    return bounds_list


def remove_background_from_svg(svg_path, output_path=None):
    """
    Remove background and crop SVG to content.
    """
    try:
        ET.register_namespace('', SVG_NS)
        ET.register_namespace('xlink', XLINK_NS)

        tree = ET.parse(svg_path)
        root = tree.getroot()

        # ── 1. Find the zoom-g group and parse its transform ──────────────
        zoom_g = root.find('.//{' + SVG_NS + '}g[@class="zoom-g"]')
        if zoom_g is None:
            # Fallback: try first <g> child
            zoom_g = root.find('{' + SVG_NS + '}g')

        tx, ty, scale = 0.0, 0.0, 1.0
        if zoom_g is not None:
            zt = zoom_g.get('transform', '')
            tx, ty, scale = _parse_zoom_g_transform(zt)

        # ── 2. Make canvas rect transparent ───────────────────────────────
        for rect in root.findall('.//{' + SVG_NS + '}rect[@id="canvas"]'):
            rect.set('style', 'fill: none; stroke: none;')
        for rect in root.findall('.//{' + SVG_NS + '}rect[@id="mouse-node"]'):
            rect.set('style', 'display: none;')

        # ── 3. Collect content bounds in zoom-g local coordinates ──────────
        bounds_list = []
        if zoom_g is not None:
            for child in zoom_g:
                _collect_bounds(child, 0.0, 0.0, bounds_list)
        else:
            _collect_bounds(root, 0.0, 0.0, bounds_list)

        if not bounds_list:
            print(f"  Warning: no content bounds found in {svg_path}")
            if output_path:
                tree.write(output_path, encoding='utf-8', xml_declaration=True)
            return False

        # ── 4. Compute tight bounds in zoom-g local space ─────────────────
        min_x = min(b[0] for b in bounds_list)
        min_y = min(b[1] for b in bounds_list)
        max_x = max(b[2] for b in bounds_list)
        max_y = max(b[3] for b in bounds_list)

        padding = 20  # local-space padding

        lx = min_x - padding
        ly = min_y - padding
        lw = (max_x - min_x) + 2 * padding
        lh = (max_y - min_y) + 2 * padding

        # ── 5. Convert to SVG space using zoom-g transform ─────────────────
        # SVG point = zoom-g point * scale + (tx, ty)
        vb_x = lx * scale + tx
        vb_y = ly * scale + ty
        vb_w = lw * scale
        vb_h = lh * scale

        root.set('viewBox', f'{vb_x:.2f} {vb_y:.2f} {vb_w:.2f} {vb_h:.2f}')
        root.set('width',  f'{vb_w:.2f}')
        root.set('height', f'{vb_h:.2f}')

        if output_path is None:
            output_path = svg_path

        tree.write(output_path, encoding='utf-8', xml_declaration=True)
        return True

    except Exception as e:
        print(f"Error processing {svg_path}: {e}")
        import traceback
        traceback.print_exc()
        return False


def process_figs_folder(figs_folder='new_svgs', output_folder=None):
    """
    Process all SVG files in figs_folder and save to output_folder.
    """
    figs_path = Path(figs_folder)

    if not figs_path.exists():
        print(f"Error: Folder '{figs_folder}' not found!")
        return

    if output_folder is None:
        output_folder = f"{figs_folder}_no_background"

    output_path = Path(output_folder)
    output_path.mkdir(exist_ok=True)
    print(f"Output folder: {output_folder}")

    svg_files = list(figs_path.glob('*.svg'))

    if not svg_files:
        print(f"No SVG files found in '{figs_folder}'")
        return

    print(f"Found {len(svg_files)} SVG file(s) in '{figs_folder}'")
    print("-" * 60)

    processed = 0
    modified  = 0

    for svg_file in svg_files:
        print(f"Processing: {svg_file.name}")
        output_file = output_path / svg_file.name
        was_modified = remove_background_from_svg(svg_file, output_file)
        if was_modified:
            print(f"  ✓ Cropped & background removed → {output_file.name}")
            modified += 1
        else:
            print(f"  → Saved to {output_file.name}")
        processed += 1
        print()

    print("-" * 60)
    print(f"Summary: Processed {processed} files, modified {modified}, "
          f"saved to '{output_folder}'")
    print(f"Original files remain unchanged in '{figs_folder}'")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Remove backgrounds from SVG files and crop to content"
    )
    parser.add_argument(
        '--folder',
        default='new_svgs',
        help='Path to folder containing SVG files (default: new_svgs)'
    )
    parser.add_argument(
        '--output',
        default=None,
        help='Output folder path (default: <input_folder>_no_background)'
    )

    args = parser.parse_args()
    process_figs_folder(args.folder, args.output)
