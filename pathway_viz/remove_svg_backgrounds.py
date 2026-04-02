#!/usr/bin/env python3
"""
Remove backgrounds from SVG files in the FIGS folder.
This script processes all SVG files and removes any rectangle elements
that appear to be backgrounds (typically large rectangles with fill colors).
"""

import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

# Define namespaces
SVG_NS = 'http://www.w3.org/2000/svg'
NAMESPACES = {'svg': SVG_NS}


def get_element_bounds(elem, parent_transform=(0, 0)):
    """Extract bounds from an SVG element, accounting for transforms."""
    bounds = None
    
    # Parse transform if present
    transform = elem.get('transform', '')
    tx, ty = parent_transform
    if 'translate' in transform:
        match = re.search(r'translate\(([^,]+),([^)]+)\)', transform)
        if match:
            tx += float(match.group(1))
            ty += float(match.group(2))
    
    # Get bounds based on element type
    tag = elem.tag.replace('{' + SVG_NS + '}', '')
    
    if tag == 'rect':
        x = float(elem.get('x', 0)) + tx
        y = float(elem.get('y', 0)) + ty
        width = float(elem.get('width', 0))
        height = float(elem.get('height', 0))
        bounds = (x, y, x + width, y + height)
        
    elif tag == 'circle':
        cx = float(elem.get('cx', 0)) + tx
        cy = float(elem.get('cy', 0)) + ty
        r = float(elem.get('r', 0))
        bounds = (cx - r, cy - r, cx + r, cy + r)
        
    elif tag == 'text':
        x = float(elem.get('x', 0)) + tx
        y = float(elem.get('y', 0)) + ty
        # Approximate text bounds (rough estimate)
        bounds = (x - 50, y - 20, x + 200, y + 20)
        
    elif tag == 'image':
        x = float(elem.get('x', 0)) + tx
        y = float(elem.get('y', 0)) + ty
        width = float(elem.get('width', 0))
        height = float(elem.get('height', 0))
        bounds = (x, y, x + width, y + height)
    
    # Account for stroke width
    style = elem.get('style', '')
    stroke_match = re.search(r'stroke-width:\s*([0-9.]+)', style)
    if stroke_match and bounds:
        stroke_width = float(stroke_match.group(1)) / 2
        bounds = (bounds[0] - stroke_width, bounds[1] - stroke_width,
                 bounds[2] + stroke_width, bounds[3] + stroke_width)
    
    return bounds, (tx, ty)


def get_content_bounds(root):
    """Find the bounding box of all visible content."""
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = float('-inf'), float('-inf')
    
    # Skip invisible elements
    invisible_tags = ['defs', 'style', 'title', 'desc']
    
    def process_element(elem, parent_transform=(0, 0)):
        nonlocal min_x, min_y, max_x, max_y
        
        tag = elem.tag.replace('{' + SVG_NS + '}', '')
        
        # Skip invisible or background elements
        if tag in invisible_tags or elem.get('id') in ['mouse-node', 'canvas']:
            return parent_transform
            
        # Check if element has visibility
        style = elem.get('style', '')
        if 'display: none' in style or 'visibility: hidden' in style:
            return parent_transform
        
        # Get bounds for this element
        bounds, current_transform = get_element_bounds(elem, parent_transform)
        
        if bounds:
            min_x = min(min_x, bounds[0])
            min_y = min(min_y, bounds[1])
            max_x = max(max_x, bounds[2])
            max_y = max(max_y, bounds[3])
        
        # Process children
        for child in elem:
            process_element(child, current_transform)
        
        return current_transform
    
    process_element(root)
    
    # Add small padding
    padding = 10
    if min_x != float('inf'):
        return (min_x - padding, min_y - padding, 
                max_x - min_x + 2 * padding, max_y - min_y + 2 * padding)
    
    return None


def remove_background_from_svg(svg_path, output_path=None):
    """
    Remove background elements from an SVG file.
    
    Args:
        svg_path: Path to the input SVG file
        output_path: Path for the output file (if None, overwrites original)
    """
    try:
        # Register namespace
        ET.register_namespace('', SVG_NS)
        
        # Parse SVG
        tree = ET.parse(svg_path)
        root = tree.getroot()
        
        # Remove canvas and mouse-node rectangles
        for rect in root.findall('.//{' + SVG_NS + '}rect[@id="canvas"]'):
            # Change fill to none instead of removing
            style = rect.get('style', '')
            style = re.sub(r'fill:\s*rgb\([^)]+\);?', 'fill: none;', style)
            style = re.sub(r'stroke:\s*rgb\([^)]+\);?', '', style)
            style = re.sub(r'stroke-width:\s*[^;]+;?', '', style)
            rect.set('style', style)
        
        for rect in root.findall('.//{' + SVG_NS + '}rect[@id="mouse-node"]'):
            rect.getparent().remove(rect) if hasattr(rect, 'getparent') else None
        
        # Get content bounds
        bounds = get_content_bounds(root)
        
        if bounds:
            # Set viewBox to crop to content
            root.set('viewBox', f'{bounds[0]:.2f} {bounds[1]:.2f} {bounds[2]:.2f} {bounds[3]:.2f}')
            root.set('width', f'{bounds[2]:.2f}')
            root.set('height', f'{bounds[3]:.2f}')
        
        # Save
        if output_path is None:
            output_path = svg_path
        
        tree.write(output_path, encoding='utf-8', xml_declaration=True)
        
        return bounds is not None
        
    except Exception as e:
        print(f"Error processing {svg_path}: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def process_figs_folder(figs_folder='FIGS', output_folder=None):
    """
    Process all SVG files in the FIGS folder and save to output folder.
    
    Args:
        figs_folder: Path to the FIGS folder (default: 'FIGS')
        output_folder: Path to output folder (default: figs_folder + '_no_background')
    """
    figs_path = Path(figs_folder)
    
    if not figs_path.exists():
        print(f"Error: Folder '{figs_folder}' not found!")
        return
    
    # Set default output folder
    if output_folder is None:
        output_folder = f"{figs_folder}_no_background"
    
    output_path = Path(output_folder)
    
    # Create output folder if it doesn't exist
    output_path.mkdir(exist_ok=True)
    print(f"Output folder: {output_folder}")
    
    # Find all SVG files
    svg_files = list(figs_path.glob('*.svg'))
    
    if not svg_files:
        print(f"No SVG files found in '{figs_folder}'")
        return
    
    print(f"Found {len(svg_files)} SVG file(s) in '{figs_folder}'")
    print("-" * 60)
    
    processed = 0
    modified = 0
    
    for svg_file in svg_files:
        print(f"Processing: {svg_file.name}")
        
        # Define output path (same filename in new folder)
        output_file = output_path / svg_file.name
        
        # Remove background and save to output folder
        was_modified = remove_background_from_svg(svg_file, output_file)
        
        if was_modified:
            print(f"  ✓ Background removed → {output_file.name}")
            modified += 1
        else:
            print(f"  → Saved to {output_file.name}")
        
        processed += 1
        print()
    
    print("-" * 60)
    print(f"Summary: Processed {processed} files, saved {processed} files to '{output_folder}'")
    print(f"Original files remain unchanged in '{figs_folder}'")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Remove backgrounds from SVG files and save to a new folder"
    )
    parser.add_argument(
        '--folder',
        default='FIGS',
        help='Path to folder containing SVG files (default: FIGS)'
    )
    parser.add_argument(
        '--output',
        default=None,
        help='Output folder path (default: <input_folder>_no_background)'
    )
    
    args = parser.parse_args()
    
    process_figs_folder(args.folder, args.output)
