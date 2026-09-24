"""
PubChem Structure Image Downloader
Downloads molecular structure images from PubChem (public domain alternative to KEGG)
- Maps KEGG IDs to PubChem CIDs via PUG REST API
- Downloads structure images (PNG format)
- Processes images for transparency
- Public domain images, no copyright restrictions
"""
import os
import json
import requests
from PIL import Image
from io import BytesIO
import time

# =============================================================================
# CONFIGURATION
# =============================================================================

OUTPUT_DIR             = './static/structure_imgs'
FAILED_DOWNLOADS_FILE  = './static/failed_downloads.json'
MAPPING_CACHE_FILE     = './static/kegg_pubchem_mapping.json'
IMAGE_DIMENSIONS_FILE  = os.path.join(OUTPUT_DIR, 'image_dimensions.json')

# PubChem API endpoints
PUBCHEM_CID_LOOKUP   = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON'
PUBCHEM_IMAGE_URL    = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/PNG?image_size=500x500'
KEGG_TO_PUBCHEM_API  = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{}/cids/JSON'

TRANSPARENT_COLORS = [(255, 255, 255), (245, 245, 245)]
IMAGE_FORMAT       = 'PNG'
MAX_RETRIES        = 3
RETRY_DELAY        = 2    # seconds between retries
REQUEST_DELAY      = 0.3  # seconds between requests (PubChem rate limit)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def setup_output_directory():
    """Create output directory if it doesn't exist."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(MAPPING_CACHE_FILE), exist_ok=True)


def load_mapping_cache():
    """Load cached KEGG → PubChem CID mappings."""
    if not os.path.exists(MAPPING_CACHE_FILE):
        return {}
    try:
        with open(MAPPING_CACHE_FILE) as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        print('Warning: Could not read mapping cache, starting fresh')
        return {}


def save_mapping_cache(mapping):
    """Save KEGG → PubChem CID mappings to disk."""
    try:
        with open(MAPPING_CACHE_FILE, 'w') as f:
            json.dump(mapping, f, indent=2)
    except IOError as e:
        print(f'Warning: Could not save mapping cache: {e}')


def load_failed_downloads():
    """Load previously failed download IDs from tracking file."""
    if not os.path.exists(FAILED_DOWNLOADS_FILE):
        return {'failed': [], 'last_updated': None}
    try:
        with open(FAILED_DOWNLOADS_FILE) as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        return {'failed': [], 'last_updated': None}


def save_failed_downloads(failed_list):
    """Save list of failed download IDs to tracking file."""
    data = {
        'failed':       failed_list,
        'last_updated': time.strftime('%Y-%m-%d %H:%M:%S'),
    }
    try:
        with open(FAILED_DOWNLOADS_FILE, 'w') as f:
            json.dump(data, f, indent=2)
    except IOError as e:
        print(f'Warning: Could not save failed downloads file: {e}')


def get_existing_images_cache():
    """Return a set of KEGG IDs that already have a downloaded PNG."""
    if not os.path.exists(OUTPUT_DIR):
        return set()
    try:
        return {
            f.replace('.png', '')
            for f in os.listdir(OUTPUT_DIR)
            if f.endswith('.png')
        }
    except Exception as e:
        print(f'Warning: Could not read existing images directory: {e}')
        return set()


def load_pathway_data(json_path):
    """Load Escher pathway JSON from disk."""
    try:
        with open(json_path) as f:
            return json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f'JSON file not found: {json_path}')
    except json.JSONDecodeError as e:
        raise ValueError(f'Invalid JSON format in {json_path}: {e}')


def extract_kegg_compounds(pathway_data):
    """
    Extract KEGG compound IDs and names from Escher pathway JSON.

    Returns:
        dict: { kegg_id: compound_name }
    """
    compounds = {}
    if len(pathway_data) > 1 and 'nodes' in pathway_data[1]:
        for node_id, info in pathway_data[1]['nodes'].items():
            bigg_id = info.get('bigg_id')
            name    = info.get('name', 'Unknown')
            if bigg_id and bigg_id.startswith('C') and len(bigg_id) >= 5:
                compounds[bigg_id] = name.replace(';', '').strip()
    return compounds


def _find_latest_output_json(json_dir):
    """
    Find the most recently modified *_output.json in json_dir.

    Returns:
        str | None: Full path to the file, or None if nothing found.
    """
    if not os.path.exists(json_dir):
        return None

    matches = [
        f for f in os.listdir(json_dir)
        if f.endswith('_output.json') and not f.startswith('.')
    ]
    if not matches:
        return None

    matches.sort(
        key=lambda f: os.path.getmtime(os.path.join(json_dir, f)),
        reverse=True,
    )
    return os.path.join(json_dir, matches[0])


# =============================================================================
# KEGG → PUBCHEM MAPPING
# =============================================================================

def kegg_to_pubchem_cid(kegg_id, mapping_cache):
    """
    Convert a KEGG compound ID to a PubChem CID.
    Results are cached in mapping_cache to avoid repeat API calls.

    Args:
        kegg_id       (str):  KEGG compound ID e.g. 'C00031'
        mapping_cache (dict): In-memory cache (mutated in place)

    Returns:
        str | None: PubChem CID string, or None if not found.
    """
    # Return cached result (including None for known failures)
    if kegg_id in mapping_cache:
        return mapping_cache[kegg_id]

    time.sleep(REQUEST_DELAY)

    # Method 1: KEGG cross-reference lookup
    try:
        resp = requests.get(KEGG_TO_PUBCHEM_API.format(kegg_id), timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                cid = str(data['IdentifierList']['CID'][0])
                mapping_cache[kegg_id] = cid
                return cid
    except Exception:
        pass

    # Method 2: Search by KEGG ID as synonym
    time.sleep(REQUEST_DELAY)
    try:
        resp = requests.get(PUBCHEM_CID_LOOKUP.format(kegg_id), timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                cid = str(data['IdentifierList']['CID'][0])
                mapping_cache[kegg_id] = cid
                return cid
    except Exception:
        pass

    # Both methods failed — cache the failure so we don't retry
    mapping_cache[kegg_id] = None
    return None


# =============================================================================
# IMAGE PROCESSING
# =============================================================================

def make_background_transparent(image):
    """Make white and near-white backgrounds transparent."""
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    pixel_data = []
    for pixel in image.getdata():
        if pixel[:3] in TRANSPARENT_COLORS:
            pixel_data.append((0, 0, 0, 0))
        else:
            pixel_data.append(pixel)
    image.putdata(pixel_data)
    return image


def crop_to_content(image, padding=4):
    """
    Crop transparent padding around visible image content.

    Args:
        image   (PIL.Image): RGBA image
        padding (int):       Pixels of padding to retain around content

    Returns:
        PIL.Image: Cropped image
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')

    alpha = image.split()[3]
    bbox  = alpha.getbbox()

    if bbox is None:
        return image  # Fully transparent — return as-is

    left   = max(0, bbox[0] - padding)
    top    = max(0, bbox[1] - padding)
    right  = min(image.width,  bbox[2] + padding)
    bottom = min(image.height, bbox[3] + padding)

    return image.crop((left, top, right, bottom))


def download_pubchem_image(cid, attempt=1):
    """
    Download a structure PNG from PubChem with retry logic.

    Args:
        cid     (str): PubChem CID
        attempt (int): Current attempt number (internal, used for recursion)

    Returns:
        PIL.Image | None
    """
    url = PUBCHEM_IMAGE_URL.format(cid)
    try:
        time.sleep(REQUEST_DELAY)
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return Image.open(BytesIO(resp.content))
    except requests.Timeout:
        if attempt < MAX_RETRIES:
            print(f'⏳ Timeout (attempt {attempt}/{MAX_RETRIES}), retrying... ', end='')
            time.sleep(RETRY_DELAY)
            return download_pubchem_image(cid, attempt + 1)
        return None
    except requests.RequestException:
        if attempt < MAX_RETRIES:
            time.sleep(RETRY_DELAY)
            return download_pubchem_image(cid, attempt + 1)
        return None
    except Exception as e:
        print(f'Error: {e}')
        return None


def save_image(image, kegg_id):
    """Save a processed image with the KEGG ID as the filename."""
    try:
        image.save(os.path.join(OUTPUT_DIR, f'{kegg_id}.png'), format=IMAGE_FORMAT)
        return True
    except Exception as e:
        print(f'Save error: {e}')
        return False


# =============================================================================
# BATCH PROCESSING
# =============================================================================

def download_and_process_compounds(compounds, existing_ids, failed_ids):
    """
    Download and process structure images for a dict of KEGG compounds.
    Skips IDs that already have images or have previously failed.

    Args:
        compounds    (dict): { kegg_id: name }
        existing_ids (set):  IDs with existing PNG files
        failed_ids   (set):  IDs that previously failed
    """
    total               = len(compounds)
    skipped             = len(existing_ids & set(compounds.keys()))
    permanently_failed  = len(failed_ids  & set(compounds.keys()))
    to_download         = total - skipped - permanently_failed

    print(f'\n=== Processing {total} compounds ===')
    print(f'Already exist:     {skipped}')
    print(f'Previously failed: {permanently_failed}')
    print(f'To download:       {to_download}\n')

    mapping_cache      = load_mapping_cache()
    downloaded_count   = 0
    failed_count       = 0
    mapping_failed     = 0
    newly_failed       = []

    for i, (kegg_id, name) in enumerate(compounds.items(), 1):
        if kegg_id in existing_ids:
            print(f'[{i}/{total}] ⏭  Skipped {kegg_id} — already exists')
            continue

        if kegg_id in failed_ids:
            print(f'[{i}/{total}] ⏭  Skipped {kegg_id} — previously failed')
            continue

        print(f'[{i}/{total}] {kegg_id} ({name})... ', end='')

        cid = kegg_to_pubchem_cid(kegg_id, mapping_cache)
        if cid is None:
            print('✗ No PubChem mapping')
            mapping_failed += 1
            newly_failed.append(kegg_id)
            continue

        print(f'CID:{cid}... ', end='')

        image = download_pubchem_image(cid)
        if image is None:
            print('✗ Download failed')
            failed_count += 1
            newly_failed.append(kegg_id)
            continue

        try:
            processed = make_background_transparent(image)
            processed = crop_to_content(processed)
            if save_image(processed, kegg_id):
                print('✓')
                downloaded_count += 1
            else:
                print('✗ Save failed')
                failed_count += 1
                newly_failed.append(kegg_id)
        except Exception as e:
            print(f'✗ Error: {e}')
            failed_count += 1
            newly_failed.append(kegg_id)

    # Persist caches
    save_mapping_cache(mapping_cache)
    save_failed_downloads(list(set(list(failed_ids) + newly_failed)))

    print(f'\n=== Summary ===')
    print(f'Downloaded:          {downloaded_count}')
    print(f'Skipped (existing):  {skipped}')
    print(f'Failed (no mapping): {mapping_failed}')
    print(f'Failed (download):   {failed_count}')
    print(f'Previously failed:   {permanently_failed}')

    if downloaded_count > 0:
        generate_image_dimensions_manifest()


def reprocess_existing_images():
    """Batch reprocess all existing PNGs to crop transparent padding."""
    if not os.path.exists(OUTPUT_DIR):
        print(f'Output directory not found: {OUTPUT_DIR}')
        return

    png_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.png')]
    total     = len(png_files)
    print(f'\n=== Reprocessing {total} existing images ===')

    cropped = 0
    for i, filename in enumerate(sorted(png_files), 1):
        filepath = os.path.join(OUTPUT_DIR, filename)
        try:
            img           = Image.open(filepath).convert('RGBA')
            original_size = img.size
            result        = crop_to_content(img)
            if result.size != original_size:
                result.save(filepath, format=IMAGE_FORMAT)
                cropped += 1
            if i % 100 == 0 or i == total:
                print(f'  [{i}/{total}] processed ({cropped} cropped so far)')
        except Exception as e:
            print(f'  Error processing {filename}: {e}')

    print(f'\n=== Done: {cropped}/{total} images cropped ===')
    generate_image_dimensions_manifest()


def generate_image_dimensions_manifest():
    """
    Write a JSON manifest mapping KEGG IDs to { w, h } image dimensions.
    Used by the frontend to scale structure images proportionally.
    """
    if not os.path.exists(OUTPUT_DIR):
        print(f'Output directory not found: {OUTPUT_DIR}')
        return

    dimensions = {}
    for filename in sorted(os.listdir(OUTPUT_DIR)):
        if not filename.endswith('.png'):
            continue
        filepath = os.path.join(OUTPUT_DIR, filename)
        try:
            img     = Image.open(filepath)
            kegg_id = filename.replace('.png', '')
            dimensions[kegg_id] = {'w': img.width, 'h': img.height}
        except Exception as e:
            print(f'  Warning: Could not read {filename}: {e}')

    try:
        with open(IMAGE_DIMENSIONS_FILE, 'w') as f:
            json.dump(dimensions, f)
        print(f'Image dimensions manifest: {len(dimensions)} entries → {IMAGE_DIMENSIONS_FILE}')
    except IOError as e:
        print(f'Warning: Could not save image dimensions manifest: {e}')

    return dimensions


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def download_structures(json_file_path=None):
    """
    Download PubChem structure images for all KEGG compounds in the pathway JSON.

    Args:
        json_file_path (str, optional): Path to the Escher output JSON file.
            If None, the most recently modified *_output.json in
            ./static/json_pathway/ is used automatically.
    """
    print('=== PubChem Structure Image Downloader ===')
    print('Using PubChem (public domain images)\n')

    try:
        setup_output_directory()

        # --- Resolve which JSON file to read compounds from ---
        if json_file_path and os.path.exists(json_file_path):
            resolved_path = json_file_path
            print(f'Using provided JSON: {resolved_path}')
        else:
            resolved_path = _find_latest_output_json('./static/json_pathway')
            if resolved_path:
                print(f'Auto-detected JSON: {resolved_path}')

        if not resolved_path:
            print('No output JSON found — skipping structure image download')
            return

        failed_tracking = load_failed_downloads()
        failed_ids      = set(failed_tracking.get('failed', []))
        print(f'Loaded {len(failed_ids)} previously failed IDs')

        pathway_data = load_pathway_data(resolved_path)
        compounds    = extract_kegg_compounds(pathway_data)

        if not compounds:
            print('No KEGG compounds found in pathway data')
            return

        existing_ids = get_existing_images_cache()
        download_and_process_compounds(compounds, existing_ids, failed_ids)

        print('\n✓ Completed!')
        print('Note: PubChem images are public domain (US Government work).')
        print('Citation: Kim S, et al. PubChem 2023 update. Nucleic Acids Res. 2023')

    except Exception as e:
        print(f'\n✗ Error: {e}')
        raise


# =============================================================================
# CLI
# =============================================================================

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == '--reprocess':
        reprocess_existing_images()
    else:
        download_structures()