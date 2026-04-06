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

# ===== CONFIGURATION CONSTANTS =====
# JSON file path - will be set dynamically based on uploaded file
JSON_FILE_PATH = './static/json_pathway/metabolite_graph_output.json'  # Default; may be overridden
OUTPUT_DIR = './static/structure_imgs'
FAILED_DOWNLOADS_FILE = './static/failed_downloads.json'
MAPPING_CACHE_FILE = './static/kegg_pubchem_mapping.json'

# PubChem API endpoints
PUBCHEM_CID_LOOKUP = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON'
PUBCHEM_IMAGE_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/PNG?image_size=500x500'
KEGG_TO_PUBCHEM_API = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{}/cids/JSON'

TRANSPARENT_COLORS = [(255, 255, 255), (245, 245, 245)]
IMAGE_FORMAT = 'PNG'
MAX_RETRIES = 3
RETRY_DELAY = 2  # PubChem recommends no more than 5 requests per second
REQUEST_DELAY = 0.3  # Delay between requests to respect rate limits

# ===== UTILITY FUNCTIONS =====
def setup_output_directory():
    """Create output directory if it doesn't exist."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(MAPPING_CACHE_FILE), exist_ok=True)

def load_mapping_cache():
    """Load cached KEGG to PubChem CID mappings."""
    if not os.path.exists(MAPPING_CACHE_FILE):
        return {}
    
    try:
        with open(MAPPING_CACHE_FILE, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        print("Warning: Could not read mapping cache, starting fresh")
        return {}

def save_mapping_cache(mapping):
    """Save KEGG to PubChem CID mappings."""
    try:
        with open(MAPPING_CACHE_FILE, 'w') as f:
            json.dump(mapping, f, indent=2)
    except IOError as e:
        print(f"Warning: Could not save mapping cache: {e}")

def load_failed_downloads():
    """Load previously failed downloads from tracking file."""
    if not os.path.exists(FAILED_DOWNLOADS_FILE):
        return {'failed': [], 'last_updated': None}
    
    try:
        with open(FAILED_DOWNLOADS_FILE, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError):
        print("Warning: Could not read failed downloads file")
        return {'failed': [], 'last_updated': None}

def save_failed_downloads(failed_list):
    """Save list of failed downloads to tracking file."""
    data = {
        'failed': failed_list,
        'last_updated': time.strftime('%Y-%m-%d %H:%M:%S')
    }
    try:
        with open(FAILED_DOWNLOADS_FILE, 'w') as f:
            json.dump(data, f, indent=2)
    except IOError as e:
        print(f"Warning: Could not save failed downloads file: {e}")

def get_existing_images_cache():
    """Build a cached set of existing image filenames."""
    if not os.path.exists(OUTPUT_DIR):
        return set()
    
    try:
        existing_files = os.listdir(OUTPUT_DIR)
        return {f.replace('.png', '') for f in existing_files if f.endswith('.png')}
    except Exception as e:
        print(f"Warning: Could not read existing images directory: {e}")
        return set()

def load_pathway_data(json_path):
    """Load pathway data from JSON file."""
    try:
        with open(json_path, 'r') as file:
            return json.load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"JSON file not found: {json_path}")
    except json.JSONDecodeError as e:
        raise json.JSONDecodeError(f"Invalid JSON format: {e}", "", 0)

def extract_kegg_compounds(pathway_data):
    """Extract KEGG compound information from pathway data."""
    compounds = {}
    
    if len(pathway_data) > 1 and 'nodes' in pathway_data[1]:
        nodes = pathway_data[1]['nodes']
        
        for node_id, node_info in nodes.items():
            bigg_id = node_info.get('bigg_id')
            name = node_info.get('name', 'Unknown')
            
            if bigg_id and bigg_id.startswith('C') and len(bigg_id) >= 5:
                compounds[bigg_id] = name.replace(';', '').strip()
    
    return compounds

def kegg_to_pubchem_cid(kegg_id, mapping_cache):
    """
    Convert KEGG ID to PubChem CID using multiple methods.
    
    Args:
        kegg_id (str): KEGG compound ID (e.g., 'C00031')
        mapping_cache (dict): Cache of previously resolved mappings
        
    Returns:
        str or None: PubChem CID if found
    """
    # Check cache first
    if kegg_id in mapping_cache:
        return mapping_cache[kegg_id]
    
    time.sleep(REQUEST_DELAY)  # Rate limiting
    
    # Method 1: Try KEGG cross-reference lookup
    try:
        url = KEGG_TO_PUBCHEM_API.format(kegg_id)
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                cid = str(data['IdentifierList']['CID'][0])
                mapping_cache[kegg_id] = cid
                return cid
    except Exception as e:
        pass  # Continue to next method
    
    # Method 2: Search by KEGG ID as synonym
    time.sleep(REQUEST_DELAY)
    try:
        url = PUBCHEM_CID_LOOKUP.format(kegg_id)
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                cid = str(data['IdentifierList']['CID'][0])
                mapping_cache[kegg_id] = cid
                return cid
    except Exception as e:
        pass
    
    # Mapping failed
    mapping_cache[kegg_id] = None
    return None

def make_background_transparent(image):
    """Make white and off-white backgrounds transparent."""
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
    """Crop transparent padding around the actual image content.
    
    Uses the alpha channel bounding box to find where the visible
    content is and crops away surrounding transparent pixels.
    
    Args:
        image (PIL.Image): RGBA image to crop
        padding (int): Pixels of padding to keep around the content
        
    Returns:
        PIL.Image: Cropped image with minimal transparent border
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    # Get bounding box of non-transparent pixels (alpha > 0)
    alpha = image.split()[3]  # Extract alpha channel
    bbox = alpha.getbbox()
    
    if bbox is None:
        # Entire image is transparent, return as-is
        return image
    
    # Add small padding around the content
    left = max(0, bbox[0] - padding)
    top = max(0, bbox[1] - padding)
    right = min(image.width, bbox[2] + padding)
    bottom = min(image.height, bbox[3] + padding)
    
    return image.crop((left, top, right, bottom))

def reprocess_existing_images():
    """Batch reprocess all existing structure images to crop transparent padding."""
    if not os.path.exists(OUTPUT_DIR):
        print(f"Output directory not found: {OUTPUT_DIR}")
        return
    
    png_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.png')]
    total = len(png_files)
    print(f"\n=== Reprocessing {total} existing images ===")
    
    cropped_count = 0
    for i, filename in enumerate(sorted(png_files), 1):
        filepath = os.path.join(OUTPUT_DIR, filename)
        try:
            img = Image.open(filepath).convert('RGBA')
            original_size = img.size
            cropped = crop_to_content(img)
            if cropped.size != original_size:
                cropped.save(filepath, format=IMAGE_FORMAT)
                cropped_count += 1
            if i % 100 == 0 or i == total:
                print(f"  [{i}/{total}] processed ({cropped_count} cropped so far)")
        except Exception as e:
            print(f"  Error processing {filename}: {e}")
    
    print(f"\n=== Done: {cropped_count}/{total} images were cropped ===")
    # Regenerate dimensions manifest after reprocessing
    generate_image_dimensions_manifest()

IMAGE_DIMENSIONS_FILE = os.path.join(OUTPUT_DIR, 'image_dimensions.json')

def generate_image_dimensions_manifest():
    """Generate a JSON manifest mapping KEGG IDs to image {width, height}.
    
    This allows the frontend to scale structure images proportionally,
    so small molecules (like water) appear smaller than complex ones.
    """
    if not os.path.exists(OUTPUT_DIR):
        print(f"Output directory not found: {OUTPUT_DIR}")
        return
    
    dimensions = {}
    png_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.png')]
    
    for filename in sorted(png_files):
        filepath = os.path.join(OUTPUT_DIR, filename)
        try:
            img = Image.open(filepath)
            kegg_id = filename.replace('.png', '')
            dimensions[kegg_id] = {'w': img.width, 'h': img.height}
        except Exception as e:
            print(f"  Warning: Could not read {filename}: {e}")
    
    try:
        with open(IMAGE_DIMENSIONS_FILE, 'w') as f:
            json.dump(dimensions, f)
        print(f"Image dimensions manifest: {len(dimensions)} entries → {IMAGE_DIMENSIONS_FILE}")
    except IOError as e:
        print(f"Warning: Could not save image dimensions manifest: {e}")
    
    return dimensions

def download_pubchem_image(cid, attempt=1):
    """
    Download structure image from PubChem with retry logic.
    
    Args:
        cid (str): PubChem Compound ID
        attempt (int): Current attempt number
        
    Returns:
        PIL.Image or None: Downloaded image or None if failed
    """
    url = PUBCHEM_IMAGE_URL.format(cid)
    
    try:
        time.sleep(REQUEST_DELAY)  # Rate limiting
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return Image.open(BytesIO(response.content))
    except requests.Timeout:
        if attempt < MAX_RETRIES:
            print(f"⏳ Timeout (attempt {attempt}/{MAX_RETRIES}), retrying... ", end="")
            time.sleep(RETRY_DELAY)
            return download_pubchem_image(cid, attempt + 1)
        else:
            return None
    except requests.RequestException:
        if attempt < MAX_RETRIES:
            time.sleep(RETRY_DELAY)
            return download_pubchem_image(cid, attempt + 1)
        else:
            return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def save_image(image, kegg_id):
    """Save processed image with KEGG ID as filename."""
    try:
        image_path = os.path.join(OUTPUT_DIR, f'{kegg_id}.png')
        image.save(image_path, format=IMAGE_FORMAT)
        return True
    except Exception as e:
        print(f"Save error: {e}")
        return False

def download_and_process_compounds(compounds, existing_ids, failed_ids):
    """
    Download and process structure images for compounds.
    
    Args:
        compounds (dict): KEGG IDs to compound names
        existing_ids (set): IDs that already have images
        failed_ids (set): IDs that previously failed
    """
    total = len(compounds)
    skipped = len(existing_ids & set(compounds.keys()))
    permanently_failed = len(failed_ids & set(compounds.keys()))
    to_download = total - skipped - permanently_failed
    
    print(f"\n=== Processing {total} compounds ===")
    print(f"Already exist: {skipped}")
    print(f"Previously failed: {permanently_failed}")
    print(f"To download: {to_download}\n")
    
    # Load mapping cache
    mapping_cache = load_mapping_cache()
    
    downloaded_count = 0
    failed_count = 0
    mapping_failed_count = 0
    newly_failed = []
    
    for i, (kegg_id, name) in enumerate(compounds.items(), 1):
        # Skip if already exists
        if kegg_id in existing_ids:
            print(f"[{i}/{total}] ⏭ Skipped {kegg_id} - exists")
            continue
        
        # Skip if permanently failed
        if kegg_id in failed_ids:
            print(f"[{i}/{total}] ⏭ Skipped {kegg_id} - previously failed")
            continue
        
        print(f"[{i}/{total}] {kegg_id} ({name})... ", end="")
        
        # Map KEGG ID to PubChem CID
        cid = kegg_to_pubchem_cid(kegg_id, mapping_cache)
        if cid is None:
            print("✗ No PubChem mapping")
            mapping_failed_count += 1
            newly_failed.append(kegg_id)
            continue
        
        print(f"CID:{cid}... ", end="")
        
        # Download image
        image = download_pubchem_image(cid)
        if image is None:
            print("✗ Download failed")
            failed_count += 1
            newly_failed.append(kegg_id)
            continue
        
        # Process and save
        try:
            processed_image = make_background_transparent(image)
            processed_image = crop_to_content(processed_image)
            if save_image(processed_image, kegg_id):
                print("✓")
                downloaded_count += 1
            else:
                print("✗ Save failed")
                failed_count += 1
                newly_failed.append(kegg_id)
        except Exception as e:
            print(f"✗ Error: {e}")
            failed_count += 1
            newly_failed.append(kegg_id)
    
    # Save updated caches
    save_mapping_cache(mapping_cache)
    all_failed = list(set(list(failed_ids) + newly_failed))
    save_failed_downloads(all_failed)
    
    # Summary
    print(f"\n=== Summary ===")
    print(f"Downloaded: {downloaded_count}")
    print(f"Skipped (existing): {skipped}")
    print(f"Failed (no mapping): {mapping_failed_count}")
    print(f"Failed (download): {failed_count}")
    print(f"Previously failed: {permanently_failed}")
    print(f"Total processed: {downloaded_count + skipped}")
    
    # Regenerate dimensions manifest after downloading new images
    if downloaded_count > 0:
        generate_image_dimensions_manifest()

# ===== MAIN EXECUTION =====
def download_structures():
    """Main execution function."""
    print("=== PubChem Structure Image Downloader ===")
    print("Using PubChem (public domain images)")
    print("API Documentation: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest\n")
    
    try:
        setup_output_directory()
        
        # Load tracking data
        failed_tracking = load_failed_downloads()
        failed_ids = set(failed_tracking.get('failed', []))
        print(f"Loaded {len(failed_ids)} previously failed IDs\n")
        
        # Load pathway data
        pathway_data = load_pathway_data(JSON_FILE_PATH)
        compounds = extract_kegg_compounds(pathway_data)
        
        if not compounds:
            print("No KEGG compounds found to process")
            return
        
        # Get existing images
        existing_ids = get_existing_images_cache()
        
        # Download and process
        download_and_process_compounds(compounds, existing_ids, failed_ids)
        
        print("\n✓ Completed!")
        print("\nNote: PubChem images are public domain (US Government work)")
        print("No attribution required, but recommended citation:")
        print("Kim S, et al. PubChem 2023 update. Nucleic Acids Res. 2023")
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        raise

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == '--reprocess':
        reprocess_existing_images()
    else:
        download_structures()