/**
 * SimplifiedPathwayApp - Minimal JavaScript for Flask forms integration
 * Handles sidebar toggle, node population, and Escher visualization initialization
 */
class SimplifiedPathwayApp {
    constructor() {
        this.escherBuilder = null;
        this.availableNodes = [];
    }

    initialize() {
        // Populate node selections
        this.populateNodeSelections();
        
        // Setup sidebar toggle
        this.setupSidebarToggle();
        
        // Setup config form event propagation
        this.setupConfigFormEventHandling();
        
        // Setup multi-node selector
        this.setupMultiNodeSelector();
        
        // Load and initialize Escher map
        if (window.initialJsonData) {
            this.initializeEscher();
        } else {
            console.error('No JSON data available for visualization');
        }
    }

    setupSidebarToggle() {
        const toggleBtn = document.getElementById('sidebar-toggle');
        const sidebar = document.getElementById('sidebar');
        
        if (toggleBtn && sidebar) {
            toggleBtn.addEventListener('click', () => {
                sidebar.classList.toggle('open');
            });
        }
    }

    populateNodeSelections() {
        fetch('/api/nodes')
            .then(response => {
                return response.json();
            })
            .then(nodes => {
                this.availableNodes = nodes;
                this.updateNodeSelects(nodes);
                this.updateMultiNodeSelect(nodes);
                this.restoreSelectionFromUrl();
            })
            .catch(error => {
                console.error('Error loading nodes:', error);
            });
    }

    /**
     * Restore start/end node selections and keep_positions checkbox from URL params.
     * Called after node selects are populated so the options exist to be selected.
     */
    restoreSelectionFromUrl() {
        const urlParams = new URLSearchParams(window.location.search);
        const urlStart = urlParams.get('start');
        const urlEnd = urlParams.get('end');
        const urlKeepPos = urlParams.get('keep_pos');

        if (urlStart) {
            const startSelect = document.getElementById('start-node');
            if (startSelect) {
                startSelect.value = urlStart;
                console.log('[RESTORE] Set start-node to:', urlStart);
            }
        }
        if (urlEnd) {
            const endSelect = document.getElementById('end-node');
            if (endSelect) {
                endSelect.value = urlEnd;
                console.log('[RESTORE] Set end-node to:', urlEnd);
            }
        }
        if (urlKeepPos !== null) {
            // Update all keep_positions checkboxes on the page
            const checkboxes = document.querySelectorAll('input[name="keep_positions"]');
            const checked = urlKeepPos === '1';
            checkboxes.forEach(cb => {
                cb.checked = checked;
            });
            console.log('[RESTORE] Set keep_positions to:', checked);
        }

        // Restore multi-node subgraph selections
        const urlSelected = urlParams.get('selected');
        const urlDist = urlParams.get('dist');

        if (urlSelected) {
            const selectedNodeIds = urlSelected.split(',').filter(n => n.trim());
            const multiSelect = document.getElementById('node-multiselect');
            const hiddenInput = document.getElementById('selected_nodes');

            if (multiSelect) {
                // Pre-select the options in the multi-select
                Array.from(multiSelect.options).forEach(opt => {
                    if (selectedNodeIds.includes(opt.value)) {
                        opt.selected = true;
                    }
                });
                console.log('[RESTORE] Set multi-node selection to:', selectedNodeIds);
            }
            if (hiddenInput) {
                hiddenInput.value = selectedNodeIds.join(',');
            }
        }
        if (urlDist) {
            const distInput = document.getElementById('connection_distance');
            if (distInput) {
                distInput.value = urlDist;
                console.log('[RESTORE] Set connection_distance to:', urlDist);
            }
        }
    }

    updateNodeSelects(nodes) {
        const startSelect = document.getElementById('start-node');
        const endSelect = document.getElementById('end-node');

        if (!startSelect || !endSelect) {
            return;
        }

        // Clear existing options (keep placeholder)
        startSelect.innerHTML = '<option value="">-- Select start node --</option>';
        endSelect.innerHTML = '<option value="">-- Select end node --</option>';

        // Add node options
        nodes.forEach((node, index) => {
            // Build display text: "Name (ID)" or just ID if no name
            const name = node.name || '';
            const displayText = (name && name !== node.id) ? `${name} (${node.id})` : node.id;
            
            const option1 = document.createElement('option');
            option1.value = node.id;
            option1.textContent = displayText;
            startSelect.appendChild(option1);

            const option2 = document.createElement('option');
            option2.value = node.id;
            option2.textContent = displayText;
            endSelect.appendChild(option2);
        });
    }

    updateMultiNodeSelect(nodes) {
        const multiSelect = document.getElementById('node-multiselect');
        if (!multiSelect) return;

        // Clear existing options (keep placeholder)
        multiSelect.innerHTML = '<option value="" disabled>-- Select nodes --</option>';

        // Add node options
        nodes.forEach(node => {
            // Build display text: "Name (ID)" or just ID if no name
            const name = node.name || '';
            const displayText = (name && name !== node.id) ? `${name} (${node.id})` : node.id;
            
            const option = document.createElement('option');
            option.value = node.id;
            option.textContent = displayText;
            multiSelect.appendChild(option);
        });
    }

    setupMultiNodeSelector() {
        const multiSelect = document.getElementById('node-multiselect');
        const hiddenInput = document.getElementById('selected_nodes');
        const searchInput = document.getElementById('node-search-input');

        if (!multiSelect || !hiddenInput) return;

        // Sync hidden input with multi-select selection on change
        // This updates the form data that gets submitted to Flask
        const syncSelectedNodes = () => {
            const selectedValues = Array.from(multiSelect.selectedOptions).map(opt => opt.value);
            hiddenInput.value = selectedValues.join(',');
        };
        
        multiSelect.addEventListener('change', syncSelectedNodes);

        // Setup typeahead search filtering if search input exists
        if (searchInput) {
            searchInput.addEventListener('input', (e) => {
                const searchTerm = e.target.value.toLowerCase();
                const options = multiSelect.querySelectorAll('option');
                
                options.forEach(option => {
                    if (option.value === '') {
                        // Always hide the placeholder
                        option.style.display = 'none';
                        return;
                    }
                    
                    const text = option.textContent.toLowerCase();
                    const matches = text.includes(searchTerm) || option.value.toLowerCase().includes(searchTerm);
                    option.style.display = (matches && searchTerm.length > 0) || searchTerm.length === 0 ? '' : 'none';
                });
            });
        }
    }

    setupConfigFormEventHandling() {
        // Block keyboard events at sidebar level to prevent Escher interference
        const sidebar = document.getElementById('sidebar');
        if (sidebar) {
            // Single consolidated listener on container level prevents all nested events from bubbling
            ['keydown', 'keyup', 'keypress'].forEach(eventType => {
                sidebar.addEventListener(eventType, (e) => {
                    e.stopPropagation();
                    e.stopImmediatePropagation();
                }, true); // Capture phase
            });
        }

        // Log width calculations from form inputs when they change
        const inputFields = document.querySelectorAll('.config-input');
        inputFields.forEach((field) => {
            field.addEventListener('change', () => {
                const smallWidth = document.querySelector('input[name="small_width"]')?.value || 'N/A';
                const smallHeight = document.querySelector('input[name="small_height"]')?.value || 'N/A';
                const mediumWidth = document.querySelector('input[name="medium_width"]')?.value || 'N/A';
                const mediumHeight = document.querySelector('input[name="medium_height"]')?.value || 'N/A';
                const largeWidth = document.querySelector('input[name="large_width"]')?.value || 'N/A';
                const largeHeight = document.querySelector('input[name="large_height"]')?.value || 'N/A';
                
                console.log('Canvas Dimensions Updated:', {
                    small: `${smallWidth}x${smallHeight}`,
                    medium: `${mediumWidth}x${mediumHeight}`,
                    large: `${largeWidth}x${largeHeight}`
                });
                
                // Also log node radius configuration
                console.log('Node Radius Configuration:', {
                    default: window.visualizationConfig.nodeRadius,
                    metabolite: window.visualizationConfig.metaboliteRadius,
                    reaction: window.visualizationConfig.reactionRadius
                });
            });
        });
    }

    initializeEscher() {
        const options = {
            enable_editing: true,
            fill_screen: false,
            reaction_styles: ["abs", "color", "size", "text"],
            never_ask_before_quit: true,
            primary_metabolite_radius: 100,
            first_load_callback: () => this.onEscherLoaded()
        };

        this.escherBuilder = escher.Builder(
            window.initialJsonData,
            null,
            null,
            d3.select("#map_container"),
            options
        );
    }

    onEscherLoaded() {
        // Initialize visualization overlays
        EscherVisualizer.initializeStructures(window.initialJsonData);
        EscherVisualizer.setupNaNLabelRemoval();
        
        // Setup export button
        this.setupExportButton();
        
        // Update subgraph status display
        this.updateSubgraphStatus();
    }

    setupExportButton() {
        const exportBtn = document.getElementById('export-svg-btn');
        if (exportBtn) {
            exportBtn.addEventListener('click', () => this.exportMapAsSVG());
        }
    }

    exportMapAsSVG() {
        try {
            // Get the SVG element from the map container
            const mapContainer = document.getElementById('map_container');
            const svgElement = mapContainer.querySelector('svg');
            
            if (!svgElement) {
                alert('No map SVG found to export');
                return;
            }
            
            // Clone the SVG FIRST before modifying anything
            const clonedSvg = svgElement.cloneNode(true);
            
            // Ensure proper namespaces are declared
            clonedSvg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
            clonedSvg.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
            
            // NOTE: Don't remove elements from the clone yet — the computed style copy
            // uses index-based matching between original and clone, so they must stay in sync.
            
            // Get all images from the CLONED SVG
            const clonedImages = clonedSvg.querySelectorAll('image');
            
            const imageConversionPromises = Array.from(clonedImages).map((imgElement, index) => {
                const href = imgElement.getAttribute('xlink:href') || imgElement.getAttribute('href');
                
                if (!href || href.startsWith('data:')) {
                    // Already a data URI or no href
                    return Promise.resolve();
                }
                
                return new Promise((resolve) => {
                    // Create a timeout promise to prevent hanging
                    const timeoutPromise = new Promise((_, reject) => {
                        setTimeout(() => reject(new Error('Image load timeout')), 5000);
                    });
                    
                    const img = new Image();
                    
                    const handleImageLoad = () => {
                        try {
                            // Convert loaded image to canvas then to data URI
                            const canvas = document.createElement('canvas');
                            canvas.width = img.width || 200;
                            canvas.height = img.height || 200;
                            const ctx = canvas.getContext('2d');
                            ctx.drawImage(img, 0, 0);
                            const dataUri = canvas.toDataURL('image/png');
                            
                            // Update the CLONED image element
                            imgElement.setAttribute('xlink:href', dataUri);
                            imgElement.removeAttribute('href'); // Remove conflicting attribute
                            resolve();
                        } catch (e) {
                            console.warn(`Error converting image ${index}:`, e);
                            resolve(); // Continue without this image
                        }
                    };
                    
                    const handleImageError = () => {
                        // Try fetch as fallback
                        fetch(href, { mode: 'cors', credentials: 'same-origin' })
                            .then(response => response.blob())
                            .then(blob => {
                                const reader = new FileReader();
                                reader.onload = () => {
                                    imgElement.setAttribute('xlink:href', reader.result);
                                    imgElement.removeAttribute('href');
                                    resolve();
                                };
                                reader.readAsDataURL(blob);
                            })
                            .catch(error => {
                                console.warn(`Could not load image ${index}:`, error);
                                resolve(); // Continue without this image
                            });
                    };
                    
                    img.onload = handleImageLoad;
                    img.onerror = handleImageError;
                    img.src = href;
                });
            });
            
            // Wait for all images to be converted
            Promise.all(imageConversionPromises)
                .then(() => {
                    // Apply computed styles from the ORIGINAL DOM elements to the cloned elements
                    const originalElements = svgElement.querySelectorAll('*');
                    const clonedElements = clonedSvg.querySelectorAll('*');
                    
                    originalElements.forEach((origElement, index) => {
                        const clonedElement = clonedElements[index];
                        if (!clonedElement) return;
                        
                        const computed = window.getComputedStyle(origElement);
                        const stylesToCopy = [
                            'display', 'visibility',
                            'fill', 'stroke', 'stroke-width', 'stroke-dasharray', 'stroke-linecap',
                            'stroke-linejoin', 'opacity', 'fill-opacity', 'stroke-opacity',
                            'font-size', 'font-family', 'font-weight', 'text-anchor',
                            'dominant-baseline', 'font-style', 'filter', 'marker-end', 'marker-start'
                        ];
                        
                        let styleStr = '';
                        stylesToCopy.forEach(prop => {
                            const value = computed.getPropertyValue(prop);
                            if (value && value !== 'initial') {
                                // Keep 'none' for display/visibility (hidden elements must stay hidden)
                                // but skip 'none' for other properties like fill, stroke, etc.
                                if (value === 'none' && prop !== 'display' && prop !== 'visibility') return;
                                styleStr += `${prop}: ${value}; `;
                            }
                        });
                        
                        if (styleStr) {
                            clonedElement.setAttribute('style', styleStr);
                        }
                    });
                    
                    // Explicitly hide Escher's default bigg_id labels in the cloned SVG
                    clonedSvg.querySelectorAll('.node-label.label').forEach(el => {
                        el.setAttribute('style', 'display: none;');
                    });
                    
                    // Remove Escher UI elements that may be inside the SVG
                    clonedSvg.querySelectorAll('.menu-bar, .button-panel, .search-bar-container, .dropdown-menu').forEach(el => el.remove());
                    // Remove any foreignObject elements (Escher uses these for HTML UI overlays)
                    clonedSvg.querySelectorAll('foreignObject').forEach(el => el.remove());
                    
                    // Add CSS styles as fallback
                    const styleSheet = document.createElement('style');
                    styleSheet.innerHTML = `
                        text { font-family: Arial, sans-serif; }
                        .stoichiometry-label { font-size: 12px; }
                        .node-label.label { display: none !important; }
                        .node-label.metabolite-name, .node-label.coproduct-name { font-size: 12px; font-weight: bold; }
                        .segment { stroke: black; stroke-width: 3px; stroke-dasharray: 5,5; }
                        path.segment { stroke: black; stroke-width: 3px; stroke-dasharray: 5,5; }
                        .node-circle { fill-opacity: 0.8; }
                        .reaction-circle { fill-opacity: 0.8; }
                        circle { fill-opacity: 0.8; }
                        .menu-bar, .button-panel, .search-bar-container { display: none !important; }
                    `;
                    clonedSvg.insertBefore(styleSheet, clonedSvg.firstChild);
                    
                    // Convert SVG to string
                    const svgString = new XMLSerializer().serializeToString(clonedSvg);
                    const imageDataCount = (svgString.match(/data:image/g) || []).length;
                    
                    // Log SVG export details
                    console.log('SVG Export:', {
                        svgStringLength: svgString.length,
                        embeddedImages: imageDataCount,
                        totalImages: clonedImages.length
                    });
                    
                    // Build descriptive filename
                    const baseName = this._buildExportFilename();
                    
                    // --- SVG download ---
                    const blob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
                    const url = URL.createObjectURL(blob);
                    const link = document.createElement('a');
                    const svgFilename = baseName + '.svg';
                    
                    link.href = url;
                    link.download = svgFilename;
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                    URL.revokeObjectURL(url);
                    
                    // --- PNG download ---
                    this._exportSvgAsPng(svgString, clonedSvg, baseName + '.png');
                    
                    alert('Map exported as SVG + PNG: ' + baseName + '\n(' + imageDataCount + ' of ' + clonedImages.length + ' images embedded)');
                })
                .catch(error => {
                    console.error('Error during export:', error);
                    alert('Error exporting map: ' + error.message);
                });
        } catch (error) {
            console.error('Error exporting map as SVG:', error);
            alert('Error exporting map: ' + error.message);
        }
    }

    /**
     * Build a descriptive filename based on current view context.
     * Format: pathway_{start}-{end}_{orientation}_{date}  (shortest path)
     *         pathway_multinode_{orientation}_{date}       (multi-node)
     *         pathway_full_{orientation}_{date}            (full graph)
     */
    /**
     * Look up the descriptive name for a node ID (e.g. "C00022" → "Pyruvate").
     * Falls back to the raw ID if no name is found.
     */
    _nodeName(nodeId) {
        const node = this.availableNodes.find(n => n.id === nodeId);
        if (node && node.name && node.name !== nodeId) {
            return node.name;
        }
        return nodeId;
    }

    /**
     * Sanitize a string for use in a filename (letters, digits, hyphens only).
     */
    _sanitize(str) {
        return str.replace(/[^a-zA-Z0-9-]/g, '_').replace(/_+/g, '_').replace(/^_|_$/g, '');
    }

    _buildExportFilename() {
        const urlParams = new URLSearchParams(window.location.search);
        const isVertical = window.visualizationConfig?.smallGraphLayoutVertical ?? false;
        const orientation = isVertical ? 'vertical' : 'horizontal';
        const timestamp = new Date().toISOString().slice(0, 10);

        const isSubgraph = urlParams.get('view') === 'subgraph';
        const startNode = urlParams.get('start') || '';
        const endNode = urlParams.get('end') || '';
        const selectedNodes = urlParams.get('selected') || '';

        if (isSubgraph && startNode && endNode) {
            // Shortest path — use descriptive names
            const s = this._sanitize(this._nodeName(startNode));
            const e = this._sanitize(this._nodeName(endNode));
            return `pathway_${s}-${e}_${orientation}_${timestamp}`;
        } else if (isSubgraph && selectedNodes) {
            // Multi-node subgraph — use descriptive names
            const nodeList = selectedNodes.split(',').map(n => this._sanitize(this._nodeName(n.trim())));
            const nodeStr = nodeList.length <= 3 ? nodeList.join('-') : nodeList.slice(0, 3).join('-') + `_plus${nodeList.length - 3}`;
            return `pathway_multi_${nodeStr}_${orientation}_${timestamp}`;
        } else if (isSubgraph) {
            return `pathway_subgraph_${orientation}_${timestamp}`;
        } else {
            return `pathway_full_${orientation}_${timestamp}`;
        }
    }

    /**
     * Render an SVG string to a canvas and download as PNG.
     */
    _exportSvgAsPng(svgString, clonedSvg, filename) {
        // Determine dimensions from the SVG viewBox or width/height
        const vb = clonedSvg.getAttribute('viewBox');
        let width, height;
        if (vb) {
            const parts = vb.split(/[\s,]+/).map(Number);
            width = parts[2];
            height = parts[3];
        } else {
            width = parseFloat(clonedSvg.getAttribute('width')) || 1920;
            height = parseFloat(clonedSvg.getAttribute('height')) || 1080;
        }

        // Scale up for high-res output (2x)
        const scale = 2;
        const canvas = document.createElement('canvas');
        canvas.width = width * scale;
        canvas.height = height * scale;
        const ctx = canvas.getContext('2d');
        ctx.scale(scale, scale);

        const img = new Image();
        const svgBlob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
        const svgUrl = URL.createObjectURL(svgBlob);

        img.onload = () => {
            // White background
            ctx.fillStyle = '#ffffff';
            ctx.fillRect(0, 0, width, height);
            ctx.drawImage(img, 0, 0, width, height);
            URL.revokeObjectURL(svgUrl);

            canvas.toBlob((blob) => {
                if (!blob) {
                    console.error('PNG export: canvas.toBlob returned null');
                    return;
                }
                const pngUrl = URL.createObjectURL(blob);
                const link = document.createElement('a');
                link.href = pngUrl;
                link.download = filename;
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
                URL.revokeObjectURL(pngUrl);
                console.log('PNG exported:', filename, `(${canvas.width}x${canvas.height})`);
            }, 'image/png');
        };

        img.onerror = (e) => {
            console.error('PNG export: failed to load SVG into image', e);
            URL.revokeObjectURL(svgUrl);
        };

        img.src = svgUrl;
    }

    rebuildEscher(newJsonData) {
        console.log('[REBUILD] Rebuilding Escher with new data...');
        
        // Update global data
        window.initialJsonData = newJsonData;
        
        // Clear the map container completely
        const mapContainer = document.getElementById('map_container');
        if (mapContainer) {
            mapContainer.innerHTML = '';
        }
        
        // Clear any D3 overlays that might be outside the container
        if (typeof d3 !== 'undefined') {
            d3.selectAll('.bar-chart').remove();
            d3.selectAll('.metabolite-name').remove();
            d3.selectAll('.coproduct-name').remove();
        }
        
        // Destroy old builder reference
        this.escherBuilder = null;
        
        // Re-initialize Escher (which triggers onEscherLoaded -> initializeStructures)
        this.initializeEscher();
        
        // Update subgraph status display
        this.updateSubgraphStatus();
        
        console.log('[REBUILD] Escher rebuild complete');
    }

    updateSubgraphStatus() {
        const statusContainer = document.getElementById('subgraph-status-container');
        const urlParams = new URLSearchParams(window.location.search);
        const isSubgraph = urlParams.get('view') === 'subgraph';
        
        if (statusContainer) {
            statusContainer.style.display = isSubgraph ? 'block' : 'none';
        }
    }
}

// Initialize app when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    const app = new SimplifiedPathwayApp();
    app.initialize();
    // Expose globally so other scripts (e.g. regenerate handler) can call rebuildEscher
    window.pathwayApp = app;
});
