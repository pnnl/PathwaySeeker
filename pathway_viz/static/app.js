// app.js
/**
 * PathwayApp
 *
 * Owns:
 *   this.initialJsonData  — Escher map JSON (from Flask via window.PathwayApp)
 *   this.config           — merged defaults + Flask config (single source of truth)
 *   this.visualizer       — EscherVisualizer instance
 *
 * window.PathwayApp is read exactly once at boot.
 * After that nothing reads window.* (except window.pathwayApp for console debugging).
 */

// =============================================================================
// CONFIG DEFAULTS  — single source of truth for fallback values
// =============================================================================
const CONFIG_DEFAULTS = Object.freeze({
    imageSize:                200,
    defaultToBiggId:          false,
    nodeRadius:               15,
    metaboliteRadius:         10,
    reactionRadius:           8,
    labelOffsetY:             20,
    coproductLabelOffsetY:    25,
    barChartOffsetY:          60,
    metaboliteLabelFontSize:  12,
    coproductLabelFontSize:   11,
    chartTitleFontSize:       10,
    chartLabelFontSize:       9,
    barChartXLabel:           '',
    barChartTitle:            '',
    barChartYLabel:           '',
    smallGraphLayoutVertical: true,
    nodeThresholdSmall:       10,
    barChartWidth:            180,
    barChartAxisPadding:      55,
    barHeight:                12,
    barChartGap:              18,
    originColours: {
        metabolomics: '#2a9d8f',
        proteomics:   '#e76f51',
        both:         '#9b5de5',
        unknown:      '#aaaaaa',
    },
    pathHighlight: {
        nodeStrokeWidth:    5,
        nodeStrokeColour:   '#f4d03f',
        segmentStrokeWidth: 6,
        segmentColour:      '#f4d03f',
        dimOpacity:         0.25,
    },
});

// =============================================================================
// UTILITIES
// =============================================================================
function showStatus(elementId, message, type = 'info', autoDismissMs = 3000) {
    const el = document.getElementById(elementId);
    if (!el) return;
    el.className    = `status-message status-${type}`;
    el.textContent  = message;
    el.style.display = 'block';
    if (autoDismissMs > 0) {
        setTimeout(() => { el.style.display = 'none'; }, autoDismissMs);
    }
}

function debounce(fn, delay) {
    let timer;
    return (...args) => {
        clearTimeout(timer);
        timer = setTimeout(() => fn(...args), delay);
    };
}

/**
 * Merge Flask config on top of defaults.
 * Handles nested objects (originColours, pathHighlight) safely.
 */
function buildConfig(rawConfig) {
    return {
        ...CONFIG_DEFAULTS,
        ...rawConfig,
        originColours: {
            ...CONFIG_DEFAULTS.originColours,
            ...(rawConfig.originColours || {}),
        },
        pathHighlight: {
            ...CONFIG_DEFAULTS.pathHighlight,
            ...(rawConfig.pathHighlight || {}),
        },
        // Normalise barChart sub-object from flat Flask keys
        barChart: {
            width:       rawConfig.barChartWidth       ?? CONFIG_DEFAULTS.barChartWidth,
            axisPadding: rawConfig.barChartAxisPadding ?? CONFIG_DEFAULTS.barChartAxisPadding,
            barHeight:   rawConfig.barHeight           ?? CONFIG_DEFAULTS.barHeight,
            gapBetween:  rawConfig.barChartGap         ?? CONFIG_DEFAULTS.barChartGap,
        },
    };
}

// =============================================================================
// APP
// =============================================================================
class PathwayApp {
    /**
     * @param {object|null} initialJsonData
     * @param {object}      rawConfig        - raw config from Flask
     * @param {string}      highlightPath    - comma-separated node IDs
     * @param {string}      viewType         - 'full' | 'subgraph'
     */
    constructor(initialJsonData, rawConfig, highlightPath, viewType) {
        this.initialJsonData = initialJsonData;
        this.config          = buildConfig(rawConfig || {});
        this.highlightPath   = highlightPath || '';
        this.viewType        = viewType || 'full';

        this.escherBuilder  = null;
        this.availableNodes = [];
        this.visualizer     = null;   // set once Escher loads

        this._configAbort   = null;   // AbortController for config AJAX
        this._els           = {};     // DOM element cache

        console.log('[PathwayApp] Created — viewType:', this.viewType,
            '| hasData:', !!initialJsonData,
            '| highlightPath:', this.highlightPath || '(none)');
    }

    // ── DOM cache ─────────────────────────────────────────────────────────
    _el(id) {
        if (!this._els[id]) this._els[id] = document.getElementById(id);
        return this._els[id];
    }

    // ─────────────────────────────────────────────────────────────────────
    initialize() {
        console.log('[PathwayApp] Initializing');
        this.setupSidebarToggle();
        this.setupKeyboardIsolation();
        this.populateNodeSelections();
        this.setupMultiNodeSelector();
        this.setupFrontendConfigAutoApply();
        this.setupBackendConfigViewState();

        if (this.initialJsonData) {
            this.initializeEscher();
        } else {
            console.warn('[PathwayApp] No JSON data — map will not render');
        }
    }

    // ─────────────────────────────────────────────────────────────────────
    // SIDEBAR
    // ─────────────────────────────────────────────────────────────────────
    setupSidebarToggle() {
        const btn  = this._el('sidebar-toggle');
        const side = this._el('sidebar');
        if (btn && side) btn.addEventListener('click', () => side.classList.toggle('open'));
    }

    setupKeyboardIsolation() {
        const sidebar = this._el('sidebar');
        if (!sidebar) return;
        ['keydown', 'keyup', 'keypress'].forEach(type => {
            sidebar.addEventListener(type, e => {
                e.stopPropagation();
                e.stopImmediatePropagation();
            }, true);
        });
    }

    // ─────────────────────────────────────────────────────────────────────
    // NODE SELECTIONS
    // ─────────────────────────────────────────────────────────────────────
    populateNodeSelections() {
        console.log('[PathwayApp] Fetching node list');
        fetch('/api/nodes')
            .then(r => {
                if (!r.ok) throw new Error(`HTTP ${r.status}`);
                return r.json();
            })
            .then(nodes => {
                this.availableNodes = nodes;
                console.log('[PathwayApp] Loaded', nodes.length, 'nodes');
                this._fillSelect('start-node', nodes, '-- Select start node --');
                this._fillSelect('end-node',   nodes, '-- Select end node --');
                this._fillMultiSelect(nodes);
                this.restoreSelectionFromUrl();
            })
            .catch(err => console.error('[PathwayApp] Error loading nodes:', err));
    }

    _fillSelect(id, nodes, placeholder) {
        const el = this._el(id);
        if (!el) return;
        el.innerHTML = `<option value="">${placeholder}</option>`;
        nodes.forEach(n => {
            const label = (n.name && n.name !== n.id) ? `${n.name} (${n.id})` : n.id;
            const opt   = document.createElement('option');
            opt.value       = n.id;
            opt.textContent = label;
            el.appendChild(opt);
        });
    }

    _fillMultiSelect(nodes) {
        const el = this._el('node-multiselect');
        if (!el) return;
        el.innerHTML = '<option value="" disabled>-- Select nodes --</option>';
        nodes.forEach(n => {
            const label = (n.name && n.name !== n.id) ? `${n.name} (${n.id})` : n.id;
            const opt   = document.createElement('option');
            opt.value       = n.id;
            opt.textContent = label;
            el.appendChild(opt);
        });
    }

    restoreSelectionFromUrl() {
        const p = new URLSearchParams(window.location.search);
        const setVal = (id, val) => {
            const el = this._el(id);
            if (el && val) el.value = val;
        };
        setVal('start-node', p.get('start'));
        setVal('end-node',   p.get('end'));

        const keepPos = p.get('keep_pos');
        if (keepPos !== null) {
            document.querySelectorAll('input[name="keep_positions"]')
                .forEach(cb => { cb.checked = keepPos === '1'; });
        }

        const selected = p.get('selected');
        if (selected) {
            const ids   = selected.split(',').filter(Boolean);
            const multi = this._el('node-multiselect');
            if (multi) {
                Array.from(multi.options).forEach(o => {
                    o.selected = ids.includes(o.value);
                });
            }
            const hidden = this._el('selected_nodes');
            if (hidden) hidden.value = ids.join(',');
        }

        const dist = p.get('dist');
        if (dist) setVal('connection_distance', dist);
    }

    setupMultiNodeSelector() {
        const multi  = this._el('node-multiselect');
        const hidden = this._el('selected_nodes');
        const search = this._el('node-search-input');
        if (!multi || !hidden) return;

        multi.addEventListener('change', () => {
            hidden.value = Array.from(multi.selectedOptions)
                .map(o => o.value).join(',');
        });

        if (search) {
            search.addEventListener('input', e => {
                const term = e.target.value.toLowerCase();
                Array.from(multi.options).forEach(opt => {
                    if (!opt.value) { opt.style.display = 'none'; return; }
                    const match = !term
                        || opt.textContent.toLowerCase().includes(term)
                        || opt.value.toLowerCase().includes(term);
                    opt.style.display = match ? '' : 'none';
                });
            });
        }
    }

    // ─────────────────────────────────────────────────────────────────────
    // BACKEND CONFIG
    // ─────────────────────────────────────────────────────────────────────
    setupBackendConfigViewState() {
        const form = this._el('backend-config-form');
        if (!form) return;

        form.querySelectorAll('.backend-config-input').forEach(input => {
            input.addEventListener('change', () => input.classList.add('backend-input-changed'));
        });

        form.addEventListener('submit', () => {
            const p          = new URLSearchParams(window.location.search);
            const isSubgraph = p.get('view') === 'subgraph';
            const set = (id, val) => {
                const el = document.getElementById(id);
                if (el) el.value = val;
            };
            set('hidden-view-type',           isSubgraph ? 'subgraph' : 'full');
            set('hidden-start-node',          p.get('start')    || '');
            set('hidden-end-node',            p.get('end')      || '');
            set('hidden-path-nodes',          p.get('nodes')    || '');
            set('hidden-selected-nodes',      p.get('selected') || '');
            set('hidden-connection-distance', p.get('dist')     || '');
            set('hidden-keep-positions',      p.get('keep_pos') || '1');
            console.log('[PathwayApp] Backend form submit — view:',
                isSubgraph ? 'subgraph' : 'full');
        });
    }

    // ─────────────────────────────────────────────────────────────────────
    // FRONTEND CONFIG  (debounced auto-apply)
    // ─────────────────────────────────────────────────────────────────────
    setupFrontendConfigAutoApply() {
        const inputs = document.querySelectorAll('.frontend-config-input');
        if (!inputs.length) return;

        const applyConfig = debounce(() => {
            const incoming  = this._collectFrontendConfig();
            const changedKeys = Object.keys(incoming).filter(
                k => incoming[k] !== this.config[k]
            );
            if (!changedKeys.length) return;
            this._sendFrontendConfig(incoming, changedKeys);
        }, 800);

        inputs.forEach(input => {
            input.addEventListener('input',  applyConfig);
            input.addEventListener('change', applyConfig);
        });
    }

    _collectFrontendConfig() {
        // Fixed parseInt: Number.isFinite handles 0 correctly
        const num = (id, fallback) => {
            const el = this._el(id);
            if (!el) return fallback;
            const parsed = parseInt(el.value, 10);
            return Number.isFinite(parsed) ? parsed : fallback;
        };
        const txt = (id, fallback = '') => {
            const el = this._el(id);
            return el ? el.value : fallback;
        };

        return {
            nodeRadius:              num('nodeRadius',              CONFIG_DEFAULTS.nodeRadius),
            metaboliteRadius:        num('metaboliteRadius',        CONFIG_DEFAULTS.metaboliteRadius),
            reactionRadius:          num('reactionRadius',          CONFIG_DEFAULTS.reactionRadius),
            imageSize:               num('imageSize',               CONFIG_DEFAULTS.imageSize),
            labelOffsetY:            num('labelOffsetY',            CONFIG_DEFAULTS.labelOffsetY),
            coproductLabelOffsetY:   num('coproductLabelOffsetY',   CONFIG_DEFAULTS.coproductLabelOffsetY),
            metaboliteLabelFontSize: num('metaboliteLabelFontSize', CONFIG_DEFAULTS.metaboliteLabelFontSize),
            coproductLabelFontSize:  num('coproductLabelFontSize',  CONFIG_DEFAULTS.coproductLabelFontSize),
            barChartWidth:           num('barChartWidth',           CONFIG_DEFAULTS.barChartWidth),
            barChartHeight:          num('barChartHeight',          CONFIG_DEFAULTS.barChartWidth),
            barHeight:               num('barHeight',               CONFIG_DEFAULTS.barHeight),
            barChartOffsetY:         num('barChartOffsetY',         CONFIG_DEFAULTS.barChartOffsetY),
            barChartAxisPadding:     num('barChartAxisPadding',     CONFIG_DEFAULTS.barChartAxisPadding),
            chartTitleFontSize:      num('chartTitleFontSize',      CONFIG_DEFAULTS.chartTitleFontSize),
            chartLabelFontSize:      num('chartLabelFontSize',      CONFIG_DEFAULTS.chartLabelFontSize),
            barChartTitle:           txt('barChartTitle'),
            barChartXLabel:          txt('barChartXLabel'),
            barChartYLabel:          txt('barChartYLabel'),
        };
    }

    async _sendFrontendConfig(incoming, changedKeys) {
        if (this._configAbort) this._configAbort.abort();
        this._configAbort = new AbortController();

        console.log('[PathwayApp] Sending config update — changed:', changedKeys.join(', '));

        try {
            const res = await fetch('/api/update-config', {
                method:  'POST',
                headers: { 'Content-Type': 'application/json' },
                body:    JSON.stringify(incoming),
                signal:  this._configAbort.signal,
            });

            if (!res.ok) {
                const body   = await res.json().catch(() => ({}));
                const detail = body.details
                    ? body.details.join('; ')
                    : (body.error || `HTTP ${res.status}`);
                throw new Error(detail);
            }

            // Update owned config — rebuild barChart sub-object too
            Object.assign(this.config, incoming);
            this.config.barChart = {
                width:       this.config.barChartWidth       ?? CONFIG_DEFAULTS.barChartWidth,
                axisPadding: this.config.barChartAxisPadding ?? CONFIG_DEFAULTS.barChartAxisPadding,
                barHeight:   this.config.barHeight           ?? CONFIG_DEFAULTS.barHeight,
                gapBetween:  this.config.barChartGap         ?? CONFIG_DEFAULTS.barChartGap,
            };

            console.log('[PathwayApp] Config accepted — redrawing');
            this._redrawVisualization(changedKeys);
            showStatus('frontend-config-status', '✓ Settings applied', 'success');

        } catch (err) {
            if (err.name === 'AbortError') return;
            console.error('[PathwayApp] Config update error:', err.message);
            showStatus('frontend-config-status', `✗ ${err.message}`, 'error', 6000);
        }
    }

    // ─────────────────────────────────────────────────────────────────────
    // TARGETED REDRAW — only re-render what actually changed
    // ─────────────────────────────────────────────────────────────────────
    _redrawVisualization(changedKeys = []) {
        if (!this.initialJsonData || !this.visualizer) return;

        const affects = (...keys) => keys.some(k => changedKeys.includes(k));

        const layoutChanged = affects('metaboliteRadius', 'reactionRadius', 'nodeRadius');
        const imagesChanged = affects('imageSize');
        const chartsChanged = affects(
            'barChartWidth', 'barChartHeight', 'barHeight',
            'barChartOffsetY', 'barChartAxisPadding', 'barChartXLabel',
            'barChartTitle', 'chartLabelFontSize', 'chartTitleFontSize'
        );
        const labelsChanged = affects(
            'labelOffsetY', 'coproductLabelOffsetY',
            'metaboliteLabelFontSize', 'coproductLabelFontSize'
        );

        console.log('[PathwayApp] Targeted redraw —',
            { layoutChanged, imagesChanged, chartsChanged, labelsChanged });

        if (layoutChanged) {
            this.visualizer.disconnectRadiusObserver();
            this.visualizer.equalizeNodeRadii(this.config);
        }

        if (imagesChanged) {
            d3.selectAll('#map_container image').remove();
            this.visualizer.loadStructureImages(this.initialJsonData, this.config);
        }

        if (chartsChanged) {
            this.visualizer.clearBarCharts();
            this.visualizer.createNodeBarCharts(this.initialJsonData, this.config);
            this.visualizer.createSegmentBarCharts(
                this.initialJsonData[1].nodes,
                this.initialJsonData[1].reactions,
                this.config
            );
        }

        if (labelsChanged) {
            d3.selectAll('.metabolite-name').remove();
            d3.selectAll('.coproduct-name').remove();
            this.visualizer.initializeLabels(this.initialJsonData, this.config);
        }

        // Re-attach tooltips if charts or layout changed (elements were recreated)
        if (chartsChanged || layoutChanged) {
            this.visualizer.attachTooltipListeners(this.initialJsonData, this.config);
        }
    }

    // Full redraw (used after rebuild)
    _fullRedraw() {
        if (!this.initialJsonData || !this.visualizer) return;
        console.log('[PathwayApp] Full redraw');
        this.visualizer.disconnectRadiusObserver();
        d3.selectAll('.bar-chart').remove();
        d3.selectAll('.metabolite-name').remove();
        d3.selectAll('.coproduct-name').remove();
        d3.selectAll('#map_container image').remove();
        this.visualizer.initializeStructures(this.initialJsonData, this.config);
    }

    // ─────────────────────────────────────────────────────────────────────
    // ESCHER
    // ─────────────────────────────────────────────────────────────────────
    initializeEscher() {
        console.log('[PathwayApp] Initializing Escher');
        this.escherBuilder = escher.Builder(
            this.initialJsonData,
            null, null,
            d3.select('#map_container'),
            {
                enable_editing:            true,
                fill_screen:               false,
                reaction_styles:           ['abs', 'color', 'size', 'text'],
                never_ask_before_quit:     true,
                primary_metabolite_radius: 100,
                first_load_callback:       () => this.onEscherLoaded(),
            }
        );
    }

    onEscherLoaded() {
        console.log('[PathwayApp] Escher loaded — creating visualizer');

        // Create a fresh visualizer instance bound to this container
        this.visualizer = new EscherVisualizer('map_container');
        this.visualizer.initializeStructures(this.initialJsonData, this.config);
        this.visualizer.setupNaNLabelRemoval();

        this.setupExportButton();
        this.updateSubgraphStatus();

        if (this.highlightPath) {
            const pathNodes = this.highlightPath
                .split(',').map(s => s.trim()).filter(Boolean);
            if (pathNodes.length) {
                console.log('[PathwayApp] Applying path highlight —',
                    pathNodes.length, 'nodes');
                this.visualizer.highlightPathInPlace(
                    pathNodes,
                    this.initialJsonData[1],
                    this.config
                );
                const c = this._el('path-highlight-container');
                if (c) c.style.display = 'block';
            }
        }
    }

    rebuildEscher(newJsonData) {
        console.log('[PathwayApp] Rebuilding Escher');

        if (this.visualizer) {
            this.visualizer.destroy();
            this.visualizer = null;
        }

        this.initialJsonData = newJsonData;
        const container = this._el('map_container');
        if (container) container.innerHTML = '';

        this.escherBuilder = null;
        this.initializeEscher();
        this.updateSubgraphStatus();
    }

    updateSubgraphStatus() {
        const el = this._el('subgraph-status-container');
        if (el) {
            const isSubgraph = new URLSearchParams(window.location.search)
                .get('view') === 'subgraph';
            el.style.display = isSubgraph ? 'block' : 'none';
        }
    }

    // Called from template: window.pathwayApp.clearPathHighlight()
    clearPathHighlight() {
        console.log('[PathwayApp] Clearing path highlight');
        if (this.visualizer) this.visualizer.clearPathHighlight(this.config);

        const p = new URLSearchParams(window.location.search);
        if (p.has('highlight_path')) {
            p.delete('highlight_path');
            const newUrl = window.location.pathname
                + (p.toString() ? '?' + p.toString() : '');
            window.history.replaceState({}, '', newUrl);
        }

        const c = this._el('path-highlight-container');
        if (c) c.style.display = 'none';
    }

    // ─────────────────────────────────────────────────────────────────────
    // EXPORT
    // ─────────────────────────────────────────────────────────────────────
    setupExportButton() {
        const btn = this._el('export-svg-btn');
        if (btn) btn.addEventListener('click', () => this.exportMapAsSVG());
    }

    exportMapAsSVG() {
        console.log('[PathwayApp] Starting SVG export');
        try {
            const svgElement = document.querySelector('#map_container svg');
            if (!svgElement) { alert('No map SVG found to export'); return; }

            const clonedSvg = svgElement.cloneNode(true);
            clonedSvg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
            clonedSvg.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');

            Promise.all(
                Array.from(clonedSvg.querySelectorAll('image'))
                    .map(img => this._embedImage(img))
            )
            .then(() => {
                this._applyComputedStyles(svgElement, clonedSvg);
                this._cleanupEscherUI(clonedSvg);
                this._injectExportStyles(clonedSvg);

                const svgString = new XMLSerializer().serializeToString(clonedSvg);
                const baseName  = this._buildExportFilename();

                console.log('[PathwayApp] Downloading:', baseName);
                this._downloadBlob(
                    new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' }),
                    baseName + '.svg'
                );
                this._exportSvgAsPng(svgString, clonedSvg, baseName + '.png');
            })
            .catch(err => {
                console.error('[PathwayApp] Export error:', err);
                alert('Export error: ' + err.message);
            });
        } catch (err) {
            console.error('[PathwayApp] Export error:', err);
            alert('Export error: ' + err.message);
        }
    }

    _embedImage(imgElement) {
        const href = imgElement.getAttribute('xlink:href')
            || imgElement.getAttribute('href');
        if (!href || href.startsWith('data:')) return Promise.resolve();
        return new Promise(resolve => {
            const done = dataUri => {
                imgElement.setAttribute('xlink:href', dataUri);
                imgElement.removeAttribute('href');
                resolve();
            };
            const img = new Image();
            img.onload = () => {
                const canvas  = document.createElement('canvas');
                canvas.width  = img.width  || 200;
                canvas.height = img.height || 200;
                canvas.getContext('2d').drawImage(img, 0, 0);
                try { done(canvas.toDataURL('image/png')); } catch { resolve(); }
            };
            img.onerror = () => {
                fetch(href, { mode: 'cors', credentials: 'same-origin' })
                    .then(r => r.blob())
                    .then(blob => {
                        const reader = new FileReader();
                        reader.onload = () => done(reader.result);
                        reader.readAsDataURL(blob);
                    })
                    .catch(() => resolve());
            };
            img.src = href;
        });
    }

    _applyComputedStyles(origSvg, clonedSvg) {
        const origEls   = origSvg.querySelectorAll('*');
        const clonedEls = clonedSvg.querySelectorAll('*');
        const props = [
            'display','visibility','fill','stroke','stroke-width','stroke-dasharray',
            'stroke-linecap','stroke-linejoin','opacity','fill-opacity','stroke-opacity',
            'font-size','font-family','font-weight','text-anchor','dominant-baseline',
            'font-style','filter','marker-end','marker-start',
        ];
        const preserveNone = new Set([
            'fill','stroke','display','visibility','marker-end','marker-start',
        ]);
        origEls.forEach((orig, i) => {
            const clone = clonedEls[i];
            if (!clone) return;
            const computed = window.getComputedStyle(orig);
            const parts    = [];
            props.forEach(p => {
                const v = computed.getPropertyValue(p);
                if (!v || v === 'initial') return;
                if (v === 'none' && !preserveNone.has(p)) return;
                parts.push(`${p}: ${v}`);
            });
            if (parts.length) clone.setAttribute('style', parts.join('; '));
        });
    }

    _cleanupEscherUI(clonedSvg) {
        clonedSvg.querySelectorAll('.node-label.label')
            .forEach(el => el.setAttribute('style', 'display: none;'));
        clonedSvg.querySelectorAll(
            '.menu-bar,.button-panel,.search-bar-container,.dropdown-menu,foreignObject'
        ).forEach(el => el.remove());
    }

    _injectExportStyles(clonedSvg) {
        const style = document.createElement('style');
        style.innerHTML = `
            text { font-family: Arial, sans-serif; }
            .node-label.label { display: none !important; }
            .node-label.metabolite-name, .node-label.coproduct-name {
                font-size: 12px; font-weight: bold;
            }
            path.segment, .segment {
                fill: none !important; stroke: grey !important;
                stroke-width: 3px !important;
                stroke-dasharray: 5,5 !important; opacity: 0.3 !important;
            }
            .node-circle { fill-opacity: 0.8; }
            .menu-bar, .button-panel { display: none !important; }
            rect#canvas, rect#mouse-node { display: none !important; }
        `;
        clonedSvg.insertBefore(style, clonedSvg.firstChild);
    }

    _downloadBlob(blob, filename) {
        const url  = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href     = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
    }

    _exportSvgAsPng(svgString, clonedSvg, filename) {
        const vb = clonedSvg.getAttribute('viewBox');
        let w, h;
        if (vb) {
            const parts = vb.split(/[\s,]+/).map(Number);
            w = parts[2]; h = parts[3];
        } else {
            w = parseFloat(clonedSvg.getAttribute('width'))  || 1920;
            h = parseFloat(clonedSvg.getAttribute('height')) || 1080;
        }
        const scale  = 2;
        const canvas = document.createElement('canvas');
        canvas.width  = w * scale;
        canvas.height = h * scale;
        const ctx = canvas.getContext('2d');
        ctx.scale(scale, scale);
        const svgUrl = URL.createObjectURL(
            new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' })
        );
        const img = new Image();
        img.onload = () => {
            ctx.fillStyle = '#ffffff';
            ctx.fillRect(0, 0, w, h);
            ctx.drawImage(img, 0, 0, w, h);
            URL.revokeObjectURL(svgUrl);
            canvas.toBlob(blob => {
                if (blob) this._downloadBlob(blob, filename);
            }, 'image/png');
        };
        img.onerror = () => {
            console.error('[PathwayApp] PNG render failed');
            URL.revokeObjectURL(svgUrl);
        };
        img.src = svgUrl;
    }

    // ─────────────────────────────────────────────────────────────────────
    // FILENAME HELPERS
    // ─────────────────────────────────────────────────────────────────────
    _nodeName(id) {
        const n = this.availableNodes.find(x => x.id === id);
        return (n && n.name && n.name !== id) ? n.name : id;
    }

    _sanitize(str) {
        return str.replace(/[^a-zA-Z0-9-]/g, '_')
                  .replace(/_+/g, '_')
                  .replace(/^_|_$/g, '');
    }

    _buildExportFilename() {
        const p      = new URLSearchParams(window.location.search);
        const vert   = this.config.smallGraphLayoutVertical ?? false;
        const orient = vert ? 'vertical' : 'horizontal';
        const date   = new Date().toISOString().slice(0, 10);
        const isSubgraph = p.get('view') === 'subgraph';
        const start      = p.get('start')    || '';
        const end        = p.get('end')      || '';
        const selected   = p.get('selected') || '';
        if (isSubgraph && start && end) {
            return `pathway_${this._sanitize(this._nodeName(start))}`
                 + `-${this._sanitize(this._nodeName(end))}_${orient}_${date}`;
        }
        if (isSubgraph && selected) {
            const names = selected.split(',')
                .map(n => this._sanitize(this._nodeName(n.trim())));
            const str = names.length <= 3
                ? names.join('-')
                : names.slice(0, 3).join('-') + `_plus${names.length - 3}`;
            return `pathway_multi_${str}_${orient}_${date}`;
        }
        if (isSubgraph) return `pathway_subgraph_${orient}_${date}`;
        return `pathway_full_${orient}_${date}`;
    }
}

// =============================================================================
// BOOT
// =============================================================================
document.addEventListener('DOMContentLoaded', () => {
    const ns = window.PathwayApp;
    if (!ns) {
        console.error('[PathwayApp] window.PathwayApp not found');
        return;
    }

    const app = new PathwayApp(
        ns.initialJsonData,
        ns.CONFIG,
        ns.HIGHLIGHT_PATH,
        ns.VIEW_TYPE
    );

    app.initialize();
    window.pathwayApp = app;  // debug access only
    console.log('[PathwayApp] Boot complete');
});