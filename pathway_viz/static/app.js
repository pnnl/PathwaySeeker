// app.js
/**
 * PathwayApp - Sidebar, node selection, Escher init,
 * frontend config auto-apply, view-state injection, SVG export.
 */
const DEBUG = false;
const log = (...args) => DEBUG && console.log('[APP]', ...args);

// =============================================================================
// UTILITIES
// =============================================================================
function showStatus(elementId, message, type = 'info', autoDismissMs = 3000) {
    const el = document.getElementById(elementId);
    if (!el) return;
    el.className = `status-message status-${type}`;
    el.textContent = message;
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

// =============================================================================
// APP
// =============================================================================
class PathwayApp {
    constructor() {
        this.escherBuilder  = null;
        this.availableNodes = [];

        // AbortController for in-flight config AJAX requests.
        // Cancelled and replaced each time a new request is made so that
        // rapid slider changes never pile up on the server.
        this._configAbort = null;

        // All MutationObservers created by this instance, keyed by a label
        // so they can be disconnected cleanly before a rebuild.
        this._observers = new Map();

        // Cache DOM references that are accessed repeatedly.
        this._els = {};
    }

    // ── DOM element cache ─────────────────────────────────────────────────
    _el(id) {
        if (!this._els[id]) this._els[id] = document.getElementById(id);
        return this._els[id];
    }

    // ── Observer registry ─────────────────────────────────────────────────
    _addObserver(label, observer) {
        // Disconnect any previous observer registered under the same label
        if (this._observers.has(label)) {
            this._observers.get(label).disconnect();
        }
        this._observers.set(label, observer);
    }

    _disconnectAllObservers() {
        this._observers.forEach(obs => obs.disconnect());
        this._observers.clear();
    }

    // ─────────────────────────────────────────────────────────────────────────
    initialize() {
        this.setupSidebarToggle();
        this.setupKeyboardIsolation();
        this.populateNodeSelections();
        this.setupMultiNodeSelector();
        this.setupFrontendConfigAutoApply();
        this.setupBackendConfigViewState();
        if (window.initialJsonData) {
            this.initializeEscher();
        } else {
            console.error('[APP] No JSON data available');
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // SIDEBAR
    // ─────────────────────────────────────────────────────────────────────────
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

    // ─────────────────────────────────────────────────────────────────────────
    // NODE SELECTIONS
    // ─────────────────────────────────────────────────────────────────────────
    populateNodeSelections() {
        fetch('/api/nodes')
            .then(r => {
                if (!r.ok) throw new Error(`HTTP ${r.status}`);
                return r.json();
            })
            .then(nodes => {
                this.availableNodes = nodes;
                this._fillSelect('start-node', nodes, '-- Select start node --');
                this._fillSelect('end-node',   nodes, '-- Select end node --');
                this._fillMultiSelect(nodes);
                this.restoreSelectionFromUrl();
            })
            .catch(err => console.error('[APP] Error loading nodes:', err));
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

        const syncHidden = () => {
            hidden.value = Array.from(multi.selectedOptions)
                .map(o => o.value).join(',');
        };

        multi.addEventListener('change', syncHidden);

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

    // ─────────────────────────────────────────────────────────────────────────
    // BACKEND CONFIG – inject view-state into hidden fields before submit
    // ─────────────────────────────────────────────────────────────────────────
    setupBackendConfigViewState() {
        const form = this._el('backend-config-form');
        if (!form) return;

        // Yellow highlight on changed-but-not-submitted inputs
        form.querySelectorAll('.backend-config-input').forEach(input => {
            input.addEventListener('change', () => {
                input.classList.add('backend-input-changed');
            });
        });

        form.addEventListener('submit', () => {
            const p          = new URLSearchParams(window.location.search);
            const isSubgraph = p.get('view') === 'subgraph';

            const set = (id, val) => {
                const el = document.getElementById(id);
                if (el) el.value = val;
            };

            set('hidden-view-type',            isSubgraph ? 'subgraph' : 'full');
            set('hidden-start-node',           p.get('start')    || '');
            set('hidden-end-node',             p.get('end')      || '');
            set('hidden-path-nodes',           p.get('nodes')    || '');
            set('hidden-selected-nodes',       p.get('selected') || '');
            set('hidden-connection-distance',  p.get('dist')     || '');
            set('hidden-keep-positions',       p.get('keep_pos') || '1');

            log('Backend form submitting with view state:', {
                view:  isSubgraph ? 'subgraph' : 'full',
                start: p.get('start'),
                end:   p.get('end'),
            });
        });
    }

    // ─────────────────────────────────────────────────────────────────────────
    // FRONTEND CONFIG – auto-apply on change (debounced, no page reload)
    // ─────────────────────────────────────────────────────────────────────────
    setupFrontendConfigAutoApply() {
        const inputs = document.querySelectorAll('.frontend-config-input');
        if (!inputs.length) return;

        const applyConfig = debounce(() => {
            const config = this._collectFrontendConfig();
            log('Auto-applying frontend config:', config);
            this._sendFrontendConfig(config);
        }, 800);

        inputs.forEach(input => {
            input.addEventListener('input',  applyConfig);
            input.addEventListener('change', applyConfig);
        });
    }

    _collectFrontendConfig() {
        const num = (id, fallback) => {
            const el = this._el(id);
            return el ? (parseInt(el.value, 10) || fallback) : fallback;
        };
        const txt = (id, fallback = '') => {
            const el = this._el(id);
            return el ? el.value : fallback;
        };

        return {
            nodeRadius:              num('nodeRadius', 10),
            metaboliteRadius:        num('metaboliteRadius', 10),
            reactionRadius:          num('reactionRadius', 8),
            imageSize:               num('imageSize', 400),
            labelOffsetY:            num('labelOffsetY', 35),
            coproductLabelOffsetY:   num('coproductLabelOffsetY', 25),
            metaboliteLabelFontSize: num('metaboliteLabelFontSize', 14),
            coproductLabelFontSize:  num('coproductLabelFontSize', 10),
            barChartWidth:           num('barChartWidth', 100),
            barChartHeight:          num('barChartHeight', 100),
            barHeight:               num('barHeight', 15),
            barChartOffsetY:         num('barChartOffsetY', 60),
            barChartAxisPadding:     num('barChartAxisPadding', 20),
            chartTitleFontSize:      num('chartTitleFontSize', 14),
            chartLabelFontSize:      num('chartLabelFontSize', 14),
            barChartTitle:           txt('barChartTitle'),
            barChartXLabel:          txt('barChartXLabel', 'Abundance'),
            barChartYLabel:          txt('barChartYLabel'),
        };
    }

    async _sendFrontendConfig(config) {
        // Cancel any in-flight request before starting a new one
        if (this._configAbort) {
            this._configAbort.abort();
        }
        this._configAbort = new AbortController();

        try {
            const res = await fetch('/api/update-config', {
                method:  'POST',
                headers: { 'Content-Type': 'application/json' },
                body:    JSON.stringify(config),
                signal:  this._configAbort.signal,
            });

            if (!res.ok) {
                const body = await res.json().catch(() => ({}));
                const detail = body.details
                    ? body.details.join('; ')
                    : (body.error || `HTTP ${res.status}`);
                throw new Error(detail);
            }

            Object.assign(window.CONFIG, config);
            this._applyConfigToVisualizer(config);
            this._redrawVisualization();
            showStatus('frontend-config-status', '✓ Settings applied', 'success');
        } catch (err) {
            if (err.name === 'AbortError') {
                // Superseded by a newer request — silent ignore
                return;
            }
            console.error('[APP] Config update error:', err);
            showStatus('frontend-config-status', `✗ ${err.message}`, 'error', 6000);
        }
    }

    _applyConfigToVisualizer(config) {
        if (typeof EscherVisualizer === 'undefined') return;
        // window.CONFIG has already been updated via Object.assign above.
        // This explicit sync keeps barChart sub-values in step for clarity.
        window.CONFIG.barChartWidth       = config.barChartWidth;
        window.CONFIG.barChartHeight      = config.barChartHeight;
        window.CONFIG.barHeight           = config.barHeight;
        window.CONFIG.barChartAxisPadding = config.barChartAxisPadding;
    }

    _redrawVisualization() {
        if (!window.initialJsonData || typeof EscherVisualizer === 'undefined') return;

        // Disconnect old radius observer so it doesn't fight new values
        if (EscherVisualizer._radiusObserver) {
            EscherVisualizer._radiusObserver.disconnect();
            EscherVisualizer._radiusObserver       = null;
            EscherVisualizer._radiusObserverActive = false;
        }

        // Clear overlays
        d3.selectAll('.bar-chart').remove();
        d3.selectAll('.metabolite-name').remove();
        d3.selectAll('.coproduct-name').remove();
        d3.selectAll('#map_container image').remove();

        // Apply new radii (observer will be re-created by equalizeNodeRadii)
        const { metaboliteRadius, reactionRadius } = EscherVisualizer.CONFIG;
        d3.selectAll('circle.node-circle.metabolite-circle')
            .each(function() { this.setAttribute('r', metaboliteRadius); });
        d3.selectAll('circle.coproduct-circle')
            .each(function() { this.setAttribute('r', reactionRadius); });

        // Redraw all overlays (also re-creates radius observer with fresh values)
        EscherVisualizer.initializeStructures(window.initialJsonData);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // ESCHER
    // ─────────────────────────────────────────────────────────────────────────
    initializeEscher() {
        this.escherBuilder = escher.Builder(
            window.initialJsonData,
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
        EscherVisualizer.initializeStructures(window.initialJsonData);
        EscherVisualizer.setupNaNLabelRemoval();
        this.setupExportButton();
        this.updateSubgraphStatus();
    }

    rebuildEscher(newJsonData) {
        // Tear down all observers registered by previous Escher/visualizer calls
        this._disconnectAllObservers();
        EscherVisualizer.disconnectAllObservers();

        window.initialJsonData = newJsonData;

        const container = this._el('map_container');
        if (container) container.innerHTML = '';

        d3.selectAll('.bar-chart').remove();
        d3.selectAll('.metabolite-name').remove();
        d3.selectAll('.coproduct-name').remove();

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

    // ─────────────────────────────────────────────────────────────────────────
    // EXPORT
    // ─────────────────────────────────────────────────────────────────────────
    setupExportButton() {
        const btn = this._el('export-svg-btn');
        if (btn) btn.addEventListener('click', () => this.exportMapAsSVG());
    }

    exportMapAsSVG() {
        try {
            const svgElement = document.querySelector('#map_container svg');
            if (!svgElement) { alert('No map SVG found to export'); return; }

            const clonedSvg = svgElement.cloneNode(true);
            clonedSvg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
            clonedSvg.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');

            const imagePromises = Array.from(clonedSvg.querySelectorAll('image'))
                .map(img => this._embedImage(img));

            Promise.all(imagePromises)
                .then(() => {
                    this._applyComputedStyles(svgElement, clonedSvg);
                    this._cleanupEscherUI(clonedSvg);
                    this._injectExportStyles(clonedSvg);

                    const svgString = new XMLSerializer().serializeToString(clonedSvg);
                    const baseName  = this._buildExportFilename();

                    this._downloadBlob(
                        new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' }),
                        baseName + '.svg'
                    );
                    this._exportSvgAsPng(svgString, clonedSvg, baseName + '.png');
                })
                .catch(err => {
                    console.error('[APP] Export error:', err);
                    alert('Export error: ' + err.message);
                });
        } catch (err) {
            console.error('[APP] Export error:', err);
            alert('Export error: ' + err.message);
        }
    }

    _embedImage(imgElement) {
        const href = imgElement.getAttribute('xlink:href') || imgElement.getAttribute('href');
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
            'fill', 'stroke', 'display', 'visibility', 'marker-end', 'marker-start',
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
                fill: none !important;
                stroke: grey !important;
                stroke-width: 3px !important;
                stroke-dasharray: 5,5 !important;
                opacity: 0.3 !important;
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
            console.error('[APP] PNG conversion failed');
            URL.revokeObjectURL(svgUrl);
        };
        img.src = svgUrl;
    }

    // ─────────────────────────────────────────────────────────────────────────
    // FILENAME HELPERS
    // ─────────────────────────────────────────────────────────────────────────
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
        const vert   = window.CONFIG?.smallGraphLayoutVertical ?? false;
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
    const app = new PathwayApp();
    app.initialize();
    window.pathwayApp = app;
});