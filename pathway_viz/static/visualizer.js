// visualizer.js
/**
 * EscherVisualizer
 *
 * Instantiated once per map load by PathwayApp.
 * Receives (jsonData, config) explicitly — reads nothing from window.
 *
 * Lifecycle:
 *   const viz = new EscherVisualizer('map_container');
 *   viz.initializeStructures(jsonData, config);
 *   viz.destroy();   // before rebuilding
 */
class EscherVisualizer {
    /**
     * @param {string} containerId  - id of the SVG container element
     */
    constructor(containerId) {
        this.containerId = containerId;
        this.container   = document.getElementById(containerId);

        this._observerRegistry   = [];
        this._radiusObserver     = null;
        this._radiusObserverActive = false;
        this._nanLabelInterval   = null;

        console.log('[EscherVisualizer] Instance created for:', containerId);
    }

    // =========================================================================
    //  OBSERVER REGISTRY
    // =========================================================================
    _registerObserver(obs) {
        this._observerRegistry.push(obs);
        return obs;
    }

    disconnectAllObservers() {
        this._observerRegistry.forEach(o => o.disconnect());
        this._observerRegistry     = [];
        this._radiusObserver       = null;
        this._radiusObserverActive = false;
        console.log('[EscherVisualizer] All observers disconnected');
    }

    disconnectRadiusObserver() {
        if (this._radiusObserver) {
            this._radiusObserver.disconnect();
            this._observerRegistry = this._observerRegistry
                .filter(o => o !== this._radiusObserver);
            this._radiusObserver       = null;
            this._radiusObserverActive = false;
        }
    }

    /** Full teardown — call before destroying the DOM */
    destroy() {
        this.disconnectAllObservers();
        this.stopNaNLabelRemoval();
        console.log('[EscherVisualizer] Destroyed');
    }

    // =========================================================================
    //  NaN LABEL REMOVAL
    // =========================================================================
    setupNaNLabelRemoval() {
        this.stopNaNLabelRemoval();  // clear any previous interval
        this._nanLabelInterval = setInterval(
            () => this.removeStoichiometryLabels(), 500
        );
    }

    stopNaNLabelRemoval() {
        if (this._nanLabelInterval) {
            clearInterval(this._nanLabelInterval);
            this._nanLabelInterval = null;
        }
    }

    removeStoichiometryLabels() {
        d3.selectAll('.stoichiometry-label').filter(function () {
            return d3.select(this).text() === 'NaN';
        }).remove();
    }

    // =========================================================================
    //  ENTRY POINT
    // =========================================================================
    initializeStructures(jsonData, config) {
        console.log('[EscherVisualizer] Initializing structures');
        this.clearBarCharts();
        this._removeTooltip();
        this.stylePathwayElements();
        this.equalizeNodeRadii(config);
        this.colourNodesByOrigin(config);
        this.loadStructureImages(jsonData, config);
        this.createNodeBarCharts(jsonData, config);
        this.createSegmentBarCharts(jsonData[1].nodes, jsonData[1].reactions, config);
        this.initializeLabels(jsonData, config);
        this.attachTooltipListeners(jsonData, config);
        console.log('[EscherVisualizer] Structures initialized');
    }

    // =========================================================================
    //  COLOUR HELPERS
    // =========================================================================
    _originColour(origin, config) {
        const map = config.originColours;
        return map[origin] || map.unknown;
    }

    _proteinColours(n) {
        const palette = [
            '#2a9d8f','#e76f51','#264653','#e9c46a',
            '#a8dadc','#f4a261','#457b9d','#e63946',
            '#06d6a0','#118ab2',
        ];
        return Array.from({ length: n }, (_, i) => palette[i % palette.length]);
    }

    // =========================================================================
    //  NODE COLOURING
    // =========================================================================
    colourNodesByOrigin(config) {
        d3.select('#' + this.containerId)
            .selectAll('.node-circle.metabolite-circle')
            .each((d, i, nodes) => {
                if (!d) return;
                const colour = this._originColour(d.origin || 'unknown', config);
                d3.select(nodes[i])
                    .style('fill',         colour)
                    .style('fill-opacity',  0.85)
                    .style('stroke',       d3.color(colour).darker(0.6).toString())
                    .style('stroke-width', '1.5px');
            });
    }

    // =========================================================================
    //  PATH HIGHLIGHTING
    // =========================================================================
    highlightPathInPlace(pathNodeIds, mapData, config) {
        console.log('[EscherVisualizer] Highlighting path —',
            pathNodeIds.length, 'nodes');

        const ph    = config.pathHighlight;
        const idSet = new Set(pathNodeIds);
        const nodes = mapData.nodes || {};
        const _key  = (a, b) => `${a}:${b}`;
        const pathSegKeys = new Set();

        Object.entries(nodes).forEach(([nid, nd]) => {
            if (nd.node_type !== 'midpoint') return;
            const fid = nd.from_node_id;
            const tid = nd.to_node_id;
            if (!fid || !tid) return;
            for (let i = 0; i < pathNodeIds.length - 1; i++) {
                const a = pathNodeIds[i];
                const b = pathNodeIds[i + 1];
                if ((fid === a && tid === b) || (fid === b && tid === a)) {
                    pathSegKeys.add(_key(fid, nid));
                    pathSegKeys.add(_key(nid, tid));
                    pathSegKeys.add(_key(b, nid));
                    pathSegKeys.add(_key(nid, a));
                }
            }
        });

        d3.select('#' + this.containerId)
            .selectAll('.node-circle.metabolite-circle')
            .each((d, i, nodes) => {
                if (!d) return;
                const onPath = idSet.has(d.bigg_id);
                const el     = d3.select(nodes[i]);
                const colour = this._originColour(d.origin || 'unknown', config);
                if (onPath) {
                    el.style('fill',         colour)
                      .style('fill-opacity',  1.0)
                      .style('stroke',        ph.nodeStrokeColour)
                      .style('stroke-width',  ph.nodeStrokeWidth + 'px');
                } else {
                    el.style('fill',         colour)
                      .style('fill-opacity',  ph.dimOpacity)
                      .style('stroke',       d3.color(colour).darker(0.6).toString())
                      .style('stroke-width', '1.5px');
                }
            });

        d3.select('#' + this.containerId)
            .selectAll('path.segment')
            .each((d, i, segs) => {
                if (!d) return;
                const key    = _key(d.from_node_id || '', d.to_node_id || '');
                const onPath = pathSegKeys.has(key);
                const el     = d3.select(segs[i]);
                if (onPath) {
                    el.style('stroke',           ph.segmentColour)
                      .style('stroke-width',      ph.segmentStrokeWidth + 'px')
                      .style('stroke-dasharray', 'none')
                      .style('opacity',           1.0);
                } else {
                    el.style('stroke',           'black')
                      .style('stroke-width',     '1px')
                      .style('stroke-dasharray', '5,5')
                      .style('opacity',           ph.dimOpacity);
                }
            });

        d3.select('#' + this.containerId)
            .classed('path-highlight-active', true);
    }

    clearPathHighlight(config) {
        console.log('[EscherVisualizer] Clearing path highlight');
        d3.select('#' + this.containerId)
            .classed('path-highlight-active', false);
        this.colourNodesByOrigin(config);
        d3.select('#' + this.containerId)
            .selectAll('path.segment')
            .style('stroke',           'black')
            .style('stroke-width',     '3px')
            .style('stroke-dasharray', '5,5')
            .style('opacity',           1.0);
    }

    // =========================================================================
    //  TOOLTIP
    // =========================================================================
    _getTooltipDiv() {
        let div = document.getElementById('escher-tooltip');
        if (!div) {
            div = document.createElement('div');
            div.id = 'escher-tooltip';
            Object.assign(div.style, {
                position:      'fixed',
                pointerEvents: 'auto',
                background:    'rgba(255,255,255,0.97)',
                border:        '1px solid #ccc',
                borderRadius:  '6px',
                padding:       '10px 14px',
                fontSize:      '12px',
                fontFamily:    'Arial, sans-serif',
                boxShadow:     '0 3px 12px rgba(0,0,0,0.18)',
                maxWidth:      '480px',
                maxHeight:     '60vh',
                overflowY:     'auto',
                zIndex:        '9999',
                display:       'none',
                lineHeight:    '1.5',
            });
            div.addEventListener('mouseleave', () => {
                div.style.display = 'none';
            });
            document.body.appendChild(div);
        }
        return div;
    }

    _removeTooltip() {
        const div = document.getElementById('escher-tooltip');
        if (div) div.style.display = 'none';
    }

    _positionTooltip(div, evt) {
        const margin = 14;
        const vpW    = window.innerWidth;
        const vpH    = window.innerHeight;
        const divW   = div.offsetWidth  || 360;
        const divH   = div.offsetHeight || 200;
        let left     = evt.clientX + margin;
        let top      = evt.clientY + margin;
        if (left + divW > vpW - margin) left = evt.clientX - divW - margin;
        if (top  + divH > vpH - margin) top  = evt.clientY - divH - margin;
        div.style.left    = Math.max(0, left) + 'px';
        div.style.top     = Math.max(0, top)  + 'px';
        div.style.display = 'block';
    }

    // ── Formatters ────────────────────────────────────────────────────────
    _fmt(v) {
        if (v === null || v === undefined) return '—';
        const n = Number(v);
        if (isNaN(n)) return String(v);
        if (n === 0)  return '0';
        return (Math.abs(n) >= 1e5 || (Math.abs(n) < 1e-2 && n !== 0))
            ? n.toExponential(3)
            : n.toPrecision(5);
    }

    _tr(label, value) {
        return `<tr>
          <td style="padding:2px 10px 2px 0;color:#666;
              white-space:nowrap;vertical-align:top">${label}</td>
          <td style="padding:2px 0;color:#111;
              word-break:break-word">${value}</td>
        </tr>`;
    }

    _replicateTable(replicates) {
        if (!replicates || !replicates.length) return '—';
        const cells = replicates.map(v =>
            `<td style="padding:1px 4px;border:1px solid #e0e0e0;
                background:#f8f8f8;color:#333;font-size:10px;
                white-space:nowrap">${this._fmt(v)}</td>`
        ).join('');
        return `<table style="border-collapse:collapse;
            display:inline-table"><tr>${cells}</tr></table>`;
    }

    _conditionBlock(entry) {
        return `
          <tr><td colspan="2" style="padding-top:5px;padding-bottom:1px;
              font-weight:bold;color:#2a6496;font-size:11px">
            ${entry.name}
          </td></tr>
          ${this._tr('Mean&nbsp;±&nbsp;SD',
              `${this._fmt(entry.mean)}&nbsp;±&nbsp;${this._fmt(entry.std_dev)}`)}
          ${this._tr('n', entry.count)}
          ${this._tr('Replicates', this._replicateTable(entry.replicates))}`;
    }

    _divRow() {
        return `<tr><td colspan="2">
          <hr style="margin:3px 0;border:none;border-top:1px solid #eee">
        </td></tr>`;
    }

    _metaboliteTooltipHTML(tt, config) {
        const origin = tt.origin || 'unknown';
        const colour = this._originColour(origin, config);
        const header = `
          <div style="border-bottom:2px solid ${colour};
              margin-bottom:6px;padding-bottom:5px;
              display:flex;align-items:baseline;gap:8px">
            <span style="font-size:13px;font-weight:bold;
                color:${colour};flex:1;word-break:break-word">
              ${tt.name || tt.id}
            </span>
            <span style="font-size:10px;color:#888;white-space:nowrap">${tt.id}</span>
            <span style="font-size:10px;background:#f0f0f0;border-radius:3px;
                padding:1px 5px;color:${colour};white-space:nowrap">${origin}</span>
          </div>`;

        if (!tt.conditions || !tt.conditions.length) {
            return header + '<div style="color:#999;font-style:italic">No abundance data</div>';
        }
        const rows = tt.conditions.map((c, i) =>
            (i > 0 ? this._divRow() : '') + this._conditionBlock(c)
        ).join('');
        return header
            + `<table style="border-collapse:collapse;width:100%;font-size:11px">${rows}</table>`;
    }

    _reactionTooltipHTML(tt) {
        const rxnId    = tt.reaction_id || 'unknown';
        const fromName = tt.from_node?.name || tt.from_node?.id || '?';
        const toName   = tt.to_node?.name   || tt.to_node?.id   || '?';
        const fromId   = tt.from_node?.id   || '';
        const toId     = tt.to_node?.id     || '';
        const header = `
          <div style="border-bottom:2px solid #e76f51;
              margin-bottom:6px;padding-bottom:5px">
            <span style="font-size:13px;font-weight:bold;color:#e76f51">${rxnId}</span>
          </div>
          <table style="border-collapse:collapse;width:100%;
              font-size:11px;margin-bottom:6px">
            ${this._tr('From',
                `<strong>${fromName}</strong>
                 <span style="color:#aaa;font-size:10px">&nbsp;(${fromId})</span>`)}
            ${this._tr('To',
                `<strong>${toName}</strong>
                 <span style="color:#aaa;font-size:10px">&nbsp;(${toId})</span>`)}
          </table>`;

        if (!tt.proteins || !tt.proteins.length) {
            return header + '<div style="color:#999;font-style:italic">No proteomics data</div>';
        }
        const colours = this._proteinColours(tt.proteins.length);
        const blocks  = tt.proteins.map((prot, pi) => {
            const colour  = colours[pi];
            const shortId = prot.protein_id.length > 45
                ? prot.protein_id.slice(0, 43) + '…'
                : prot.protein_id;
            const condRows = (prot.conditions || []).map((c, i) =>
                (i > 0 ? this._divRow() : '') + this._conditionBlock(c)
            ).join('');
            return `
              <div style="margin-top:8px;border-left:3px solid ${colour};padding-left:8px">
                <div title="${prot.protein_id}"
                    style="font-size:11px;font-weight:bold;color:${colour};
                    margin-bottom:3px;word-break:break-all">${shortId}</div>
                <table style="border-collapse:collapse;width:100%;
                    font-size:11px">${condRows}</table>
              </div>`;
        }).join('');
        return header + blocks;
    }

    _attachHover(el, getHTML) {
        const div = this._getTooltipDiv();
        el.addEventListener('mouseenter', evt => {
            div.innerHTML = getHTML();
            this._positionTooltip(div, evt);
        });
        el.addEventListener('mousemove', evt => {
            if (div.style.display !== 'none') this._positionTooltip(div, evt);
        });
        el.addEventListener('mouseleave', evt => {
            const related = evt.relatedTarget;
            if (related && (div.contains(related) || div === related)) return;
            div.style.display = 'none';
        });
    }

    attachTooltipListeners(jsonData, config) {
        const nodeTooltips = {};
        Object.entries(jsonData[1]?.nodes || {}).forEach(([nid, nd]) => {
            if (nd.tooltip) nodeTooltips[nid] = nd.tooltip;
        });

        const segTooltips = {};
        Object.values(jsonData[1]?.reactions || {}).forEach(rxn => {
            Object.values(rxn.segments || {}).forEach(seg => {
                if (seg.tooltip && seg.edge_type === 'reactant_edge') {
                    segTooltips[
                        (seg.from_node_id || '') + ':' + (seg.to_node_id || '')
                    ] = seg.tooltip;
                }
            });
        });

        d3.select('#' + this.containerId)
            .selectAll('.node-circle.metabolite-circle')
            .each((d, i, nodes) => {
                if (!d?.tooltip) return;
                const tt = d.tooltip;
                this._attachHover(nodes[i],
                    () => this._metaboliteTooltipHTML(tt, config));
            });

        d3.select('#' + this.containerId)
            .selectAll('.node-circle')
            .each((d, i, nodes) => {
                if (!d?.tooltip || d.node_type !== 'midpoint') return;
                const tt = d.tooltip;
                this._attachHover(nodes[i], () => this._reactionTooltipHTML(tt));
            });

        d3.select('#' + this.containerId)
            .selectAll('path.segment')
            .each((d, i, segs) => {
                if (!d) return;
                const key = (d.from_node_id || '') + ':' + (d.to_node_id || '');
                const tt  = segTooltips[key];
                if (!tt) return;
                this._attachHover(segs[i], () => this._reactionTooltipHTML(tt));
            });

        // Bar chart tooltips — attached after render, no rAF needed
        // (charts are created synchronously before this runs)
        d3.select('#' + this.containerId)
            .selectAll('.metabolite-bar-chart')
            .each((d, i, charts) => {
                const biggId = charts[i].dataset?.biggId;
                if (!biggId) return;
                const tt = nodeTooltips[biggId];
                if (!tt) return;
                this._attachHover(charts[i],
                    () => this._metaboliteTooltipHTML(tt, config));
            });

        d3.select('#' + this.containerId)
            .selectAll('.proteomics-bar-chart')
            .each((d, i, charts) => {
                const key = charts[i].dataset?.segKey;
                if (!key) return;
                const tt = segTooltips[key];
                if (!tt) return;
                this._attachHover(charts[i], () => this._reactionTooltipHTML(tt));
            });
    }

    // =========================================================================
    //  LABELS
    // =========================================================================
    addMetaboliteLabels(config) {
        d3.selectAll('g.node').each((d, i, nodes) => {
            const group  = d3.select(nodes[i]);
            const circle = group.select('.node-circle.metabolite-circle');
            if (circle.empty()) return;
            const data = circle.data()[0];
            if (!data) return;

            group.selectAll('.label').style('display', 'none');

            let label = group.select('.node-label.metabolite-name');
            if (label.empty()) {
                label = group.append('text')
                    .attr('class',           'node-label metabolite-name')
                    .style('font-family',    'Arial, sans-serif')
                    .style('font-weight',    'bold')
                    .style('fill',           '#333')
                    .style('text-anchor',    'middle')
                    .style('pointer-events', 'none');
            }

            const update = () => {
                label
                    .style('font-size', config.metaboliteLabelFontSize + 'px')
                    .text((data.name || data.bigg_id || 'Unknown').replace(/;\s*$/, ''));
                const transform = circle.attr('transform');
                if (transform) {
                    const m = transform.match(/translate\(([^,]+),([^)]+)\)/);
                    if (m) {
                        label.attr('transform',
                            `translate(${m[1]},${parseFloat(m[2]) + config.labelOffsetY})`);
                    }
                }
            };

            update();
            this._registerObserver(
                new MutationObserver(update)
            ).observe(circle.node(), { attributes: true, attributeFilter: ['transform'] });
        });
    }

    addCoproductLabels(config) {
        d3.selectAll('.coproduct-circle').each((d, i, circles) => {
            const circle     = d3.select(circles[i]);
            const data       = circle.data()[0];
            if (!data) return;
            const parentNode = d3.select(circles[i].parentNode);

            parentNode.selectAll('.label').style('display', 'none');

            let label = parentNode.select('.node-label.coproduct-name');
            if (label.empty()) {
                label = parentNode.append('text')
                    .attr('class',           'node-label coproduct-name')
                    .style('font-family',    'Arial, sans-serif')
                    .style('fill',           '#666')
                    .style('text-anchor',    'middle')
                    .style('pointer-events', 'none');
            }

            const update = () => {
                const text = (config.defaultToBiggId
                    ? (data.bigg_id || data.name || 'Coproduct')
                    : (data.name    || data.bigg_id || 'Coproduct')
                ).replace(/;\s*$/, '');
                label.style('font-size', config.coproductLabelFontSize + 'px').text(text);
                const transform = circle.attr('transform');
                if (transform) {
                    const m = transform.match(/translate\(([^,]+),([^)]+)\)/);
                    if (m) {
                        label.attr('transform',
                            `translate(${m[1]},${parseFloat(m[2]) - config.coproductLabelOffsetY})`);
                    }
                }
            };

            update();
            this._registerObserver(
                new MutationObserver(update)
            ).observe(circles[i], { attributes: true, attributeFilter: ['transform'] });
        });
    }

    initializeLabels(jsonData, config) {
        this.removeStoichiometryLabels();
        this.addMetaboliteLabels(config);
        this.addCoproductLabels(config);
    }

    // =========================================================================
    //  STRUCTURE IMAGES
    // =========================================================================
    loadStructureImages(jsonData, config) {
        const maxDisplaySize = config.imageSize;
        const isVertical     = config.smallGraphLayoutVertical ?? true;
        const nodeThreshold  = config.nodeThresholdSmall ?? 10;

        fetch('static/structure_imgs/image_dimensions.json')
            .then(r => r.ok ? r.json() : {})
            .catch(() => ({}))
            .then(dimensions => {
                const circles  = d3.select('#' + this.containerId)
                    .selectAll('.metabolite-circle');
                const numNodes = circles.size();
                const useLeft  = isVertical && numNodes < nodeThreshold;

                let maxNatW = Object.values(dimensions)
                    .reduce((m, d) => Math.max(m, d.w), 0);
                if (maxNatW === 0) maxNatW = 1;
                const scale = maxDisplaySize / maxNatW;

                console.log('[EscherVisualizer] Loading structure images —',
                    numNodes, 'nodes, useLeft:', useLeft);

                circles.each((data, i, nodes) => {
                    const circle     = d3.select(nodes[i]);
                    const parentNode = d3.select(nodes[i].parentNode);

                    if (data.highlight) circle.classed('highlighted', true);

                    const imgPath = `static/structure_imgs/${data.bigg_id}.png`;
                    const dim     = dimensions[data.bigg_id];
                    const dispW   = dim ? dim.w * scale : maxDisplaySize;
                    const dispH   = dim ? dim.h * scale : maxDisplaySize;

                    const getOffset = () => useLeft
                        ? { x: -dispW - 10, y: -dispH }
                        : { x: -dispW / 2,  y: -dispH - 10 };

                    const img = new Image();
                    img.onload = () => {
                        const o = getOffset();
                        parentNode.insert('image', 'text')
                            .attr('class',      'structure-image')
                            .attr('transform',
                                `translate(${data.x + o.x},${data.y + o.y})`)
                            .attr('width',  dispW)
                            .attr('height', dispH)
                            .attr('xlink:href', imgPath);
                    };
                    img.onerror = () => {};  // missing images silently ignored
                    img.src = imgPath;

                    const updatePos = () => {
                        const transform = circle.attr('transform');
                        if (!transform) return;
                        const o = getOffset();
                        parentNode.select('image').attr('transform',
                            `${transform} translate(${o.x},${o.y})`);
                    };
                    updatePos();

                    this._registerObserver(
                        new MutationObserver(updatePos)
                    ).observe(nodes[i], { attributes: true, attributeFilter: ['transform'] });
                });
            });
    }

    // =========================================================================
    //  BAR CHARTS – METABOLOMICS
    // =========================================================================
    _renderMetaboliteBarChart(parentNode, element, data, config, options = {}) {
        const bc = config.barChart;
        const cfg = Object.assign({
            chartWidth:  bc.width,
            barHeight:   bc.barHeight,
            axisPadding: bc.axisPadding,
            barColor:    this._originColour(data.origin || 'unknown', config),
            hoverColor:  '#1f7a67',
            getPosition: null,
        }, options);

        const gi = data.graph_info;
        if (!gi || typeof gi !== 'object' || Array.isArray(gi)) return;

        const conditions = Object.keys(gi).filter(
            k => k !== 'metabolite' && k !== 'KEGG_C_number'
        );
        if (!conditions.length) return;

        const values = conditions.map(k => gi[k]?.average ?? 0);
        const maxVal = Math.max(...values, 1e-9);
        const chartH = cfg.barHeight * conditions.length;
        const drawW  = cfg.chartWidth - cfg.axisPadding;
        const xScale = d3.scaleLinear().domain([0, maxVal]).range([0, drawW]);

        const pos = cfg.getPosition
            ? cfg.getPosition(data, cfg)
            : { x: (data.x || 0) - cfg.chartWidth, y: (data.y || 0) - chartH / 2 };

        const group = parentNode.append('g')
            .attr('class',     'bar-chart metabolite-bar-chart')
            .attr('transform', `translate(${pos.x},${pos.y})`);

        if (data.bigg_id) group.node().dataset.biggId = data.bigg_id;

        const updatePos = () => {
            const np = cfg.getPosition ? cfg.getPosition(data, cfg) : {
                x: (data.x || 0) - cfg.chartWidth,
                y: (data.y || 0) - chartH / 2,
            };
            group.attr('transform', `translate(${np.x},${np.y})`);
        };

        this._registerObserver(
            new MutationObserver(updatePos)
        ).observe(element.node(), { attributes: true, attributeFilter: ['transform'] });

        group.insert('rect', ':first-child')
            .attr('x',      -cfg.axisPadding - 4)
            .attr('y',      -8)
            .attr('width',   cfg.chartWidth + 8)
            .attr('height',  chartH + 16)
            .attr('fill',   'white')
            .attr('opacity', 0.88)
            .attr('rx', 4)
            .style('cursor', 'default');

        conditions.forEach((cond, i) => {
            group.append('text')
                .attr('x',  -4)
                .attr('y',   i * cfg.barHeight + cfg.barHeight / 2)
                .attr('dy', '0.35em')
                .style('text-anchor', 'end')
                .style('font-size',   config.chartLabelFontSize + 'px')
                .style('fill',        '#555')
                .text(cond);
        });

        group.append('g')
            .attr('class',     'x-axis')
            .attr('transform', `translate(0,${chartH})`)
            .call(d3.axisBottom(xScale).ticks(3).tickFormat(d3.format('.1e')))
            .selectAll('text')
            .style('font-size', config.chartLabelFontSize + 'px')
            .style('fill',      '#555');

        conditions.forEach((cond, i) => {
            const avg  = gi[cond]?.average ?? 0;
            const std  = gi[cond]?.std_dev ?? 0;
            const barW = xScale(avg);

            group.append('rect')
                .attr('class',  'bar')
                .attr('x',       0)
                .attr('y',       i * cfg.barHeight)
                .attr('width',   Math.max(barW, 0))
                .attr('height',  cfg.barHeight - 2)
                .attr('fill',    cfg.barColor)
                .on('mouseover', function () {
                    d3.select(this).attr('fill', cfg.hoverColor);
                    group.append('text').attr('class', 'value-label')
                        .attr('x',  barW + 3)
                        .attr('y',  i * cfg.barHeight + cfg.barHeight / 2)
                        .attr('dy', '0.35em')
                        .style('font-size',   config.chartLabelFontSize + 'px')
                        .style('fill',        '#111')
                        .style('font-weight', 'bold')
                        .text(`${d3.format('.1e')(avg)} ± ${d3.format('.1e')(std)}`);
                })
                .on('mouseout', function () {
                    d3.select(this).attr('fill', cfg.barColor);
                    group.selectAll('.value-label').remove();
                });
        });
    }

    // =========================================================================
    //  BAR CHARTS – PROTEOMICS
    // =========================================================================
    _renderProteomicsBarCharts(parentNode, segData, nodeData, config) {
        const bc = config.barChart;
        const proteinList = segData.graph_info;
        if (!Array.isArray(proteinList) || !proteinList.length) return;

        const conditions = Object.keys(proteinList[0]?.stats || {});
        if (!conditions.length) return;

        const nProteins = proteinList.length;
        const colours   = this._proteinColours(nProteins);
        const barH      = bc.barHeight;
        const axisPad   = bc.axisPadding;
        const chartW    = bc.width;
        const drawW     = chartW - axisPad;
        const singleH   = barH * conditions.length;
        const gap       = bc.gapBetween;

        let maxVal = 1e-9;
        proteinList.forEach(p =>
            conditions.forEach(c => {
                const v = p.stats?.[c]?.average ?? 0;
                if (v > maxVal) maxVal = v;
            })
        );
        const xScale = d3.scaleLinear().domain([0, maxVal]).range([0, drawW]);

        const fromNode = nodeData[segData.to_node_id];
        const anchorX  = (fromNode?.x ?? 0) - chartW - 40;
        const anchorY  = (fromNode?.y ?? 0)
            - (nProteins * singleH + (nProteins - 1) * gap) / 2;

        const segKey     = (segData.from_node_id || '') + ':' + (segData.to_node_id || '');
        const outerGroup = parentNode.append('g')
            .attr('class',     'bar-chart proteomics-bar-chart')
            .attr('transform', `translate(${anchorX},${anchorY})`);
        outerGroup.node().dataset.segKey = segKey;

        const totalH = nProteins * singleH + (nProteins - 1) * gap;
        outerGroup.insert('rect', ':first-child')
            .attr('x',      -axisPad - 4)
            .attr('y',      -8)
            .attr('width',   chartW + 8)
            .attr('height',  totalH + 16)
            .attr('fill',   'white')
            .attr('opacity', 0.88)
            .attr('rx', 4)
            .style('cursor', 'default');

        proteinList.forEach((prot, pi) => {
            const colour  = colours[pi];
            const offsetY = pi * (singleH + gap);
            const pGroup  = outerGroup.append('g')
                .attr('transform', `translate(0,${offsetY})`);

            const shortId = prot.protein_id.length > 28
                ? prot.protein_id.slice(0, 26) + '…'
                : prot.protein_id;

            pGroup.append('text')
                .attr('x',  drawW / 2)
                .attr('y',  -4)
                .style('text-anchor', 'middle')
                .style('font-size',   config.chartTitleFontSize + 'px')
                .style('font-weight', 'bold')
                .style('fill',        colour)
                .text(shortId)
                .append('title').text(prot.protein_id);

            if (pi > 0) {
                pGroup.append('line')
                    .attr('x1', -axisPad).attr('x2', chartW - axisPad + 4)
                    .attr('y1', -gap / 2 - singleH)
                    .attr('y2', -gap / 2 - singleH)
                    .attr('stroke', '#ddd').attr('stroke-width', 1);
            }

            conditions.forEach((cond, ci) => {
                pGroup.append('text')
                    .attr('x',  -4)
                    .attr('y',   ci * barH + barH / 2)
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'end')
                    .style('font-size',   config.chartLabelFontSize + 'px')
                    .style('fill',        '#555')
                    .text(cond);
            });

            if (pi === nProteins - 1) {
                pGroup.append('g')
                    .attr('class',     'x-axis')
                    .attr('transform', `translate(0,${singleH})`)
                    .call(d3.axisBottom(xScale).ticks(3).tickFormat(d3.format('.1e')))
                    .selectAll('text')
                    .style('font-size', config.chartLabelFontSize + 'px')
                    .style('fill',      '#555');

                if (config.barChartXLabel) {
                    pGroup.append('text')
                        .attr('x',  drawW / 2)
                        .attr('y',  singleH + 28)
                        .style('text-anchor', 'middle')
                        .style('font-size',   config.chartLabelFontSize + 'px')
                        .style('fill',        '#666')
                        .text(config.barChartXLabel);
                }
            }

            conditions.forEach((cond, ci) => {
                const avg = prot.stats?.[cond]?.average ?? null;
                const std = prot.stats?.[cond]?.std_dev ?? 0;
                if (avg === null) return;
                const barW = xScale(avg);

                pGroup.append('rect')
                    .attr('x',       0)
                    .attr('y',       ci * barH)
                    .attr('width',   Math.max(barW, 0))
                    .attr('height',  barH - 2)
                    .attr('fill',    colour)
                    .on('mouseover', function () {
                        d3.select(this).attr('opacity', 0.7);
                        pGroup.append('text').attr('class', 'value-label')
                            .attr('x',  barW + 3)
                            .attr('y',  ci * barH + barH / 2)
                            .attr('dy', '0.35em')
                            .style('font-size',   config.chartLabelFontSize + 'px')
                            .style('fill',        '#111')
                            .style('font-weight', 'bold')
                            .text(`${d3.format('.1e')(avg)} ± ${d3.format('.1e')(std)}`);
                    })
                    .on('mouseout', function () {
                        d3.select(this).attr('opacity', 1);
                        pGroup.selectAll('.value-label').remove();
                    });
            });
        });
    }

    // =========================================================================
    //  NODE BAR CHARTS
    // =========================================================================
    createNodeBarCharts(jsonData, config) {
        const isVert    = config.smallGraphLayoutVertical ?? true;
        const threshold = config.nodeThresholdSmall ?? 10;
        const numNodes  = d3.select('#' + this.containerId)
            .selectAll('.metabolite-circle').size();
        const useLeft   = isVert && numNodes < threshold;

        console.log('[EscherVisualizer] Creating node bar charts — useLeft:', useLeft);

        d3.select('#' + this.containerId)
            .selectAll('.node-circle.metabolite-circle')
            .each((data, i, nodes) => {
                if (!data?.graph_info || Array.isArray(data.graph_info)) return;
                const element    = d3.select(nodes[i]);
                const parentNode = d3.select(nodes[i].parentNode);
                this._renderMetaboliteBarChart(
                    parentNode, element, data, config, {
                        getPosition: (d, cfg) => useLeft
                            ? { x: d.x + config.barChartOffsetY,
                                y: d.y - cfg.chartWidth }
                            : { x: d.x - config.barChart.width / 2,
                                y: d.y + config.barChartOffsetY },
                    }
                );
            });
    }

    // =========================================================================
    //  SEGMENT BAR CHARTS
    // =========================================================================
    createSegmentBarCharts(nodeData, reactions, config) {
        const segments = {};
        Object.values(reactions || {}).forEach(rxn => {
            Object.values(rxn.segments || {}).forEach(seg => {
                if (
                    seg.edge_type === 'reactant_edge' &&
                    Array.isArray(seg.graph_info) &&
                    seg.graph_info.length > 0
                ) {
                    segments[
                        (seg.from_node_id || '') + ':' + (seg.to_node_id || '')
                    ] = seg;
                }
            });
        });

        console.log('[EscherVisualizer] Creating segment bar charts —',
            Object.keys(segments).length, 'proteomics segments');

        d3.select('#' + this.containerId)
            .selectAll('.segment')
            .each((data, i, segs) => {
                if (!data) return;
                const key = (data.from_node_id || '') + ':' + (data.to_node_id || '');
                const seg = segments[key];
                if (!seg) return;
                const parentNode = d3.select(segs[i].parentNode);
                this._renderProteomicsBarCharts(parentNode, seg, nodeData, config);
            });
    }

    // =========================================================================
    //  NODE STYLING
    // =========================================================================
    stylePathwayElements() {
        d3.selectAll('.segment')
            .style('stroke',           'black')
            .style('stroke-width',     '3px')
            .style('stroke-dasharray', '5,5');
    }

    equalizeNodeRadii(config) {
        const metR  = config.metaboliteRadius ?? config.nodeRadius;
        const reacR = config.reactionRadius   ?? config.nodeRadius;

        const update = () => {
            d3.selectAll('circle.node-circle.metabolite-circle')
                .each(function () { this.setAttribute('r', metR); });
            d3.selectAll('circle.coproduct-circle')
                .each(function () { this.setAttribute('r', reacR); });
        };

        update();

        if (this.container && !this._radiusObserverActive) {
            const obs = new MutationObserver(() => {
                let needs = false;
                d3.selectAll('circle.node-circle.metabolite-circle')
                    .each(function () {
                        if (parseFloat(this.getAttribute('r')) !== metR) needs = true;
                    });
                if (needs) update();
            });
            obs.observe(this.container, {
                subtree: true, attributes: true, attributeFilter: ['r'],
            });
            this._radiusObserver       = obs;
            this._radiusObserverActive = true;
            this._registerObserver(obs);
        }
    }

    clearBarCharts() {
        d3.selectAll('.bar-chart').remove();
    }

    // =========================================================================
    //  LEGACY HIGHLIGHTING
    // =========================================================================
    highlightPath(pathNodes) {
        d3.selectAll('circle.node-circle.metabolite-circle')
            .classed('highlighted-path', function () {
                const d = d3.select(this).data()[0];
                return d && pathNodes.includes(d.bigg_id);
            });
        d3.selectAll('path.segment')
            .classed('highlighted-path', function () {
                const d = d3.select(this).data()[0];
                return d && pathNodes.includes(d.to_node_id);
            });
    }

    highlightMultiNodes(selectedNodeIds) {
        d3.selectAll('circle.node-circle.metabolite-circle')
            .classed('highlighted-multi', function () {
                const d = d3.select(this).data()[0];
                return d && selectedNodeIds.includes(d.bigg_id);
            });
        d3.selectAll('path.segment')
            .classed('highlighted-multi', function () {
                const d = d3.select(this).data()[0];
                return d && selectedNodeIds.includes(d.to_node_id);
            });
    }
}