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
    /**
     * Colour metabolite nodes and midpoint nodes by their omics origin.
     *
     * @param {object}      config           - vis config (has originColours)
     * @param {object|null} originOverrides  - { bigg_id: 'metabolomics'|'proteomics'|'both'|'unknown' }
     *                                         Overrides d.origin for metabolite nodes.
     * @param {object|null} midpointOrigins  - { node_id: 'proteomics' }
     *                                         Origin for midpoint (reaction) nodes.
     */
    colourNodesByOrigin(config, originOverrides, midpointOrigins) {
        // ── Metabolite nodes ──────────────────────────────────────────────
        d3.select('#' + this.containerId)
            .selectAll('.node-circle.metabolite-circle')
            .each((d, i, nodes) => {
                if (!d) return;
                const origin = (originOverrides && originOverrides[d.bigg_id])
                    || d.origin
                    || 'unknown';
                const colour = this._originColour(origin, config);
                d3.select(nodes[i])
                    .style('fill',         colour)
                    .style('fill-opacity',  0.85)
                    .style('stroke',       d3.color(colour).darker(0.6).toString())
                    .style('stroke-width', '1.5px');
            });

        // ── Midpoint (reaction) nodes ─────────────────────────────────────
        if (midpointOrigins && Object.keys(midpointOrigins).length) {
            d3.select('#' + this.containerId)
                .selectAll('.node-circle')
                .each((d, i, nodes) => {
                    if (!d || d.node_type !== 'midpoint') return;
                    // d.bigg_id for midpoints is "midpoint_<id>"; use the Escher node id
                    // The D3 datum for a midpoint circle has node_id or we can derive it
                    // from bigg_id: "midpoint_<id>" -> id is the map node key
                    // However, midpointOrigins is keyed by the map node ID (string).
                    // d.bigg_id = "midpoint_<nodeId>" so strip the prefix.
                    const rawId = d.bigg_id
                        ? d.bigg_id.replace(/^midpoint_/, '')
                        : null;
                    const origin = (rawId && midpointOrigins[rawId])
                        || (d.bigg_id && midpointOrigins[d.bigg_id])
                        || null;
                    if (!origin) return;
                    const colour = this._originColour(origin, config);
                    d3.select(nodes[i])
                        .style('fill',         colour)
                        .style('fill-opacity',  0.85)
                        .style('stroke',       d3.color(colour).darker(0.6).toString())
                        .style('stroke-width', '1.5px');
                });
        }
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

        // New list format: tooltip.rows = [{metabolite_name, conditions: [...]}]
        if (Array.isArray(tt.rows) && tt.rows.length) {
            const rowBlocks = tt.rows.map((row, ri) => {
                const rowName = row.metabolite_name || '';
                const rowHeader = rowName
                    ? `<div style="font-size:11px;font-weight:bold;color:#444;
                           margin-top:${ri > 0 ? 10 : 0}px;margin-bottom:3px;
                           border-top:${ri > 0 ? '1px solid #ddd' : 'none'};
                           padding-top:${ri > 0 ? 6 : 0}px">${rowName}</div>`
                    : '';
                const conds = (row.conditions || []).filter(c =>
                    c.replicates && c.replicates.length > 0
                );
                if (!conds.length) return '';
                const condRows = conds.map((c, i) =>
                    (i > 0 ? this._divRow() : '') + this._conditionBlock(c)
                ).join('');
                return rowHeader
                    + `<table style="border-collapse:collapse;width:100%;font-size:11px">${condRows}</table>`;
            }).filter(Boolean).join('');
            return header + (rowBlocks || '<div style="color:#999;font-style:italic">No abundance data</div>');
        }

        // Legacy flat format: tooltip.conditions = [...]
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

    /**
     * Render one panel of a metabolite bar chart for a single CSV row.
     *
     * @param {d3.Selection} parentNode  - SVG group to append into
     * @param {object}       rowEntry    - {metabolite_name, conditions: {cond:{average,std_dev,count}}}
     * @param {object}       data        - node data (for position, bigg_id, origin)
     * @param {object}       config      - vis config
     * @param {object}       cfgOverride - merged on top of default bar-chart cfg
     * @param {number}       yOffset     - vertical offset for stacking multiple panels
     * @returns {number}  height of this panel (so caller can stack)
     */
    _renderMetaboliteBarChartPanel(parentNode, rowEntry, data, config, cfgOverride, yOffset) {
        const bc  = config.barChart;
        const cfg = Object.assign({
            chartWidth:  bc.width,
            barHeight:   bc.barHeight,
            axisPadding: bc.axisPadding,
            barColor:    this._originColour(data.origin || 'unknown', config),
            hoverColor:  '#1f7a67',
        }, cfgOverride);

        const conditionsObj = rowEntry.conditions || {};
        const minCount = config.barMinCount ?? 1;
        // Only include conditions that have a real average (not null/undefined/NaN)
        // AND whose replicate count meets the minimum threshold
        const conditions = Object.keys(conditionsObj).filter(k => {
            const entry = conditionsObj[k];
            const avg   = entry?.average;
            if (avg === null || avg === undefined || !isFinite(avg)) return false;
            const cnt = entry?.count ?? entry?.n ?? Infinity;
            return cnt >= minCount;
        });
        if (!conditions.length) return 0;

        const values  = conditions.map(k => conditionsObj[k]?.average ?? 0);
        const stdDevs = conditions.map(k => conditionsObj[k]?.std_dev ?? 0);
        const maxVal  = Math.max(...values.map((v, i) => v + stdDevs[i]), 1e-9);
        const chartH = cfg.barHeight * conditions.length;
        const drawW  = cfg.chartWidth - cfg.axisPadding;
        const xScale = d3.scaleLinear().domain([0, maxVal]).range([0, drawW]);

        const panel = parentNode.append('g')
            .attr('transform', `translate(0,${yOffset})`);

        // Row label (metabolite name) above the panel
        // Show only the method part if name is "MetaboliteName (method)"
        const rowName = rowEntry.metabolite_name || '';
        const methodMatch = rowName.match(/\(([^)]+)\)$/);
        const displayLabel = methodMatch ? methodMatch[1] : rowName;
        const labelH  = rowName ? cfg.barHeight : 0;
        if (rowName) {
            panel.append('text')
                .attr('x',  drawW / 2)
                .attr('y',  labelH / 2)
                .attr('dy', '0.35em')
                .style('text-anchor', 'middle')
                .style('font-size',   (config.chartLabelFontSize + 1) + 'px')
                .style('font-weight', 'bold')
                .style('fill',        '#444')
                .text(displayLabel)
                .append('title').text(rowName);
        }

        const barsGroup = panel.append('g')
            .attr('transform', `translate(0,${labelH})`);

        conditions.forEach((cond, i) => {
            barsGroup.append('text')
                .attr('x',  -4)
                .attr('y',   i * cfg.barHeight + cfg.barHeight / 2)
                .attr('dy', '0.35em')
                .style('text-anchor', 'end')
                .style('font-size',   config.chartLabelFontSize + 'px')
                .style('fill',        '#555')
                .text(cond);
        });

        barsGroup.append('g')
            .attr('class',     'x-axis')
            .attr('transform', `translate(0,${chartH})`)
            .call(d3.axisBottom(xScale).ticks(3).tickFormat(d => d3.format('.1e')(d).replace(/\.0(e)/, '$1')))
            .selectAll('text')
            .style('font-size', config.chartLabelFontSize + 'px')
            .style('fill',      '#555');

        conditions.forEach((cond, i) => {
            const avg  = conditionsObj[cond]?.average ?? 0;
            const std  = conditionsObj[cond]?.std_dev ?? 0;
            const barW = xScale(avg);
            const barColor = cfg.barColor;
            const hoverColor = cfg.hoverColor;

            barsGroup.append('rect')
                .attr('class',  'bar')
                .attr('x',       0)
                .attr('y',       i * cfg.barHeight)
                .attr('width',   Math.max(barW, 0))
                .attr('height',  cfg.barHeight - 2)
                .attr('fill',    barColor)
                .on('mouseover', function () {
                    d3.select(this).attr('fill', hoverColor);
                    barsGroup.append('text').attr('class', 'value-label')
                        .attr('x',  barW + 3)
                        .attr('y',  i * cfg.barHeight + cfg.barHeight / 2)
                        .attr('dy', '0.35em')
                        .style('font-size',   config.chartLabelFontSize + 'px')
                        .style('fill',        '#111')
                        .style('font-weight', 'bold')
                        .text(`${d3.format('.1e')(avg)} ± ${d3.format('.1e')(std)}`);
                })
                .on('mouseout', function () {
                    d3.select(this).attr('fill', barColor);
                    barsGroup.selectAll('.value-label').remove();
                });

            // Error bar
            if (std > 0) {
                const cy       = i * cfg.barHeight + (cfg.barHeight - 2) / 2;
                const capH     = Math.min(cfg.barHeight - 4, 6);
                const errEnd   = xScale(avg + std);
                const errStart = xScale(Math.max(avg - std, 0));
                // horizontal line
                barsGroup.append('line')
                    .attr('x1', errStart).attr('x2', errEnd)
                    .attr('y1', cy).attr('y2', cy)
                    .attr('stroke', '#333').attr('stroke-width', 1.5)
                    .style('pointer-events', 'none');
                // left cap
                barsGroup.append('line')
                    .attr('x1', errStart).attr('x2', errStart)
                    .attr('y1', cy - capH / 2).attr('y2', cy + capH / 2)
                    .attr('stroke', '#333').attr('stroke-width', 1.5)
                    .style('pointer-events', 'none');
                // right cap
                barsGroup.append('line')
                    .attr('x1', errEnd).attr('x2', errEnd)
                    .attr('y1', cy - capH / 2).attr('y2', cy + capH / 2)
                    .attr('stroke', '#333').attr('stroke-width', 1.5)
                    .style('pointer-events', 'none');
            }
        });

        return labelH + chartH + bc.barHeight; // panel height incl. axis space
    }

    /**
     * Render all metabolite bar chart panels for a node.
     * graph_info is now a list: [{metabolite_name, conditions: {...}}, ...]
     * One panel per list entry; entries with all-NaN conditions are silently skipped.
     */
    _renderMetaboliteBarChart(parentNode, element, data, config, options = {}) {
        const bc = config.barChart;
        const cfgOverride = {
            chartWidth:  bc.width,
            barHeight:   bc.barHeight,
            axisPadding: bc.axisPadding,
            barColor:    this._originColour(data.origin || 'unknown', config),
            hoverColor:  '#1f7a67',
            ...options,
        };
        delete cfgOverride.getPosition; // handled below

        const gi = data.graph_info;
        // New format: list of row-dicts
        // Old dict format still accepted for legacy JSON (silently wrapped)
        let rowList;
        if (Array.isArray(gi)) {
            rowList = gi;
        } else if (gi && typeof gi === 'object' && Object.keys(gi).length) {
            // Legacy dict: wrap in a single-element list
            rowList = [{ metabolite_name: '', conditions: gi }];
        } else {
            return;
        }

        const minCount = config.barMinCount ?? 1;

        // Filter out rows with no renderable conditions (any valid average)
        const renderable = rowList.filter(row => {
            const conds = row.conditions || {};
            return Object.values(conds).some(v => {
                const avg = v?.average;
                return avg !== null && avg !== undefined && isFinite(avg);
            });
        });
        if (!renderable.length) return;

        // Compute total height across all panels (for background rect and positioning)
        const barHeight   = cfgOverride.barHeight;
        const chartWidth  = cfgOverride.chartWidth;
        const axisPadding = cfgOverride.axisPadding;
        const gap         = barHeight; // gap between panels

        // Helper: compute panel height for a row given a condition-filter predicate
        const _panelH = (row, condPred) => {
            const conds = Object.keys(row.conditions || {}).filter(condPred);
            if (!conds.length) return 0;
            const labelH = row.metabolite_name ? barHeight : 0;
            return labelH + barHeight * conds.length + barHeight; // label + bars + axis
        };

        // "Full" height — all rows that have any valid average (ignoring count threshold)
        // Used to anchor the bottom of the chart at the same Y regardless of filtering.
        const fullPanelHeights = renderable.map(row =>
            _panelH(row, k => {
                const avg = row.conditions[k]?.average;
                return avg !== null && avg !== undefined && isFinite(avg);
            })
        );
        const fullTotalH = fullPanelHeights.reduce((s, h) => s + h, 0)
            + gap * Math.max(renderable.length - 1, 0);

        // "Visible" rows — those that have at least one bar passing the count threshold
        const visibleRows = renderable.filter(row =>
            Object.entries(row.conditions || {}).some(([k, v]) => {
                const avg = v?.average;
                if (avg === null || avg === undefined || !isFinite(avg)) return false;
                const cnt = v?.count ?? v?.n ?? Infinity;
                return cnt >= minCount;
            })
        );
        if (!visibleRows.length) return;

        // Visible panel heights (only conditions passing count threshold)
        const visiblePanelHeights = visibleRows.map(row =>
            _panelH(row, k => {
                const entry = row.conditions[k];
                const avg   = entry?.average;
                if (avg === null || avg === undefined || !isFinite(avg)) return false;
                const cnt = entry?.count ?? entry?.n ?? Infinity;
                return cnt >= minCount;
            })
        );
        const visibleTotalH = visiblePanelHeights.reduce((s, h) => s + h, 0)
            + gap * Math.max(visibleRows.length - 1, 0);

        // Y-shift so the bottom of the visible chart aligns with the full chart's bottom
        const yShift = fullTotalH - visibleTotalH;

        const getPosition = options.getPosition || ((d, c) => ({
            x: (d.x || 0) + config.barChartOffsetX,
            y: (d.y || 0) + config.barChartOffsetY - fullTotalH,
        }));

        const pos = getPosition(data, cfgOverride);

        const group = parentNode.append('g')
            .attr('class',     'bar-chart metabolite-bar-chart')
            .attr('transform', `translate(${pos.x},${pos.y})`);

        if (data.bigg_id) group.node().dataset.biggId = data.bigg_id;

        // Background rect — sized to visible content only, shifted down to bottom-align
        group.insert('rect', ':first-child')
            .attr('x',      -axisPadding - 4)
            .attr('y',      yShift - 8)
            .attr('width',   chartWidth + 8)
            .attr('height',  visibleTotalH + 16)
            .attr('fill',   'white')
            .attr('opacity', 0.88)
            .attr('rx', 4)
            .style('cursor', 'default');

        // Render visible panels (shifted down by yShift to keep bottom aligned)
        let yOff = yShift;
        visibleRows.forEach((row, ri) => {
            if (ri > 0) {
                // Separator line between panels
                group.append('line')
                    .attr('x1', -axisPadding)
                    .attr('x2', chartWidth - axisPadding + 4)
                    .attr('y1', yOff - gap / 2)
                    .attr('y2', yOff - gap / 2)
                    .attr('stroke', '#ddd')
                    .attr('stroke-width', 1);
            }
            const ph = this._renderMetaboliteBarChartPanel(
                group, row, data, config, cfgOverride, yOff
            );
            yOff += ph + gap;
        });

        // Keep position in sync with node movement
        const updatePos = () => {
            const np = getPosition(data, cfgOverride);
            group.attr('transform', `translate(${np.x},${np.y})`);
        };
        this._registerObserver(
            new MutationObserver(updatePos)
        ).observe(element.node(), { attributes: true, attributeFilter: ['transform'] });
    }

    // =========================================================================
    //  BAR CHARTS – PROTEOMICS
    // =========================================================================
    _renderProteomicsBarCharts(parentNode, segData, nodeData, config) {
        const bc       = config.barChart;
        const minCount = config.barMinCount ?? 1;
        const proteinList = segData.graph_info;
        if (!Array.isArray(proteinList) || !proteinList.length) return;

        // Filter conditions by minCount (use first protein's stats as reference keys)
        const allConditions = Object.keys(proteinList[0]?.stats || {});
        const conditions = allConditions.filter(c => {
            // A condition passes if at least one protein has count >= minCount for it
            return proteinList.some(p => {
                const entry = p.stats?.[c];
                const avg   = entry?.average;
                if (avg === null || avg === undefined || !isFinite(avg)) return false;
                const cnt = entry?.count ?? entry?.n ?? Infinity;
                return cnt >= minCount;
            });
        });
        if (!conditions.length) return;

        // Filter proteins: only keep those with at least one renderable condition
        const renderableProteins = proteinList.filter(p =>
            conditions.some(c => {
                const entry = p.stats?.[c];
                const avg   = entry?.average;
                if (avg === null || avg === undefined || !isFinite(avg)) return false;
                const cnt = entry?.count ?? entry?.n ?? Infinity;
                return cnt >= minCount;
            })
        );
        if (!renderableProteins.length) return;

        const nProteins = renderableProteins.length;
        const colours   = this._proteinColours(nProteins);
        const barH      = bc.barHeight;
        const axisPad   = bc.axisPadding;
        const chartW    = bc.width;
        const drawW     = chartW - axisPad;
        const singleH   = barH * conditions.length;
        const gap       = bc.gapBetween;

        let maxVal = 1e-9;
        renderableProteins.forEach(p =>
            conditions.forEach(c => {
                const avg = p.stats?.[c]?.average ?? 0;
                const std = p.stats?.[c]?.std_dev ?? 0;
                if (avg + std > maxVal) maxVal = avg + std;
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

        renderableProteins.forEach((prot, pi) => {
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
                    .call(d3.axisBottom(xScale).ticks(3).tickFormat(d => d3.format('.1e')(d).replace(/\.0(e)/, '$1')))
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
                const entry = prot.stats?.[cond];
                const avg   = entry?.average ?? null;
                const std   = entry?.std_dev ?? 0;
                const cnt   = entry?.count ?? entry?.n ?? Infinity;
                // Skip this bar if count is below threshold
                if (avg === null || cnt < minCount) return;
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

                // Error bar
                if (std > 0) {
                    const cy       = ci * barH + (barH - 2) / 2;
                    const capH     = Math.min(barH - 4, 6);
                    const errEnd   = xScale(avg + std);
                    const errStart = xScale(Math.max(avg - std, 0));
                    pGroup.append('line')
                        .attr('x1', errStart).attr('x2', errEnd)
                        .attr('y1', cy).attr('y2', cy)
                        .attr('stroke', '#333').attr('stroke-width', 1.5)
                        .style('pointer-events', 'none');
                    pGroup.append('line')
                        .attr('x1', errStart).attr('x2', errStart)
                        .attr('y1', cy - capH / 2).attr('y2', cy + capH / 2)
                        .attr('stroke', '#333').attr('stroke-width', 1.5)
                        .style('pointer-events', 'none');
                    pGroup.append('line')
                        .attr('x1', errEnd).attr('x2', errEnd)
                        .attr('y1', cy - capH / 2).attr('y2', cy + capH / 2)
                        .attr('stroke', '#333').attr('stroke-width', 1.5)
                        .style('pointer-events', 'none');
                }
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
                if (!data?.graph_info) return;
                // Accept both list format (new) and dict format (legacy)
                const gi = data.graph_info;
                const hasData = Array.isArray(gi)
                    ? gi.length > 0
                    : (gi && typeof gi === 'object' && Object.keys(gi).length > 0);
                if (!hasData) return;

                const element    = d3.select(nodes[i]);
                const parentNode = d3.select(nodes[i].parentNode);
                this._renderMetaboliteBarChart(
                    parentNode, element, data, config
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

    // =========================================================================
    //  CANVAS BAR CHARTS  (Vega-Lite embedded via <foreignObject> on the SVG)
    // =========================================================================
    /**
     * Render Vega-Lite bar charts directly on the Escher SVG canvas next to
     * every node/midpoint that has data in barchart_data.
     *
     * Charts are placed in <foreignObject> elements so they use the same
     * Vega-Lite spec as the sidebar panel.
     *
     * @param {object} jsonData  - Escher map JSON
     * @param {object} config    - vis config
     * @param {object} bcd       - { metabolites, reactions, kegg_names }
     */
    createCanvasBarCharts(jsonData, config, bcd) {
        if (!bcd || (!Object.keys(bcd.metabolites || {}).length && !Object.keys(bcd.reactions || {}).length)) return;

        const metByKegg = bcd.metabolites || {};
        const protByRxn = bcd.reactions   || {};
        const keggNames = bcd.kegg_names  || {};

        const svgEl = document.querySelector('#' + this.containerId + ' svg');
        if (!svgEl) return;

        // Remove any previously rendered canvas charts
        d3.select('#' + this.containerId).selectAll('.canvas-vega-chart').remove();

        // The Escher JSON nodes object — same data that Escher binds as D3 datums.
        // We use it as a fallback to look up reaction_kegg_ids by node ID when the
        // D3 datum on a midpoint circle doesn't carry that field directly.
        const jsonNodes = (jsonData[1] && jsonData[1].nodes) || {};

        const CHART_W  = 360;   // foreignObject width — must fit full chart incl. y-axis
        const CHART_H  = 180;
        const OFFSET_X = 20;
        const OFFSET_Y = -(CHART_H + 20);

        // Append foreignObjects into the Escher canvas group (the panned/zoomed <g>)
        // so that node coordinates are already in the right space.
        // Escher wraps everything in a <g class="escher-3d"> or similar; we find the
        // deepest <g> that contains the node circles.
        // Find the Escher canvas group (panned/zoomed <g>).
        // Walk up from a node circle to find the group that is a direct child of the SVG —
        // that is the zoom/pan group whose coordinate space matches the node transforms.
        // Find the Escher canvas group (panned/zoomed <g>).
        // Walk up from a node circle to find the group that is a direct child of the SVG.
        const firstCircle = d3.select('#' + this.containerId + ' svg .node-circle').node();
        let canvasG = null;
        if (firstCircle) {
            // Log the parent chain for debugging
            const chain = [];
            let dbg = firstCircle;
            while (dbg && dbg !== document.body) {
                chain.push(dbg.tagName + (dbg.className?.baseVal ? '.' + dbg.className.baseVal.split(' ').join('.') : ''));
                dbg = dbg.parentNode;
            }
            console.log('[EscherVisualizer] Node circle parent chain:', chain.join(' > '));

            // Walk up until we find a <g> whose parent is the SVG
            let el = firstCircle.parentNode;
            while (el && el !== svgEl) {
                if (el.parentNode === svgEl && el.tagName === 'g') {
                    canvasG = el;
                    break;
                }
                el = el.parentNode;
            }
        }
        // Fallback: try known Escher class names
        if (!canvasG) {
            canvasG = svgEl.querySelector('g.escher-3d')
                || svgEl.querySelector('g.zoom-g')
                || svgEl.querySelector('g.canvas-g')
                || svgEl.querySelector('g');
        }
        console.log('[EscherVisualizer] Canvas group found:', canvasG?.tagName,
            canvasG?.className?.baseVal || canvasG?.className);

        // We append foreignObjects into the canvas group so coordinates match.
        // If we couldn't find a canvas group, fall back to a dedicated overlay <g>
        // that we manually keep in sync with the canvas transform.
        let canvasSel;
        if (canvasG) {
            canvasSel = d3.select(canvasG);
        } else {
            // Last resort: create an overlay group on the SVG and copy the transform
            // from the first <g> child of the SVG (which is the zoom group).
            const zoomG = svgEl.querySelector('g');
            const overlayG = d3.select(svgEl).append('g')
                .attr('class', 'canvas-chart-overlay');
            if (zoomG) {
                const copyTransform = () => {
                    const tf = d3.select(zoomG).attr('transform');
                    if (tf) overlayG.attr('transform', tf);
                };
                copyTransform();
                new MutationObserver(copyTransform)
                    .observe(zoomG, { attributes: true, attributeFilter: ['transform'] });
            }
            canvasSel = overlayG;
        }

        // Count how many node circles exist
        const circleCount = d3.select('#' + this.containerId).selectAll('.node-circle').size();
        const metCircleCount = d3.select('#' + this.containerId).selectAll('.node-circle.metabolite-circle').size();
        console.log('[EscherVisualizer] Node circles:', circleCount, '| metabolite circles:', metCircleCount);
        console.log('[EscherVisualizer] metByKegg keys:', Object.keys(metByKegg).length,
            '| protByRxn keys:', Object.keys(protByRxn).length);

        // Helper: append a foreignObject Vega chart at (nx, ny) in canvas space
        const _appendChart = (nx, ny, safeId, spec) => {
            const fo = canvasSel.append('foreignObject')
                .attr('class', 'canvas-vega-chart')
                .attr('x', nx + OFFSET_X)
                .attr('y', ny + OFFSET_Y)
                .attr('width',  CHART_W)
                .attr('height', CHART_H)
                .style('overflow', 'visible')
                .style('pointer-events', 'none');

            fo.append('xhtml:div')
                .attr('id', safeId)
                .style('width',         CHART_W + 'px')
                .style('height',        CHART_H + 'px')
                .style('background',    'rgba(255,255,255,0.93)')
                .style('border',        '1px solid #ccc')
                .style('border-radius', '4px')
                .style('overflow',      'hidden')
                .style('pointer-events','auto');

            if (typeof vegaEmbed !== 'undefined') {
                // Strip the legend to save horizontal space on the canvas.
                // Deep-clone layers and remove color encoding legend.
                const stripLegend = layers => (layers || []).map(layer => {
                    if (!layer.encoding?.color) return layer;
                    return {
                        ...layer,
                        encoding: {
                            ...layer.encoding,
                            color: { ...layer.encoding.color, legend: null },
                        },
                    };
                });
                // For layered specs, autosize 'fit' doesn't work well.
                // Instead, set a small view width and let the y-axis labels
                // extend into the padding area. The foreignObject has overflow:visible.
                const compact = {
                    ...spec,
                    width:    160,   // plot area only; y-axis labels extend left into padding
                    height:   CHART_H - 50,
                    autosize: { type: 'none' },
                    title:    { ...(spec.title || {}), fontSize: 9 },
                    padding:  { left: 90, right: 8, top: 14, bottom: 20 },
                    layer:    stripLegend(spec.layer),
                    config: {
                        ...(spec.config || {}),
                        axis:   { labelFontSize: 7, titleFontSize: 0, labelLimit: 80 },
                        legend: { disable: true },
                        background: 'transparent',
                    },
                };
                // Use rAF to ensure foreignObject div is in the DOM before vegaEmbed
                requestAnimationFrame(() => {
                    vegaEmbed('#' + safeId, compact, { actions: false }).catch(console.error);
                });
            }
        };

        // ── Metabolite nodes ──────────────────────────────────────────────
        let metCount = 0;
        d3.select('#' + this.containerId)
            .selectAll('.node-circle.metabolite-circle')
            .each((d, i, nodes) => {
                if (!d) return;
                const keggId = d.bigg_id;
                const entries = metByKegg[keggId];
                if (!entries || !entries.length) return;

                const firstEntry = entries[0];
                const conds = (firstEntry.conditions || []).filter(
                    c => c.mean !== null && c.mean !== undefined && isFinite(c.mean)
                );
                if (!conds.length) return;

                const name  = keggNames[keggId] || d.name || keggId;
                const spec  = EscherVisualizer._buildVegaSpec(conds, name, 'Abundance');
                const safeId = 'cvega-met-' + keggId.replace(/[^a-zA-Z0-9]/g, '_');

                // Position: use parent <g> transform
                const parentG = nodes[i].parentNode;
                const tf = d3.select(parentG).attr('transform') || '';
                const m  = tf.match(/translate\(\s*([^,\s]+)[,\s]+([^)\s]+)\s*\)/);
                const nx = m ? parseFloat(m[1]) : (d.x || 0);
                const ny = m ? parseFloat(m[2]) : (d.y || 0);

                _appendChart(nx, ny, safeId, spec);
                metCount++;
            });

        // ── Midpoint (reaction) nodes ─────────────────────────────────────
        let rxnCount = 0;
        const CHART_GAP = 8;  // vertical gap between stacked protein charts
        d3.select('#' + this.containerId)
            .selectAll('.node-circle')
            .each((d, i, nodes) => {
                if (!d || d.node_type !== 'midpoint') return;

                // Try D3 datum first, then fall back to JSON lookup
                const nodeEl = nodes[i];
                const parentG = nodeEl.parentNode;

                // Get the reaction KEGG id from the midpoint's tooltip.
                // The graph builder stores it as d.tooltip.reaction_id (e.g. "R09293").
                // Fall back to jsonNodes lookup in case Escher strips the tooltip field.
                const rxnId = (d.tooltip && d.tooltip.reaction_id)
                    || (() => {
                        const gId = parentG.id || '';
                        const nid = gId.replace(/^n/, '');
                        const jn  = jsonNodes[nid];
                        return (jn && jn.tooltip && jn.tooltip.reaction_id) || null;
                    })();

                if (!rxnId) return;
                const protEntry = protByRxn[rxnId];
                if (!protEntry) return;

                const proteins = protEntry.proteins || [];
                if (!proteins.length) return;

                const tf = d3.select(parentG).attr('transform') || '';
                const m  = tf.match(/translate\(\s*([^,\s]+)[,\s]+([^)\s]+)\s*\)/);
                const nx = m ? parseFloat(m[1]) : (d.x || 0);
                const ny = m ? parseFloat(m[2]) : (d.y || 0);

                // Stack one chart per protein vertically.
                // The first chart sits at OFFSET_Y above the node; subsequent
                // charts are placed below it with a small gap.
                let stackOffset = 0;
                proteins.forEach((prot, pi) => {
                    const conds = (prot.conditions || []).filter(
                        c => c.mean !== null && c.mean !== undefined && isFinite(c.mean)
                    );
                    if (!conds.length) return;

                    const desc  = prot.description || '';
                    const title = desc
                        ? rxnId + ': ' + desc.slice(0, 35) + (desc.length > 35 ? '…' : '')
                        : rxnId;
                    const spec   = EscherVisualizer._buildVegaSpec(conds, title, 'Abundance');
                    const safeId = ('cvega-rxn-' + rxnId + '_p' + pi)
                        .replace(/[^a-zA-Z0-9]/g, '_');

                    // Compute per-chart height from the spec so stacking is accurate
                    const chartH = (spec.height || CHART_H) + 50; // +50 for title + padding

                    const fo = canvasSel.append('foreignObject')
                        .attr('class', 'canvas-vega-chart')
                        .attr('x', nx + OFFSET_X)
                        .attr('y', ny + OFFSET_Y - stackOffset)
                        .attr('width',  CHART_W)
                        .attr('height', chartH)
                        .style('overflow', 'visible')
                        .style('pointer-events', 'none');

                    fo.append('xhtml:div')
                        .attr('id', safeId)
                        .style('width',         CHART_W + 'px')
                        .style('height',        chartH + 'px')
                        .style('background',    'rgba(255,255,255,0.93)')
                        .style('border',        '1px solid #ccc')
                        .style('border-radius', '4px')
                        .style('overflow',      'hidden')
                        .style('pointer-events','auto');

                    if (typeof vegaEmbed !== 'undefined') {
                        const stripLegend = layers => (layers || []).map(layer => {
                            if (!layer.encoding?.color) return layer;
                            return {
                                ...layer,
                                encoding: {
                                    ...layer.encoding,
                                    color: { ...layer.encoding.color, legend: null },
                                },
                            };
                        });
                        const compact = {
                            ...spec,
                            width:    160,
                            height:   chartH - 50,
                            autosize: { type: 'none' },
                            title:    { ...(spec.title || {}), fontSize: 9 },
                            padding:  { left: 90, right: 8, top: 14, bottom: 20 },
                            layer:    stripLegend(spec.layer),
                            config: {
                                ...(spec.config || {}),
                                axis:   { labelFontSize: 7, titleFontSize: 0, labelLimit: 80 },
                                legend: { disable: true },
                                background: 'transparent',
                            },
                        };
                        requestAnimationFrame(() => {
                            vegaEmbed('#' + safeId, compact, { actions: false }).catch(console.error);
                        });
                    }

                    stackOffset += chartH + CHART_GAP;
                    rxnCount++;
                });
            });

        console.log(`[EscherVisualizer] Canvas charts: ${metCount} metabolite, ${rxnCount} reaction`);
    }

    // =========================================================================
    //  SIDEBAR CHART PANEL  (click-triggered, reference design with Vega-Lite)
    // =========================================================================
    /**
     * Show the reference-style sidebar panel for a clicked edge or metabolite node.
     *
     * For edges (reactions):
     *   - Header: rxnId + equation label
     *   - Proteomics section: one protein-block per protein with KO, description,
     *     entry_name, UniProt, KEGG gene, and a Vega-Lite bar chart
     *   - Metabolomics section: one sub-section per C-number in the equation
     *
     * For metabolite nodes (no rxnId):
     *   - Header: compound name
     *   - Metabolomics section only
     *
     * @param {string|null} rxnId         - KEGG reaction ID (e.g. "R00771") or null
     * @param {string}      equationLabel - full equation string or display label
     * @param {object|null} protEntry     - raw proteomics entry {reaction_id, proteins:[...]}
     * @param {string[]}    cNumbers      - C-numbers to show metabolomics for
     * @param {object}      metByKegg     - kegg_id -> [{kegg_id, metabolite, conditions:[...]}]
     * @param {object}      keggNames     - kegg_id -> human-readable name
     */
    showSidebarCharts(rxnId, equationLabel, protEntry, cNumbers, metByKegg, keggNames) {
        const panel       = document.getElementById('chart-panel');
        const titleEl     = document.getElementById('chart-panel-title');
        const subtitleEl  = document.getElementById('chart-panel-subtitle');
        const placeholder = document.getElementById('chart-panel-placeholder');
        const content     = document.getElementById('chart-panel-content');
        if (!panel || !content) return;

        // Update header
        if (titleEl)    titleEl.textContent    = rxnId || equationLabel || '—';
        if (subtitleEl) subtitleEl.textContent = equationLabel || '';

        // Show panel, hide placeholder, clear content
        if (placeholder) placeholder.style.display = 'none';
        content.style.display = 'block';
        content.innerHTML = '';

        const pendingCharts = [];  // {divId, conditions, title, yLabel}

        // ── Proteomics section ────────────────────────────────────────────
        const protSection = document.createElement('div');
        protSection.className = 'cp-chart-section';
        const protH3 = document.createElement('h3');
        protH3.textContent = rxnId ? `Proteomics — ${rxnId}` : 'Proteomics';
        protSection.appendChild(protH3);

        if (!protEntry) {
            const p = document.createElement('p');
            p.className = 'no-data';
            p.textContent = rxnId
                ? `No proteomics data for reaction ${rxnId}.`
                : 'No reaction ID found.';
            protSection.appendChild(p);
        } else {
            const proteins = protEntry.proteins || [];
            const rxnWrapper = document.createElement('div');
            rxnWrapper.style.marginBottom = '18px';
            const rxnLbl = document.createElement('div');
            rxnLbl.className = 'cp-reaction-label';
            rxnLbl.textContent = `${rxnId}  (${proteins.length} protein${proteins.length !== 1 ? 's' : ''})`;
            rxnWrapper.appendChild(rxnLbl);

            proteins.forEach((prot, pi) => {
                const block = document.createElement('div');
                block.className = 'cp-protein-block';

                const protLbl = document.createElement('div');
                protLbl.className = 'cp-protein-label';
                protLbl.textContent = prot.protein_id || prot.title || `Protein ${pi + 1}`;
                block.appendChild(protLbl);

                if (prot.ko) {
                    const d = document.createElement('div');
                    d.className = 'cp-meta-label';
                    d.textContent = 'KO: ' + prot.ko;
                    block.appendChild(d);
                }
                if (prot.description) {
                    const d = document.createElement('div');
                    d.className = 'cp-meta-label italic';
                    d.textContent = prot.description;
                    block.appendChild(d);
                }
                if (prot.entry_name) {
                    const d = document.createElement('div');
                    d.className = 'cp-meta-label';
                    d.textContent = 'Entry name: ' + prot.entry_name;
                    block.appendChild(d);
                }
                if (prot.entry) {
                    const d = document.createElement('div');
                    d.className = 'cp-meta-label';
                    d.textContent = 'UniProt: ' + prot.entry;
                    block.appendChild(d);
                }
                if (prot.kegg) {
                    const d = document.createElement('div');
                    d.className = 'cp-meta-label';
                    d.textContent = 'KEGG gene: ' + prot.kegg;
                    block.appendChild(d);
                }

                // Vega chart div
                const safeId = ('prot_' + (rxnId || 'x') + '_p' + pi).replace(/[^a-zA-Z0-9]/g, '_');
                const vegaDiv = document.createElement('div');
                vegaDiv.id = 'vega-' + safeId;
                block.appendChild(vegaDiv);
                rxnWrapper.appendChild(block);

                // Normalise conditions: support both {condition,mean,std,n} and {name,mean,std_dev,count}
                const rawConds = prot.conditions || [];
                pendingCharts.push({
                    divId:      '#vega-' + safeId,
                    conditions: rawConds,
                    title:      prot.protein_id || prot.title || `Protein ${pi + 1}`,
                    yLabel:     'Abundance',
                });
            });
            protSection.appendChild(rxnWrapper);
        }
        content.appendChild(protSection);

        // ── Metabolomics section ──────────────────────────────────────────
        const metSection = document.createElement('div');
        metSection.className = 'cp-chart-section';
        const metH3 = document.createElement('h3');
        metH3.textContent = cNumbers.length
            ? `Metabolomics — ${cNumbers.length} compound${cNumbers.length !== 1 ? 's' : ''} in equation`
            : 'Metabolomics';
        metSection.appendChild(metH3);

        if (!cNumbers.length) {
            const p = document.createElement('p');
            p.className = 'no-data';
            p.textContent = 'No C-numbers found.';
            metSection.appendChild(p);
        } else {
            let anyData = false;
            cNumbers.forEach(cId => {
                const entries = (metByKegg[cId] || []);
                if (!entries.length) return;
                anyData = true;

                const cHeader = document.createElement('div');
                cHeader.className = 'cp-compound-header';
                const name = (keggNames && keggNames[cId]) || cId;
                cHeader.textContent = name !== cId ? `${name} (${cId})` : cId;
                metSection.appendChild(cHeader);

                entries.forEach((metEntry, mi) => {
                    const rawConds = metEntry.conditions || [];
                    const validConds = rawConds.filter(c => {
                        const m = c.mean ?? c.mean;
                        return m !== null && m !== undefined && isFinite(m);
                    });
                    if (!validConds.length) return;

                    // Row label when multiple entries for same C-number
                    if (entries.length > 1 || metEntry.method) {
                        const lbl = document.createElement('div');
                        lbl.className = 'cp-row-label';
                        lbl.textContent = metEntry.metabolite || metEntry.name || name;
                        metSection.appendChild(lbl);
                        if (metEntry.method) {
                            const ml = document.createElement('div');
                            ml.className = 'cp-meta-label italic';
                            ml.textContent = 'Method: ' + metEntry.method;
                            metSection.appendChild(ml);
                        }
                    }

                    const safeId = 'met_' + cId.replace(/[^a-zA-Z0-9]/g, '_') + '_r' + mi;
                    const vegaDiv = document.createElement('div');
                    vegaDiv.id = safeId;
                    metSection.appendChild(vegaDiv);
                    pendingCharts.push({
                        divId:      '#' + safeId,
                        conditions: validConds,
                        title:      name !== cId ? `${name} (${cId})` : cId,
                        yLabel:     'Abundance',
                    });
                });
            });

            if (!anyData) {
                const p = document.createElement('p');
                p.className = 'no-data';
                p.textContent = 'No metabolomics data available for any compound in this reaction.';
                metSection.appendChild(p);
            }
        }
        content.appendChild(metSection);

        // Show panel before rendering Vega (needs DOM dimensions)
        panel.style.display = 'flex';

        // Render all Vega charts after DOM insertion
        pendingCharts.forEach(({ divId, conditions, title, yLabel }) => {
            const spec = EscherVisualizer._buildVegaSpec(conditions, title, yLabel);
            if (typeof vegaEmbed !== 'undefined') {
                vegaEmbed(divId, spec, { actions: false }).catch(console.error);
            }
        });
    }

    // =========================================================================
    //  VEGA-LITE SPEC BUILDER  (static helper)
    // =========================================================================
    /**
     * Build a Vega-Lite horizontal bar chart with error bars and individual dots.
     * Bars are coloured by n (number of replicates), matching the reference design.
     *
     * @param {Array}  conditions - [{condition|name, mean, std|std_dev, n|count, values?, null_columns?}]
     * @param {string} title
     * @param {string} yLabel
     */
    static _buildVegaSpec(conditions, title, yLabel) {
        const N_COLORS = ['#e53935', '#fb8c00', '#c0ca33', '#43a047', '#1e88e5', '#8e24aa'];

        // Normalise field names: support both raw (condition/std/n) and pre-processed (name/std_dev/count)
        const norm = c => ({
            condition:    c.condition    ?? c.name     ?? '',
            mean:         c.mean         ?? 0,
            std:          c.std          ?? c.std_dev  ?? null,
            n:            c.n            ?? c.count    ?? 1,
            values:       c.values       || [],
            null_columns: c.null_columns || [],
            columns:      c.columns      || [],
            subgroups:    c.subgroups    || [],   // [{subgroup, mean, std, n, values}, ...]
        });

        const normed = conditions.map(norm).filter(c => c.mean !== null && isFinite(c.mean));

        // ── Check whether any condition has subgroups ──────────────────────
        const hasSubgroups = normed.some(c => Array.isArray(c.subgroups) && c.subgroups.length > 0);

        if (hasSubgroups) {
            // ── Subgroup mode: faceted by condition, bioreps as grouped bars ─
            // Each condition gets its own facet row; bioreps are the y-axis within it.
            const sgValues = [];
            const sgDots   = [];
            const conditionOrder = [];

            normed.forEach(c => {
                conditionOrder.push(c.condition);
                if (Array.isArray(c.subgroups) && c.subgroups.length > 0) {
                    c.subgroups.forEach(sg => {
                        const hasReps = (sg.n || 1) > 1;
                        const stdVal  = hasReps && sg.std !== null && isFinite(sg.std) ? sg.std : null;
                        sgValues.push({
                            condition: c.condition,
                            subgroup:  sg.subgroup,
                            mean:      sg.mean,
                            std:       stdVal,
                            n:         sg.n,
                            lo:        hasReps ? sg.mean - (stdVal || 0) : null,
                            hi:        hasReps ? sg.mean + (stdVal || 0) : null,
                        });
                        (sg.values || []).forEach(v => {
                            if (v !== null && v !== undefined && isFinite(v))
                                sgDots.push({ condition: c.condition, subgroup: sg.subgroup, value: v });
                        });
                    });
                } else {
                    // condition has no subgroups — show as single bar
                    const hasReps = (c.n || 1) > 1;
                    const stdVal  = hasReps && c.std !== null && isFinite(c.std) ? c.std : null;
                    sgValues.push({
                        condition: c.condition, subgroup: c.condition,
                        mean: c.mean, std: stdVal, n: c.n,
                        lo: hasReps ? c.mean - (stdVal || 0) : null,
                        hi: hasReps ? c.mean + (stdVal || 0) : null,
                    });
                    (c.values || []).forEach(v => {
                        if (v !== null && v !== undefined && isFinite(v))
                            sgDots.push({ condition: c.condition, subgroup: c.condition, value: v });
                    });
                }
            });

            // Unique subgroup names for colour scale
            const sgNames      = [...new Set(sgValues.map(d => d.subgroup).filter(Boolean))];
            const sgColorRange = ['#1e88e5', '#43a047', '#fb8c00', '#e53935', '#8e24aa', '#00acc1'];
            const numSg        = sgNames.length || 1;
            const facetRowH    = numSg * 18 + 28;   // height per facet panel
            const chartH       = conditionOrder.length * facetRowH + 40;

            const colorEncoding = sgNames.length > 0 ? {
                field: 'subgroup', type: 'nominal',
                scale: { domain: sgNames, range: sgColorRange.slice(0, sgNames.length) },
                legend: { title: 'Biorep', labelFontSize: 9 },
            } : { value: '#1e88e5' };

            const tooltipFields = [
                { field: 'condition', type: 'nominal',      title: 'Condition' },
                { field: 'subgroup',  type: 'nominal',      title: 'Biorep' },
                { field: 'mean',      type: 'quantitative', title: 'Mean',  format: '.4g' },
                { field: 'std',       type: 'quantitative', title: 'Std',   format: '.4g' },
                { field: 'n',         type: 'quantitative', title: 'n' },
            ];

            return {
                $schema: 'https://vega.github.io/schema/vega-lite/v5.json',
                title:  { text: title, fontSize: 11, color: '#333' },
                width:  320,
                data:   { values: sgValues },
                config: {
                    axis:       { labelLimit: 120 },
                    view:       { stroke: '#ddd' },
                    mark:       { tooltip: true },
                    background: '#fafafa',
                    header:     { labelFontSize: 9, labelLimit: 160, titleFontSize: 0 },
                    facet:      { spacing: 4 },
                },
                facet: {
                    row: {
                        field: 'condition', type: 'nominal',
                        sort:  conditionOrder,
                        scale: { domain: conditionOrder },
                        header: { labelFontSize: 9, labelAngle: 0, labelAlign: 'left', labelLimit: 160, titleFontSize: 0 },
                    },
                },
                spec: {
                    height: facetRowH,
                    layer: [
                        {
                            mark: { type: 'bar', opacity: 0.85, cornerRadiusTopRight: 3, cornerRadiusBottomRight: 3 },
                            encoding: {
                                y: {
                                    field: 'subgroup', type: 'nominal',
                                    axis:  { labelFontSize: 8, title: null, labelPadding: 4 },
                                    sort:  sgNames,
                                },
                                x: {
                                    field: 'mean', type: 'quantitative',
                                    axis:  { title: yLabel || 'Abundance', titleFontSize: 9, labelFontSize: 8 },
                                },
                                color: colorEncoding,
                                tooltip: tooltipFields,
                            },
                        },
                        {
                            transform: [{ filter: 'datum.lo !== null && datum.hi !== null' }],
                            mark: { type: 'errorbar', color: '#555', ticks: true },
                            encoding: {
                                y:  { field: 'subgroup', type: 'nominal' },
                                x:  { field: 'lo',       type: 'quantitative' },
                                x2: { field: 'hi' },
                            },
                        },
                        ...(sgDots.length > 0 ? [{
                            data: { values: sgDots },
                            mark: { type: 'point', color: '#333', opacity: 0.6, size: 20, filled: true },
                            encoding: {
                                y: { field: 'subgroup', type: 'nominal' },
                                x: { field: 'value',    type: 'quantitative' },
                                tooltip: [
                                    { field: 'condition', type: 'nominal',      title: 'Condition' },
                                    { field: 'subgroup',  type: 'nominal',      title: 'Biorep' },
                                    { field: 'value',     type: 'quantitative', title: 'Value', format: '.4g' },
                                ],
                            },
                        }] : []),
                    ],
                },
            };
        }

        // ── Standard mode (no subgroups): one bar per condition ────────────
        const nDomain = [1, 2, 3, 4, 5, 6];
        
        // Build conditionOrder to preserve input order
        const conditionOrder = normed.map(c => c.condition);

        const values = normed.map(c => {
            const nCapped = Math.min(c.n || 1, 6);
            const hasReps = (c.n || 1) > 1;
            const stdVal  = hasReps && c.std !== null && isFinite(c.std) ? c.std : null;
            const cols    = (c.columns.length ? c.columns : c.null_columns).map(col => {
                if (col == null) return '(null)';
                return (c.null_columns || []).includes(col) ? col + ' (null)' : col;
            });
            return {
                condition:   c.condition,
                mean:        c.mean,
                std:         stdVal,
                n:           c.n,
                nCapped,
                lo:          hasReps ? c.mean - (stdVal || 0) : null,
                hi:          hasReps ? c.mean + (stdVal || 0) : null,
                columns:     cols.join('; '),
                indivValues: c.values.join(', '),
            };
        });

        const dotValues = [];
        normed.forEach(c => {
            (c.values || []).forEach(v => {
                if (v !== null && v !== undefined && isFinite(v)) {
                    dotValues.push({ condition: c.condition, value: v });
                }
            });
        });

        const barHeight   = Math.min(28, Math.max(16, Math.floor(160 / Math.max(values.length, 1))));
        const chartHeight = values.length * barHeight + 45;

        // Build a domain sort spec that preserves the input condition order
        // by explicitly specifying each condition value in order
        const conditionSortSpec = conditionOrder.length > 0 ? conditionOrder : undefined;

        return {
            $schema: 'https://vega.github.io/schema/vega-lite/v5.json',
            title:   { text: title, fontSize: 11, color: '#333' },
            width:   380,
            height:  chartHeight,
            config: {
                axis: { labelLimit: 0 },
                view: { stroke: 'transparent', continuousWidth: 300, continuousHeight: chartHeight },
                mark: { tooltip: true },
                background: '#fafafa',
            },
            layer: [
                {
                    data: { values },
                    mark: { type: 'bar', opacity: 0.85, cornerRadiusTopRight: 3, cornerRadiusBottomRight: 3, height: { band: 0.7 } },
                    encoding: {
                        y: {
                            field: 'condition', type: 'nominal',
                            axis:  { labelFontSize: 9, titleFontSize: 10, title: 'Condition', labelPadding: 8 },
                            sort:  conditionSortSpec,
                            scale: { domain: conditionOrder },
                        },
                        x: {
                            field: 'mean', type: 'quantitative',
                            axis:  { title: yLabel, titleFontSize: 10, labelFontSize: 9 },
                        },
                        color: {
                            field: 'nCapped', type: 'ordinal',
                            scale: { domain: nDomain, range: N_COLORS },
                            legend: { title: 'n (replicates)', labelExpr: "datum.label == '6' ? '≥6' : datum.label" },
                        },
                        tooltip: [
                            { field: 'condition',   type: 'nominal',      title: 'Condition' },
                            { field: 'mean',        type: 'quantitative', title: 'Mean',   format: '.4g' },
                            { field: 'std',         type: 'quantitative', title: 'Std',    format: '.4g' },
                            { field: 'n',           type: 'quantitative', title: 'n' },
                            { field: 'indivValues', type: 'nominal',      title: 'Individual values' },
                            { field: 'columns',     type: 'nominal',      title: 'Replicates' },
                        ],
                    },
                },
                {
                    data: { values: values.filter(v => v.lo !== null && v.hi !== null) },
                    mark: { type: 'errorbar', color: '#333', ticks: true },
                    encoding: {
                        y:  { field: 'condition', type: 'nominal' },
                        x:  { field: 'lo', type: 'quantitative' },
                        x2: { field: 'hi' },
                    },
                },
                ...(dotValues.length > 0 ? [{
                    data: { values: dotValues },
                    mark: { type: 'point', color: '#333', opacity: 0.7, size: 30, filled: true },
                    encoding: {
                        y: { field: 'condition', type: 'nominal' },
                        x: { field: 'value',     type: 'quantitative' },
                    },
                }] : []),
            ],
        };
    }
}
