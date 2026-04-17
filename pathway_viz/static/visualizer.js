/**
 * EscherVisualizer
 * Reads config from window.CONFIG (set by Flask template).
 *
 * Node colouring:
 *   Metabolite nodes are coloured by their 'origin' field:
 *     metabolomics → teal    (#2a9d8f)
 *     proteomics   → orange  (#e76f51)
 *     both         → purple  (#9b5de5)
 *     unknown      → grey    (#aaaaaa)
 *
 * Path highlighting (keep_positions=true):
 *   Nodes on the path get a thick coloured ring.
 *   Segments on the path get a thick solid coloured stroke.
 *   The rest of the graph is dimmed.
 *
 * Proteomics graph_info: [{protein_id, stats:{condition:{average,std_dev,count}}}]
 * Metabolomics graph_info: {condition:{average,std_dev,count}}
 */
class EscherVisualizer {

    // =========================================================================
    //  CONFIG
    // =========================================================================
    static get CONFIG() {
        const C = window.CONFIG || {};
        return {
            imageSize:               C.imageSize               ?? 200,
            defaultToBiggId:         false,
            nodeRadius:              C.nodeRadius              ?? 15,
            metaboliteRadius:        C.metaboliteRadius        ?? 10,
            reactionRadius:          C.reactionRadius          ?? 8,
            labelOffsetY:            C.labelOffsetY            ?? 20,
            coproductLabelOffsetY:   C.coproductLabelOffsetY   ?? 25,
            barChartOffsetY:         C.barChartOffsetY         ?? 60,
            metaboliteLabelFontSize: C.metaboliteLabelFontSize ?? 12,
            coproductLabelFontSize:  C.coproductLabelFontSize  ?? 11,
            chartTitleFontSize:      C.chartTitleFontSize      ?? 10,
            chartLabelFontSize:      C.chartLabelFontSize      ?? 9,
            barChartXLabel:          C.barChartXLabel          ?? '',
            originColours: Object.assign({
                metabolomics: '#2a9d8f',
                proteomics:   '#e76f51',
                both:         '#9b5de5',
                unknown:      '#aaaaaa',
            }, C.originColours || {}),
            // Path highlight style
            pathHighlight: {
                nodeStrokeWidth:    5,
                nodeStrokeColour:   '#f4d03f',
                segmentStrokeWidth: 6,
                segmentColour:      '#f4d03f',
                dimOpacity:         0.25,
            },
            barChart: {
                width:       C.barChartWidth       ?? 180,
                axisPadding: C.barChartAxisPadding ?? 55,
                barHeight:   C.barHeight           ?? 12,
                gapBetween:  C.barChartGap         ?? 18,
            },
        };
    }

    // =========================================================================
    //  COLOUR HELPERS
    // =========================================================================
    static _originColour(origin) {
        const map = EscherVisualizer.CONFIG.originColours;
        return map[origin] || map.unknown;
    }

    static _proteinColours(n) {
        const palette = [
            '#2a9d8f', '#e76f51', '#264653', '#e9c46a',
            '#a8dadc', '#f4a261', '#457b9d', '#e63946',
            '#06d6a0', '#118ab2',
        ];
        return Array.from({ length: n }, (_, i) => palette[i % palette.length]);
    }

    // =========================================================================
    //  ENTRY POINT
    // =========================================================================
    static initializeStructures(jsonData) {
        this.clearBarCharts();
        this._removeTooltip();
        this.stylePathwayElements();
        this.equalizeNodeRadii();
        this.colourNodesByOrigin();
        this.loadStructureImages(jsonData);
        this.createNodeBarCharts(jsonData);
        this.createSegmentBarCharts(jsonData[1].nodes, jsonData[1].reactions);
        this.initializeLabels(jsonData);
        this._attachTooltipListeners(jsonData);

        // Apply path highlight if the URL carries one
        const rawPath = (window.HIGHLIGHT_PATH || '').trim();
        if (rawPath) {
            const pathNodes = rawPath.split(',').map(s => s.trim()).filter(Boolean);
            if (pathNodes.length) {
                this.highlightPathInPlace(pathNodes, jsonData[1]);
            }
        }
    }

    static setupNaNLabelRemoval() {
        setInterval(() => this.removeStoichiometryLabels(), 500);
    }

    static removeStoichiometryLabels() {
        d3.selectAll('.stoichiometry-label').filter(function () {
            return d3.select(this).text() === 'NaN';
        }).remove();
    }

    // =========================================================================
    //  NODE COLOURING BY ORIGIN
    // =========================================================================
    static colourNodesByOrigin() {
        d3.select('#map_container')
            .selectAll('.node-circle.metabolite-circle')
            .each(function (d) {
                if (!d) return;
                const colour = EscherVisualizer._originColour(d.origin || 'unknown');
                d3.select(this)
                    .style('fill',         colour)
                    .style('fill-opacity',  0.85)
                    .style('stroke',       d3.color(colour).darker(0.6).toString())
                    .style('stroke-width', '1.5px');
            });
    }

    // =========================================================================
    //  PATH HIGHLIGHTING IN PLACE  (keep_positions=true)
    //
    //  - Dims everything not on the path.
    //  - Draws a thick bright ring around path nodes.
    //  - Draws a thick solid line along path segments.
    //  - Does NOT create a subgraph.
    // =========================================================================
    static highlightPathInPlace(pathNodeIds, mapData) {
        const C    = EscherVisualizer.CONFIG;
        const ph   = C.pathHighlight;
        const idSet = new Set(pathNodeIds);

        // Build set of segment IDs that connect consecutive path nodes
        // We need to match (from_node_id -> midpoint -> to_node_id)
        const nodes    = mapData.nodes    || {};
        const segments = {};
        Object.values(mapData.reactions || {}).forEach(rxn => {
            Object.values(rxn.segments || {}).forEach(seg => {
                segments[seg.from_node_id + '_' + seg.to_node_id] = seg;
            });
        });

        // For each consecutive pair in the path, find the midpoint(s)
        // between them and mark all segments in that chain.
        const pathSegKeys = new Set();

        // Walk nodes dict to find midpoints that connect path pairs
        Object.entries(nodes).forEach(([nid, nd]) => {
            if (nd.node_type !== 'midpoint') return;
            const fid = nd.from_node_id;
            const tid = nd.to_node_id;
            if (!fid || !tid) return;
            // This midpoint is on the path if both endpoints are in pathNodeIds
            // AND they are consecutive in the path array
            for (let i = 0; i < pathNodeIds.length - 1; i++) {
                const a = pathNodeIds[i];
                const b = pathNodeIds[i + 1];
                if ((fid === a && tid === b) || (fid === b && tid === a)) {
                    // Mark the two sub-segments
                    pathSegKeys.add(fid + '_' + nid);
                    pathSegKeys.add(nid + '_' + tid);
                    pathSegKeys.add(b   + '_' + nid);
                    pathSegKeys.add(nid + '_' + a);
                }
            }
        });

        // ── Dim / highlight nodes ─────────────────────────────────────
        d3.select('#map_container')
            .selectAll('.node-circle.metabolite-circle')
            .each(function (d) {
                if (!d) return;
                const onPath = idSet.has(d.bigg_id);
                const el     = d3.select(this);
                const colour = EscherVisualizer._originColour(d.origin || 'unknown');

                if (onPath) {
                    el
                        .style('fill',         colour)
                        .style('fill-opacity',  1.0)
                        .style('stroke',        ph.nodeStrokeColour)
                        .style('stroke-width',  ph.nodeStrokeWidth + 'px');
                } else {
                    el
                        .style('fill',         colour)
                        .style('fill-opacity',  ph.dimOpacity)
                        .style('stroke',        d3.color(colour).darker(0.6).toString())
                        .style('stroke-width',  '1.5px');
                }
            });

        // ── Dim / highlight segments ──────────────────────────────────
        d3.select('#map_container')
            .selectAll('path.segment')
            .each(function (d) {
                if (!d) return;
                const key    = (d.from_node_id || '') + '_' + (d.to_node_id || '');
                const onPath = pathSegKeys.has(key);
                const el     = d3.select(this);

                if (onPath) {
                    el
                        .style('stroke',            ph.segmentColour)
                        .style('stroke-width',       ph.segmentStrokeWidth + 'px')
                        .style('stroke-dasharray',  'none')
                        .style('opacity',            1.0);
                } else {
                    el
                        .style('stroke',           'black')
                        .style('stroke-width',     '1px')
                        .style('stroke-dasharray', '5,5')
                        .style('opacity',           ph.dimOpacity);
                }
            });

        // ── Add a CSS class so it can be cleared later ────────────────
        d3.select('#map_container')
            .classed('path-highlight-active', true);
    }

    static clearPathHighlight() {
        d3.select('#map_container')
            .classed('path-highlight-active', false);

        // Restore origin colours
        EscherVisualizer.colourNodesByOrigin();

        // Restore default segment style
        d3.select('#map_container')
            .selectAll('path.segment')
            .style('stroke',           'black')
            .style('stroke-width',     '3px')
            .style('stroke-dasharray', '5,5')
            .style('opacity',           1.0);
    }

    // =========================================================================
    //  TOOLTIP
    // =========================================================================
    static _getTooltipDiv() {
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

    static _removeTooltip() {
        const div = document.getElementById('escher-tooltip');
        if (div) div.style.display = 'none';
    }

    static _positionTooltip(div, evt) {
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

    // ── Formatters ───────────────────────────────────────────────────
    static _fmt(v) {
        if (v === null || v === undefined) return '—';
        const n = Number(v);
        if (isNaN(n)) return String(v);
        if (n === 0)  return '0';
        return (Math.abs(n) >= 1e5 || (Math.abs(n) < 1e-2 && n !== 0))
            ? n.toExponential(3)
            : n.toPrecision(5);
    }

    static _tr(label, value) {
        return `<tr>
          <td style="padding:2px 10px 2px 0;color:#666;
              white-space:nowrap;vertical-align:top">${label}</td>
          <td style="padding:2px 0;color:#111;
              word-break:break-word">${value}</td>
        </tr>`;
    }

    static _replicateTable(replicates) {
        if (!replicates || !replicates.length) return '—';
        const cells = replicates.map(v =>
            `<td style="padding:1px 4px;border:1px solid #e0e0e0;
                background:#f8f8f8;color:#333;font-size:10px;
                white-space:nowrap">${EscherVisualizer._fmt(v)}</td>`
        ).join('');
        return `<table style="border-collapse:collapse;
            display:inline-table"><tr>${cells}</tr></table>`;
    }

    static _conditionBlock(entry) {
        return `
          <tr><td colspan="2" style="padding-top:5px;padding-bottom:1px;
              font-weight:bold;color:#2a6496;font-size:11px">
            ${entry.name}
          </td></tr>
          ${EscherVisualizer._tr('Mean&nbsp;±&nbsp;SD',
              `${EscherVisualizer._fmt(entry.mean)}&nbsp;±&nbsp;`
              + `${EscherVisualizer._fmt(entry.std_dev)}`)}
          ${EscherVisualizer._tr('n', entry.count)}
          ${EscherVisualizer._tr('Replicates',
              EscherVisualizer._replicateTable(entry.replicates))}`;
    }

    static _divRow() {
        return `<tr><td colspan="2">
          <hr style="margin:3px 0;border:none;border-top:1px solid #eee">
        </td></tr>`;
    }

    static _metaboliteTooltipHTML(tt) {
        const origin  = tt.origin || 'unknown';
        const colour  = EscherVisualizer._originColour(origin);
        const header  = `
          <div style="border-bottom:2px solid ${colour};
              margin-bottom:6px;padding-bottom:5px;
              display:flex;align-items:baseline;gap:8px">
            <span style="font-size:13px;font-weight:bold;
                color:${colour};flex:1;word-break:break-word">
              ${tt.name || tt.id}
            </span>
            <span style="font-size:10px;color:#888;
                white-space:nowrap">${tt.id}</span>
            <span style="font-size:10px;background:#f0f0f0;
                border-radius:3px;padding:1px 5px;
                color:${colour};white-space:nowrap">${origin}</span>
          </div>`;
        if (!tt.conditions || !tt.conditions.length) {
            return header
                + '<div style="color:#999;font-style:italic">'
                + 'No abundance data</div>';
        }
        const rows = tt.conditions.map((c, i) =>
            (i > 0 ? EscherVisualizer._divRow() : '')
            + EscherVisualizer._conditionBlock(c)
        ).join('');
        return header
            + `<table style="border-collapse:collapse;width:100%;
               font-size:11px">${rows}</table>`;
    }

    static _reactionTooltipHTML(tt) {
        const rxnId    = tt.reaction_id || 'unknown';
        const fromName = tt.from_node?.name || tt.from_node?.id || '?';
        const toName   = tt.to_node?.name   || tt.to_node?.id   || '?';
        const fromId   = tt.from_node?.id   || '';
        const toId     = tt.to_node?.id     || '';
        const header   = `
          <div style="border-bottom:2px solid #e76f51;
              margin-bottom:6px;padding-bottom:5px">
            <span style="font-size:13px;font-weight:bold;
                color:#e76f51">${rxnId}</span>
          </div>
          <table style="border-collapse:collapse;width:100%;
              font-size:11px;margin-bottom:6px">
            ${EscherVisualizer._tr('From',
                `<strong>${fromName}</strong>
                 <span style="color:#aaa;font-size:10px">
                 &nbsp;(${fromId})</span>`)}
            ${EscherVisualizer._tr('To',
                `<strong>${toName}</strong>
                 <span style="color:#aaa;font-size:10px">
                 &nbsp;(${toId})</span>`)}
          </table>`;
        if (!tt.proteins || !tt.proteins.length) {
            return header
                + '<div style="color:#999;font-style:italic">'
                + 'No proteomics data</div>';
        }
        const nP      = tt.proteins.length;
        const colours = EscherVisualizer._proteinColours(nP);
        const blocks  = tt.proteins.map((prot, pi) => {
            const colour  = colours[pi];
            const shortId = prot.protein_id.length > 45
                ? prot.protein_id.slice(0, 43) + '…'
                : prot.protein_id;
            const condRows = (prot.conditions || []).map((c, i) =>
                (i > 0 ? EscherVisualizer._divRow() : '')
                + EscherVisualizer._conditionBlock(c)
            ).join('');
            return `
              <div style="margin-top:8px;border-left:3px solid ${colour};
                  padding-left:8px">
                <div title="${prot.protein_id}"
                    style="font-size:11px;font-weight:bold;
                    color:${colour};margin-bottom:3px;
                    word-break:break-all">${shortId}</div>
                <table style="border-collapse:collapse;
                    width:100%;font-size:11px">${condRows}</table>
              </div>`;
        }).join('');
        return header + blocks;
    }

    // ── Core hover helper ────────────────────────────────────────────
    static _attachHover(el, getHTML) {
        const self = EscherVisualizer;
        const div  = self._getTooltipDiv();
        el.addEventListener('mouseenter', evt => {
            div.innerHTML = getHTML();
            self._positionTooltip(div, evt);
        });
        el.addEventListener('mousemove', evt => {
            if (div.style.display !== 'none') {
                self._positionTooltip(div, evt);
            }
        });
        el.addEventListener('mouseleave', evt => {
            const related = evt.relatedTarget;
            if (related && (div.contains(related) || div === related)) return;
            div.style.display = 'none';
        });
    }

    static _attachTooltipListeners(jsonData) {
        const self = EscherVisualizer;

        const nodeTooltips = {};
        Object.entries(jsonData[1]?.nodes || {}).forEach(([nid, nd]) => {
            if (nd.tooltip) nodeTooltips[nid] = nd.tooltip;
        });

        const segTooltips = {};
        Object.values(jsonData[1]?.reactions || {}).forEach(rxn => {
            Object.values(rxn.segments || {}).forEach(seg => {
                if (seg.tooltip && seg.edge_type === 'reactant_edge') {
                    segTooltips[
                        (seg.from_node_id || '') + '_' + (seg.to_node_id || '')
                    ] = seg.tooltip;
                }
            });
        });

        // Metabolite node circles
        d3.select('#map_container')
            .selectAll('.node-circle.metabolite-circle')
            .each(function (d) {
                if (!d?.tooltip) return;
                const tt = d.tooltip;
                self._attachHover(this,
                    () => self._metaboliteTooltipHTML(tt));
            });

        // Midpoint node circles
        d3.select('#map_container')
            .selectAll('.node-circle')
            .each(function (d) {
                if (!d?.tooltip || d.node_type !== 'midpoint') return;
                const tt = d.tooltip;
                self._attachHover(this,
                    () => self._reactionTooltipHTML(tt));
            });

        // Segment paths
        d3.select('#map_container')
            .selectAll('path.segment')
            .each(function (d) {
                if (!d) return;
                const key = (d.from_node_id || '')
                    + '_' + (d.to_node_id || '');
                const tt = segTooltips[key];
                if (!tt) return;
                self._attachHover(this,
                    () => self._reactionTooltipHTML(tt));
            });

        // Bar chart groups (after render)
        requestAnimationFrame(() => {
            d3.select('#map_container')
                .selectAll('.metabolite-bar-chart')
                .each(function () {
                    const biggId = this.dataset?.biggId;
                    if (!biggId) return;
                    const tt = nodeTooltips[biggId];
                    if (!tt) return;
                    self._attachHover(this,
                        () => self._metaboliteTooltipHTML(tt));
                });

            d3.select('#map_container')
                .selectAll('.proteomics-bar-chart')
                .each(function () {
                    const key = this.dataset?.segKey;
                    if (!key) return;
                    const tt = segTooltips[key];
                    if (!tt) return;
                    self._attachHover(this,
                        () => self._reactionTooltipHTML(tt));
                });
        });
    }

    // =========================================================================
    //  LABELS
    // =========================================================================
    static addMetaboliteLabels() {
        d3.selectAll('g.node').each(function () {
            const circle = d3.select(this)
                .select('.node-circle.metabolite-circle');
            if (circle.empty()) return;
            const data = circle.data()[0];
            if (!data) return;

            const parentNode = d3.select(this);
            parentNode.selectAll('.label').style('display', 'none');

            let label = parentNode.select('.node-label.metabolite-name');
            if (label.empty()) {
                label = parentNode.append('text')
                    .attr('class',           'node-label metabolite-name')
                    .style('font-family',    'Arial, sans-serif')
                    .style('font-weight',    'bold')
                    .style('fill',           '#333')
                    .style('text-anchor',    'middle')
                    .style('pointer-events', 'none');
            }

            const update = () => {
                const C = EscherVisualizer.CONFIG;
                label
                    .style('font-size', C.metaboliteLabelFontSize + 'px')
                    .text((data.name || data.bigg_id || 'Unknown')
                          .replace(/;\s*$/, ''));
                const transform = circle.attr('transform');
                if (transform) {
                    const m = transform.match(/translate\(([^,]+),([^)]+)\)/);
                    if (m) {
                        label.attr('transform',
                            `translate(${m[1]},`
                            + `${parseFloat(m[2]) + C.labelOffsetY})`);
                    }
                }
            };
            update();
            new MutationObserver(update).observe(circle.node(), {
                attributes: true, attributeFilter: ['transform'],
            });
        });
    }

    static addCoproductLabels() {
        d3.selectAll('.coproduct-circle').each(function () {
            const circle = d3.select(this);
            const data   = circle.data()[0];
            if (!data) return;

            const parentNode = d3.select(this.parentNode);
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
                const C    = EscherVisualizer.CONFIG;
                const text = (C.defaultToBiggId
                    ? (data.bigg_id || data.name || 'Coproduct')
                    : (data.name    || data.bigg_id || 'Coproduct')
                ).replace(/;\s*$/, '');
                label
                    .style('font-size', C.coproductLabelFontSize + 'px')
                    .text(text);
                const transform = circle.attr('transform');
                if (transform) {
                    const m = transform.match(/translate\(([^,]+),([^)]+)\)/);
                    if (m) {
                        label.attr('transform',
                            `translate(${m[1]},`
                            + `${parseFloat(m[2]) - C.coproductLabelOffsetY})`);
                    }
                }
            };
            update();
            new MutationObserver(update).observe(circle.node(), {
                attributes: true, attributeFilter: ['transform'],
            });
        });
    }

    static initializeLabels(jsonData) {
        this.removeStoichiometryLabels();
        this.addMetaboliteLabels();
        this.addCoproductLabels();
    }

    // =========================================================================
    //  STRUCTURE IMAGES
    // =========================================================================
    static loadStructureImages(jsonData) {
        const C              = EscherVisualizer.CONFIG;
        const maxDisplaySize = C.imageSize;
        const isVertical     = window.CONFIG?.smallGraphLayoutVertical ?? true;
        const nodeThreshold  = window.CONFIG?.nodeThresholdSmall ?? 10;

        fetch('static/structure_imgs/image_dimensions.json')
            .then(r => r.ok ? r.json() : {})
            .catch(() => ({}))
            .then(dimensions => {
                const circles  = d3.select('#map_container')
                    .selectAll('.metabolite-circle');
                const numNodes = circles.size();
                const useLeft  = isVertical && numNodes < nodeThreshold;

                let maxNatW = Object.values(dimensions)
                    .reduce((m, d) => Math.max(m, d.w), 0);
                if (maxNatW === 0) maxNatW = 1;
                const scale = maxDisplaySize / maxNatW;

                circles.each(function (data) {
                    const circle     = d3.select(this);
                    const parentNode = d3.select(circle.node().parentNode);
                    if (data.highlight) circle.classed('highlighted', true);

                    const imgPath   = `static/structure_imgs/${data.bigg_id}.png`;
                    const dim       = dimensions[data.bigg_id];
                    const dispW     = dim ? dim.w * scale : maxDisplaySize;
                    const dispH     = dim ? dim.h * scale : maxDisplaySize;
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
                    img.onerror = () => {};
                    img.src = imgPath;

                    const updatePos = () => {
                        const transform = circle.attr('transform');
                        if (!transform) return;
                        const o = getOffset();
                        parentNode.select('image').attr('transform',
                            `${transform} translate(${o.x},${o.y})`);
                    };
                    updatePos();
                    new MutationObserver(updatePos).observe(circle.node(), {
                        attributes: true, attributeFilter: ['transform'],
                    });
                });
            });
    }

    // =========================================================================
    //  BAR CHARTS – METABOLOMICS
    // =========================================================================
    static _renderMetaboliteBarChart(parentNode, element, data, options = {}) {
        const C  = EscherVisualizer.CONFIG;
        const bc = C.barChart;
        const cfg = Object.assign({
            chartWidth:  bc.width,
            barHeight:   bc.barHeight,
            axisPadding: bc.axisPadding,
            barColor:    EscherVisualizer._originColour(data.origin || 'unknown'),
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
            : {
                x: (data.x || 0) - cfg.chartWidth,
                y: (data.y || 0) - chartH / 2,
            };

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
        new MutationObserver(updatePos).observe(element.node(), {
            attributes: true, attributeFilter: ['transform'],
        });

        group.insert('rect', ':first-child')
            .attr('class',   'chart-bg')
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
                .style('font-size',   C.chartLabelFontSize + 'px')
                .style('fill',        '#555')
                .text(cond);
        });

        group.append('g')
            .attr('class',     'x-axis')
            .attr('transform', `translate(0,${chartH})`)
            .call(d3.axisBottom(xScale).ticks(3).tickFormat(d3.format('.1e')))
            .selectAll('text')
            .style('font-size', C.chartLabelFontSize + 'px')
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
                        .style('font-size',   C.chartLabelFontSize + 'px')
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
    //  BAR CHARTS – PROTEOMICS  (one chart per protein, stacked)
    // =========================================================================
    static _renderProteomicsBarCharts(parentNode, segData, nodeData) {
        const C  = EscherVisualizer.CONFIG;
        const bc = C.barChart;

        const proteinList = segData.graph_info;
        if (!Array.isArray(proteinList) || proteinList.length === 0) return;

        const conditions = Object.keys(proteinList[0]?.stats || {});
        if (!conditions.length) return;

        const nProteins = proteinList.length;
        const colours   = EscherVisualizer._proteinColours(nProteins);
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

        const segKey = (segData.from_node_id || '')
            + '_' + (segData.to_node_id || '');

        const outerGroup = parentNode.append('g')
            .attr('class',     'bar-chart proteomics-bar-chart')
            .attr('transform', `translate(${anchorX},${anchorY})`);

        outerGroup.node().dataset.segKey = segKey;

        const totalH = nProteins * singleH + (nProteins - 1) * gap;

        outerGroup.insert('rect', ':first-child')
            .attr('class',   'chart-bg')
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

            const pGroup = outerGroup.append('g')
                .attr('transform', `translate(0,${offsetY})`);

            const shortId = prot.protein_id.length > 28
                ? prot.protein_id.slice(0, 26) + '…'
                : prot.protein_id;
            pGroup.append('text')
                .attr('x',  drawW / 2)
                .attr('y',  -4)
                .style('text-anchor', 'middle')
                .style('font-size',   C.chartTitleFontSize + 'px')
                .style('font-weight', 'bold')
                .style('fill',        colour)
                .text(shortId)
                .append('title').text(prot.protein_id);

            if (pi > 0) {
                pGroup.append('line')
                    .attr('x1', -axisPad)
                    .attr('x2',  chartW - axisPad + 4)
                    .attr('y1', -gap / 2 - singleH)
                    .attr('y2', -gap / 2 - singleH)
                    .attr('stroke',       '#ddd')
                    .attr('stroke-width',  1);
            }

            conditions.forEach((cond, ci) => {
                pGroup.append('text')
                    .attr('x',  -4)
                    .attr('y',   ci * barH + barH / 2)
                    .attr('dy', '0.35em')
                    .style('text-anchor', 'end')
                    .style('font-size',   C.chartLabelFontSize + 'px')
                    .style('fill',        '#555')
                    .text(cond);
            });

            if (pi === nProteins - 1) {
                pGroup.append('g')
                    .attr('class',     'x-axis')
                    .attr('transform', `translate(0,${singleH})`)
                    .call(
                        d3.axisBottom(xScale).ticks(3).tickFormat(d3.format('.1e'))
                    )
                    .selectAll('text')
                    .style('font-size', C.chartLabelFontSize + 'px')
                    .style('fill',      '#555');

                if (C.barChartXLabel) {
                    pGroup.append('text')
                        .attr('x',  drawW / 2)
                        .attr('y',  singleH + 28)
                        .style('text-anchor', 'middle')
                        .style('font-size',   C.chartLabelFontSize + 'px')
                        .style('fill',        '#666')
                        .text(C.barChartXLabel);
                }
            }

            conditions.forEach((cond, ci) => {
                const avg = prot.stats?.[cond]?.average ?? null;
                const std = prot.stats?.[cond]?.std_dev ?? 0;
                if (avg === null) return;
                const barW = xScale(avg);
                pGroup.append('rect')
                    .attr('class',  'bar')
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
                            .style('font-size',   C.chartLabelFontSize + 'px')
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
    //  NODE BAR CHARTS  (metabolomics)
    // =========================================================================
    static createNodeBarCharts(jsonData) {
        const C         = EscherVisualizer.CONFIG;
        const isVert    = window.CONFIG?.smallGraphLayoutVertical ?? true;
        const threshold = window.CONFIG?.nodeThresholdSmall ?? 10;
        const numNodes  = d3.select('#map_container')
            .selectAll('.metabolite-circle').size();
        const useLeft   = isVert && numNodes < threshold;

        d3.select('#map_container')
            .selectAll('.node-circle.metabolite-circle')
            .each(function (data) {
                if (!data?.graph_info || Array.isArray(data.graph_info)) return;
                const element    = d3.select(this);
                const parentNode = d3.select(element.node().parentNode);
                EscherVisualizer._renderMetaboliteBarChart(
                    parentNode, element, data, {
                        getPosition: (d, cfg) => useLeft
                            ? { x: d.x + C.barChartOffsetY,
                                y: d.y - cfg.chartWidth }
                            : { x: d.x - C.barChart.width / 2,
                                y: d.y + C.barChartOffsetY },
                    }
                );
            });
    }

    // =========================================================================
    //  SEGMENT BAR CHARTS  (proteomics – stacked per protein)
    // =========================================================================
    static createSegmentBarCharts(nodeData, reactions) {
        const segments = {};
        Object.values(reactions || {}).forEach(rxn => {
            Object.values(rxn.segments || {}).forEach(seg => {
                if (
                    seg.edge_type === 'reactant_edge' &&
                    Array.isArray(seg.graph_info) &&
                    seg.graph_info.length > 0
                ) {
                    segments[
                        (seg.from_node_id || '') + '_' + (seg.to_node_id || '')
                    ] = seg;
                }
            });
        });

        d3.select('#map_container').selectAll('.segment')
            .each(function (data) {
                if (!data) return;
                const key = (data.from_node_id || '')
                    + '_' + (data.to_node_id || '');
                const seg = segments[key];
                if (!seg) return;
                const parentNode = d3.select(
                    d3.select(this).node().parentNode
                );
                EscherVisualizer._renderProteomicsBarCharts(
                    parentNode, seg, nodeData
                );
            });
    }

    // =========================================================================
    //  NODE STYLING
    // =========================================================================
    static stylePathwayElements() {
        d3.selectAll('.segment')
            .style('stroke',           'black')
            .style('stroke-width',     '3px')
            .style('stroke-dasharray', '5,5');
    }

    static equalizeNodeRadii() {
        const C     = EscherVisualizer.CONFIG;
        const metR  = C.metaboliteRadius ?? C.nodeRadius;
        const reacR = C.reactionRadius   ?? C.nodeRadius;
        const update = () => {
            d3.selectAll('circle.node-circle.metabolite-circle')
                .each(function () { this.setAttribute('r', metR); });
            d3.selectAll('circle.coproduct-circle')
                .each(function () { this.setAttribute('r', reacR); });
        };
        update();
        const container = document.getElementById('map_container');
        if (container && !this._radiusObserverActive) {
            new MutationObserver(() => {
                let needs = false;
                d3.selectAll('circle.node-circle.metabolite-circle')
                    .each(function () {
                        if (parseFloat(this.getAttribute('r')) !== metR)
                            needs = true;
                    });
                if (needs) update();
            }).observe(container, {
                subtree: true, attributes: true, attributeFilter: ['r'],
            });
            this._radiusObserverActive = true;
        }
    }

    static clearBarCharts() {
        d3.selectAll('.bar-chart').remove();
    }

    // =========================================================================
    //  LEGACY HIGHLIGHTING  (subgraph mode — kept for compatibility)
    // =========================================================================
    static highlightPath(pathNodes) {
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

    static highlightMultiNodes(selectedNodeIds) {
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