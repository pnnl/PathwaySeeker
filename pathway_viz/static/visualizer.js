/**
 * EscherVisualizer
 * Reads config from window.CONFIG (set by Flask template).
 *
 * Proteomics graph_info is now a LIST:
 *   [{protein_id: str, stats: {condition: {average, std_dev, count}}}]
 * → one grouped bar chart per segment, bars grouped by condition,
 *   one bar-group per protein, with a colour legend.
 *
 * Metabolomics graph_info remains a flat dict:
 *   {condition: {average, std_dev, count}}
 * → one bar chart per node, one bar per condition.
 */
class EscherVisualizer {

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
            coproductLabelFontSize:  C.coproductLabelFontSize  ?? 12,
            chartTitleFontSize:      C.chartTitleFontSize      ?? 12,
            chartLabelFontSize:      C.chartLabelFontSize      ?? 10,
            barChartTitle:           C.barChartTitle           ?? '',
            barChartXLabel:          C.barChartXLabel          ?? '',
            barChartYLabel:          C.barChartYLabel          ?? '',
            barChart: {
                width:       C.barChartWidth       ?? 160,
                height:      C.barChartHeight      ?? 120,
                axisPadding: C.barChartAxisPadding ?? 20,
                barHeight:   C.barHeight           ?? 14,
            },
        };
    }

    // ── Colour palette for proteins ──────────────────────────────────
    static _proteinColours(n) {
        const palette = [
            '#2a9d8f', '#e76f51', '#264653', '#e9c46a',
            '#a8dadc', '#f4a261', '#457b9d', '#e63946',
            '#06d6a0', '#118ab2',
        ];
        return Array.from({length: n}, (_, i) => palette[i % palette.length]);
    }

    // =========================================================================
    //  ENTRY POINT
    // =========================================================================
    static initializeStructures(jsonData) {
        this.clearBarCharts();
        this.stylePathwayElements();
        this.equalizeNodeRadii();
        this.loadStructureImages(jsonData);
        this.createNodeBarCharts(jsonData);
        this.createSegmentBarCharts(jsonData[1].nodes, jsonData[1].reactions);
        this.initializeLabels(jsonData);
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
    //  LABELS
    // =========================================================================
    static addMetaboliteLabels() {
        d3.selectAll('g.node').each(function () {
            const circle = d3.select(this).select('.node-circle.metabolite-circle');
            if (circle.empty()) return;
            const data = circle.data()[0];
            if (!data) return;
            const parentNode = d3.select(this);
            parentNode.selectAll('.label').style('display', 'none');
            let label = parentNode.select('.node-label.metabolite-name');
            if (label.empty()) {
                label = parentNode.append('text')
                    .attr('class', 'node-label metabolite-name')
                    .style('font-family', 'Arial, sans-serif')
                    .style('font-weight', 'bold')
                    .style('fill', '#333')
                    .style('text-anchor', 'middle')
                    .style('pointer-events', 'none');
            }
            const update = () => {
                const C = EscherVisualizer.CONFIG;
                label.style('font-size', C.metaboliteLabelFontSize + 'px');
                const transform = circle.attr('transform');
                if (transform) {
                    const m = transform.match(/translate\(([^,]+),([^)]+)\)/);
                    if (m) {
                        label.attr('transform',
                            `translate(${m[1]},${parseFloat(m[2]) + C.labelOffsetY})`);
                    }
                }
                label.text((data.name || data.bigg_id || 'Unknown').replace(/;\s*$/, ''));
            };
            update();
            new MutationObserver(update)
                .observe(circle.node(), {attributes: true, attributeFilter: ['transform']});
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
                    .attr('class', 'node-label coproduct-name')
                    .style('font-family', 'Arial, sans-serif')
                    .style('fill', '#666')
                    .style('text-anchor', 'middle')
                    .style('pointer-events', 'none');
            }
            const update = () => {
                const C = EscherVisualizer.CONFIG;
                label.style('font-size', C.coproductLabelFontSize + 'px');
                const transform = circle.attr('transform');
                if (transform) {
                    const m = transform.match(/translate\(([^,]+),([^)]+)\)/);
                    if (m) {
                        label.attr('transform',
                            `translate(${m[1]},${parseFloat(m[2]) - C.coproductLabelOffsetY})`);
                    }
                }
                const text = C.defaultToBiggId
                    ? (data.bigg_id || data.name || 'Coproduct')
                    : (data.name    || data.bigg_id || 'Coproduct');
                label.text(text.replace(/;\s*$/, ''));
            };
            update();
            new MutationObserver(update)
                .observe(circle.node(), {attributes: true, attributeFilter: ['transform']});
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
                const circles  = d3.select('#map_container').selectAll('.metabolite-circle');
                const numNodes = circles.size();
                const useLeft  = isVertical && numNodes < nodeThreshold;
                let maxNatW    = Object.values(dimensions).reduce((m, d) => Math.max(m, d.w), 0);
                if (maxNatW === 0) maxNatW = 1;
                const scale = maxDisplaySize / maxNatW;

                circles.each(function (data) {
                    const circle     = d3.select(this);
                    const parentNode = d3.select(circle.node().parentNode);
                    if (data.highlight) circle.classed('highlighted', true);
                    circle.style('fill', '#999').style('stroke', '#666').style('fill-opacity', 0.7);

                    const imgPath = `static/structure_imgs/${data.bigg_id}.png`;
                    const dim     = dimensions[data.bigg_id];
                    const dispW   = dim ? dim.w * scale : maxDisplaySize;
                    const dispH   = dim ? dim.h * scale : maxDisplaySize;
                    const getOffset = () => useLeft
                        ? {x: -dispW - 10, y: -dispH}
                        : {x: -dispW / 2,  y: -dispH - 10};

                    const img = new Image();
                    img.onload = () => {
                        const o = getOffset();
                        parentNode.insert('image', 'text')
                            .attr('class', 'structure-image')
                            .attr('transform', `translate(${data.x + o.x},${data.y + o.y})`)
                            .attr('width', dispW).attr('height', dispH)
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
                    new MutationObserver(updatePos)
                        .observe(circle.node(), {attributes: true, attributeFilter: ['transform']});
                });
            });
    }

    // =========================================================================
    //  BAR CHARTS – METABOLOMICS  (flat dict → one bar per condition)
    // =========================================================================
    static _renderMetaboliteBarChart(parentNode, element, data, options = {}) {
        const C   = EscherVisualizer.CONFIG;
        const bc  = C.barChart;
        const cfg = Object.assign({
            chartWidth:     bc.width,
            barHeight:      bc.barHeight,
            axisPadding:    bc.axisPadding,
            barColor:       '#2a9d8f',
            hoverColor:     '#1f7a67',
            positionOffset: {x: 0, y: 0},
            getPosition:    null,
        }, options);

        const gi = data.graph_info;
        if (!gi || typeof gi !== 'object' || Array.isArray(gi)) return;

        const conditions = Object.keys(gi).filter(k => k !== 'metabolite' && k !== 'KEGG_C_number');
        if (!conditions.length) return;

        const values  = conditions.map(k => gi[k]?.average ?? 0);
        const maxVal  = Math.max(...values, 1e-9);
        const chartH  = cfg.barHeight * conditions.length;

        const xScale = d3.scaleLinear()
            .domain([0, maxVal])
            .range([0, cfg.chartWidth - cfg.axisPadding]);

        const pos = cfg.getPosition
            ? cfg.getPosition(data, cfg)
            : {
                x: (data.x || 0) - cfg.chartWidth + cfg.positionOffset.x,
                y: (data.y || 0) - chartH / 2    + cfg.positionOffset.y,
            };

        const group = parentNode.append('g')
            .attr('class', 'bar-chart metabolite-bar-chart')
            .attr('transform', `translate(${pos.x},${pos.y})`);

        const updatePos = () => {
            const np = cfg.getPosition ? cfg.getPosition(data, cfg) : {
                x: (data.x || 0) - cfg.chartWidth + cfg.positionOffset.x,
                y: (data.y || 0) - chartH / 2    + cfg.positionOffset.y,
            };
            group.attr('transform', `translate(${np.x},${np.y})`);
        };
        new MutationObserver(updatePos)
            .observe(element.node(), {attributes: true, attributeFilter: ['transform']});

        // Background
        group.insert('rect', ':first-child')
            .attr('x', -cfg.axisPadding).attr('y', -10)
            .attr('width', cfg.chartWidth + 10).attr('height', chartH + 20)
            .attr('fill', 'white').attr('opacity', 0.85).attr('rx', 4);

        // Y-axis labels
        const yScale = d3.scalePoint()
            .domain(conditions)
            .range([cfg.barHeight / 2, chartH - cfg.barHeight / 2]);
        group.append('g').attr('class', 'y-axis')
            .call(d3.axisLeft(yScale).tickSize(0))
            .selectAll('text')
            .style('font-size', C.chartLabelFontSize + 'px')
            .style('fill', '#555').style('text-anchor', 'end')
            .attr('x', -4).style('dy', '0.35em');

        // X-axis
        group.append('g').attr('class', 'x-axis')
            .attr('transform', `translate(0,${chartH})`)
            .call(d3.axisBottom(xScale).ticks(3).tickFormat(d3.format('.2e')))
            .selectAll('text')
            .style('font-size', C.chartLabelFontSize + 'px').style('fill', '#555');

        // Bars
        group.selectAll('.bar')
            .data(conditions)
            .enter().append('rect')
            .attr('class', 'bar')
            .attr('x', 0)
            .attr('y', (d, i) => i * cfg.barHeight)
            .attr('width',  d => xScale(gi[d]?.average ?? 0))
            .attr('height', cfg.barHeight - 2)
            .attr('fill',   cfg.barColor)
            .on('mouseover', function (d) {
                d3.select(this).attr('fill', cfg.hoverColor);
                const avg = gi[d]?.average ?? 0;
                const std = gi[d]?.std_dev ?? 0;
                group.append('text').attr('class', 'value-label')
                    .attr('x', xScale(avg) + 3)
                    .attr('y', conditions.indexOf(d) * cfg.barHeight + cfg.barHeight / 2)
                    .attr('dy', '0.35em')
                    .style('font-size', C.chartLabelFontSize + 'px')
                    .style('fill', '#222').style('font-weight', 'bold')
                    .text(`${d3.format('.2e')(avg)} ± ${d3.format('.1e')(std)}`);
            })
            .on('mouseout', function (d) {
                d3.select(this).attr('fill', cfg.barColor);
                group.selectAll('.value-label').remove();
            });
    }

    // =========================================================================
    //  BAR CHARTS – PROTEOMICS  (list → grouped bars, one group per condition,
    //                             one bar per protein, colour legend)
    // =========================================================================
    static _renderProteomicsBarChart(parentNode, element, data, nodeData, options = {}) {
        const C  = EscherVisualizer.CONFIG;
        const bc = C.barChart;

        const proteinList = data.graph_info;
        if (!Array.isArray(proteinList) || proteinList.length === 0) return;

        // Collect all conditions in consistent order
        const conditionSet = new Set();
        proteinList.forEach(p => Object.keys(p.stats || {}).forEach(c => conditionSet.add(c)));
        const conditions = Array.from(conditionSet);
        if (!conditions.length) return;

        const nProteins   = proteinList.length;
        const colours     = EscherVisualizer._proteinColours(nProteins);
        const barH        = bc.barHeight;
        const groupH      = barH * nProteins + 2;          // height of one condition group
        const chartH      = groupH * conditions.length;
        const chartW      = bc.width;
        const axisPad     = bc.axisPadding;

        // Max value across all proteins and conditions
        let maxVal = 1e-9;
        proteinList.forEach(p =>
            conditions.forEach(c => {
                const v = p.stats?.[c]?.average ?? 0;
                if (v > maxVal) maxVal = v;
            })
        );

        const xScale = d3.scaleLinear()
            .domain([0, maxVal])
            .range([0, chartW - axisPad]);

        // Position near the segment midpoint
        const fromNode = nodeData[data.to_node_id];
        const px = (fromNode?.x ?? 0) - chartW - 10;
        const py = (fromNode?.y ?? 0) - chartH / 2;

        const group = parentNode.append('g')
            .attr('class', 'bar-chart proteomics-bar-chart')
            .attr('transform', `translate(${px},${py})`);

        // Background
        const legendH = nProteins * 14 + 6;
        group.insert('rect', ':first-child')
            .attr('x', -axisPad - 4).attr('y', -legendH - 14)
            .attr('width',  chartW + axisPad + 8)
            .attr('height', chartH + legendH + 24)
            .attr('fill', 'white').attr('opacity', 0.88).attr('rx', 4);

        // ── Legend ───────────────────────────────────────────────────
        const legend = group.append('g')
            .attr('class', 'proteomics-legend')
            .attr('transform', `translate(0,${-legendH - 4})`);

        proteinList.forEach((p, i) => {
            const row = legend.append('g')
                .attr('transform', `translate(0,${i * 14})`);
            row.append('rect')
                .attr('width', 10).attr('height', 10)
                .attr('fill', colours[i]);
            row.append('text')
                .attr('x', 14).attr('y', 9)
                .style('font-size', (C.chartLabelFontSize - 1) + 'px')
                .style('fill', '#333')
                .text(p.protein_id.length > 30
                    ? p.protein_id.slice(0, 28) + '…'
                    : p.protein_id);
        });

        // ── Condition groups ─────────────────────────────────────────
        conditions.forEach((cond, ci) => {
            const gY = ci * groupH;

            // Condition label on Y-axis
            group.append('text')
                .attr('x', -4).attr('y', gY + groupH / 2)
                .attr('dy', '0.35em')
                .style('text-anchor', 'end')
                .style('font-size', C.chartLabelFontSize + 'px')
                .style('fill', '#555')
                .text(cond);

            // One bar per protein within this condition
            proteinList.forEach((p, pi) => {
                const avg = p.stats?.[cond]?.average ?? null;
                const std = p.stats?.[cond]?.std_dev ?? 0;
                if (avg === null) return;   // not detected – skip bar

                const barY = gY + pi * barH;
                const barW = xScale(avg);

                group.append('rect')
                    .attr('class', 'bar')
                    .attr('x', 0).attr('y', barY)
                    .attr('width',  Math.max(barW, 0))
                    .attr('height', barH - 1)
                    .attr('fill', colours[pi])
                    .on('mouseover', function () {
                        d3.select(this).attr('opacity', 0.7);
                        group.append('text').attr('class', 'value-label')
                            .attr('x', barW + 3).attr('y', barY + barH / 2)
                            .attr('dy', '0.35em')
                            .style('font-size', C.chartLabelFontSize + 'px')
                            .style('fill', '#111').style('font-weight', 'bold')
                            .text(`${d3.format('.2e')(avg)} ± ${d3.format('.1e')(std)}`);
                    })
                    .on('mouseout', function () {
                        d3.select(this).attr('opacity', 1);
                        group.selectAll('.value-label').remove();
                    });
            });

            // Thin separator between condition groups
            if (ci < conditions.length - 1) {
                group.append('line')
                    .attr('x1', 0).attr('x2', chartW - axisPad)
                    .attr('y1', gY + groupH).attr('y2', gY + groupH)
                    .attr('stroke', '#ddd').attr('stroke-width', 1);
            }
        });

        // X-axis at bottom
        group.append('g').attr('class', 'x-axis')
            .attr('transform', `translate(0,${chartH})`)
            .call(d3.axisBottom(xScale).ticks(3).tickFormat(d3.format('.2e')))
            .selectAll('text')
            .style('font-size', C.chartLabelFontSize + 'px').style('fill', '#555');

        if (C.barChartXLabel) {
            group.append('text')
                .attr('x', (chartW - axisPad) / 2).attr('y', chartH + 30)
                .style('text-anchor', 'middle')
                .style('font-size', C.chartLabelFontSize + 'px').style('fill', '#666')
                .text(C.barChartXLabel);
        }
    }

    // =========================================================================
    //  NODE BAR CHARTS  (metabolomics)
    // =========================================================================
    static createNodeBarCharts(jsonData) {
        const C         = EscherVisualizer.CONFIG;
        const isVert    = window.CONFIG?.smallGraphLayoutVertical ?? true;
        const threshold = window.CONFIG?.nodeThresholdSmall ?? 10;
        const numNodes  = d3.select('#map_container').selectAll('.metabolite-circle').size();
        const useLeft   = isVert && numNodes < threshold;

        d3.select('#map_container').selectAll('.node-circle.metabolite-circle')
            .each(function (data) {
                if (!data?.graph_info || Array.isArray(data.graph_info)) return;
                const element    = d3.select(this);
                const parentNode = d3.select(element.node().parentNode);
                EscherVisualizer._renderMetaboliteBarChart(
                    parentNode, element, data, {
                        getPosition: (d, cfg) => useLeft
                            ? {x: d.x + C.barChartOffsetY, y: d.y - cfg.chartWidth}
                            : {x: d.x - C.barChart.width / 2, y: d.y + C.barChartOffsetY},
                    }
                );
            });
    }

    // =========================================================================
    //  SEGMENT BAR CHARTS  (proteomics)
    // =========================================================================
    static createSegmentBarCharts(nodeData, reactions) {
        // Collect all reactant-edge segment data objects from the JSON
        // (D3 binds segment data to .segment elements, but we also
        //  need the graph_info list from the JSON directly)
        const segments = {};
        Object.values(reactions || {}).forEach(rxn => {
            Object.values(rxn.segments || {}).forEach(seg => {
                if (seg.edge_type === 'reactant_edge' &&
                    Array.isArray(seg.graph_info) &&
                    seg.graph_info.length > 0) {
                    segments[seg.from_node_id + '_' + seg.to_node_id] = seg;
                }
            });
        });

        d3.select('#map_container').selectAll('.segment')
            .each(function (data) {
                if (!data) return;
                // Match by from/to node IDs
                const key = (data.from_node_id || '') + '_' + (data.to_node_id || '');
                const seg = segments[key];
                if (!seg) return;

                const element    = d3.select(this);
                const parentNode = d3.select(element.node().parentNode);
                EscherVisualizer._renderProteomicsBarChart(
                    parentNode, element, seg, nodeData
                );
            });
    }

    // =========================================================================
    //  NODE STYLING
    // =========================================================================
    static stylePathwayElements() {
        d3.selectAll('.segment')
            .style('stroke', 'black')
            .style('stroke-width', '3px')
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
                d3.selectAll('circle.node-circle.metabolite-circle').each(function () {
                    if (parseFloat(this.getAttribute('r')) !== metR) needs = true;
                });
                if (needs) update();
            }).observe(container, {subtree: true, attributes: true, attributeFilter: ['r']});
            this._radiusObserverActive = true;
        }
    }

    static clearBarCharts() {
        d3.selectAll('.bar-chart').remove();
    }

    // =========================================================================
    //  HIGHLIGHTING
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