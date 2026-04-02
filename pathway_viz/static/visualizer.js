/**
 * EscherVisualizer - D3 rendering overlays on top of Escher maps.
 * Reads all configuration from window.CONFIG (set by Flask template).
 */
class EscherVisualizer {

    static get CONFIG() {
        const C = window.CONFIG || {};
        return {
            imageSize:        C.imageSize  ?? 200,
            defaultToBiggId:  false,
            nodeRadius:        C.nodeRadius        ?? 15,
            metaboliteRadius:  C.metaboliteRadius  ?? 10,
            reactionRadius:    C.reactionRadius    ?? 8,
            labelOffsetY:      C.labelOffsetY      ?? 20,
            coproductLabelOffsetY: C.coproductLabelOffsetY ?? 25,
            barChartOffsetY:   C.barChartOffsetY   ?? 60,
            metaboliteLabelFontSize: C.metaboliteLabelFontSize ?? 12,
            coproductLabelFontSize:  C.coproductLabelFontSize  ?? 12,
            chartTitleFontSize: C.chartTitleFontSize ?? 12,
            chartLabelFontSize: C.chartLabelFontSize ?? 10,
            barChartTitle:  C.barChartTitle  ?? '',
            barChartXLabel: C.barChartXLabel ?? '',
            barChartYLabel: C.barChartYLabel ?? '',
            barChart: {
                width:       C.barChartWidth       ?? 100,
                height:      C.barChartHeight      ?? 100,
                axisPadding: C.barChartAxisPadding ?? 20,
                barHeight:   C.barHeight           ?? 15,
            },
        };
    }

    static initializeStructures(jsonData) {
        this.clearBarCharts();
        this.stylePathwayElements();
        this.equalizeNodeRadii();
        this.loadStructureImages(jsonData);
        this.createNodeBarCharts(jsonData);
        this.createSegmentBarCharts(jsonData[1].nodes);
        this.initializeLabels(jsonData);
    }

    static setupNaNLabelRemoval() {
        setInterval(() => this.removeStoichiometryLabels(), 500);
    }

    static removeStoichiometryLabels() {
        d3.selectAll('.stoichiometry-label').filter(function() {
            return d3.select(this).text() === 'NaN';
        }).remove();
    }

    // ─────────────────────────────────────────────────────────────────────────
    // LABELS
    // ─────────────────────────────────────────────────────────────────────────

    static addMetaboliteLabels(jsonData) {
        d3.selectAll('g.node').each(function() {
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
                .observe(circle.node(), { attributes: true, attributeFilter: ['transform'] });
        });
    }

    static addCoproductLabels(jsonData) {
        d3.selectAll('.coproduct-circle').each(function() {
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
                .observe(circle.node(), { attributes: true, attributeFilter: ['transform'] });
        });
    }

    static initializeLabels(jsonData) {
        this.removeStoichiometryLabels();
        this.addMetaboliteLabels(jsonData);
        this.addCoproductLabels(jsonData);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // STRUCTURE IMAGES
    // ─────────────────────────────────────────────────────────────────────────

    static loadStructureImages(jsonData) {
        const C             = EscherVisualizer.CONFIG;
        const maxDisplaySize = C.imageSize;
        const isVertical    = window.CONFIG?.smallGraphLayoutVertical ?? true;
        const nodeThreshold = window.CONFIG?.nodeThresholdSmall ?? 10;

        fetch('static/structure_imgs/image_dimensions.json')
            .then(r => r.ok ? r.json() : {})
            .catch(() => ({}))
            .then(dimensions => {
                const circles  = d3.select('#map_container').selectAll('.metabolite-circle');
                const numNodes = circles.size();
                const useLeft  = isVertical && numNodes < nodeThreshold;

                let maxNatW = Object.values(dimensions).reduce((m, d) => Math.max(m, d.w), 0);
                if (maxNatW === 0) maxNatW = 1;
                const scale = maxDisplaySize / maxNatW;

                circles.each(function(data) {
                    const circle     = d3.select(this);
                    const parentNode = d3.select(circle.node().parentNode);

                    if (data.highlight) circle.classed('highlighted', true);
                    circle.style('fill', '#999').style('stroke', '#666').style('fill-opacity', 0.7);

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
                            .attr('class', 'structure-image')
                            .attr('transform', `translate(${data.x + o.x},${data.y + o.y})`)
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
                    new MutationObserver(updatePos)
                        .observe(circle.node(), { attributes: true, attributeFilter: ['transform'] });
                });
            });
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BAR CHARTS
    // ─────────────────────────────────────────────────────────────────────────

    static createBarCharts(selector, options = {}) {
        const C = EscherVisualizer.CONFIG;
        const defaults = {
            chartWidth:       C.barChart.width,
            chartHeight:      C.barChart.height,
            axisPadding:      C.barChart.axisPadding,
            fixedBarHeight:   C.barChart.barHeight,
            barColor:         '#2a9d8f',
            hoverColor:       '#1f7a67',
            showAxes:         true,
            showBackground:   true,
            showValueOnHover: true,
            positionOffset:   { x: 0, y: 0 },
            getPosition:      null,
            getParentNode:    null,
            filterKeys: key => key !== 'metabolite' && key !== 'KEGG_C_number',
        };
        const cfg = Object.assign({}, defaults, options);

        d3.select('#map_container').selectAll(selector).each(function(data) {
            if (!data?.graph_info) return;
            const C2 = EscherVisualizer.CONFIG;

            const element    = d3.select(this);
            const parentNode = cfg.getParentNode
                ? cfg.getParentNode(element)
                : d3.select(element.node().parentNode);

            const keys   = Object.keys(data.graph_info).filter(cfg.filterKeys);
            const values = keys.map(k => {
                const avg = data.graph_info[k]?.average;
                return isNaN(avg) ? 0 : avg;
            });
            if (!values.length) return;

            const maxVal    = Math.max(...values);
            const barHeight = cfg.fixedBarHeight;
            const chartH    = barHeight * values.length;
            const showTitle = C2.barChartTitle && C2.barChartTitle.length > 0;
            const titleH    = showTitle && data.name ? 5 : 0;

            const xScale = d3.scaleLinear()
                .domain([0, maxVal])
                .range([0, cfg.chartWidth - (cfg.showAxes ? cfg.axisPadding : 0)]);

            let pos = cfg.getPosition
                ? cfg.getPosition(data, cfg)
                : {
                    x: (data.x || data.b1?.x || 0) - cfg.chartWidth + cfg.positionOffset.x,
                    y: (data.y || data.b1?.y || 0) - chartH / 2 + cfg.positionOffset.y,
                };
            if (titleH) pos = { x: pos.x, y: pos.y - titleH / 2 };

            const group = parentNode.append('g')
                .attr('class', 'bar-chart')
                .attr('transform', `translate(${pos.x},${pos.y})`);

            const updatePos = () => {
                let np = cfg.getPosition ? cfg.getPosition(data, cfg) : {
                    x: (data.x || data.b1?.x || 0) - cfg.chartWidth + cfg.positionOffset.x,
                    y: (data.y || data.b1?.y || 0) - chartH / 2 + cfg.positionOffset.y,
                };
                if (titleH) np = { x: np.x, y: np.y - titleH / 2 };
                group.attr('transform', `translate(${np.x},${np.y})`);
            };
            new MutationObserver(updatePos)
                .observe(element.node(), { attributes: true, attributeFilter: ['transform'] });

            if (showTitle && data.name) {
                group.append('text')
                    .attr('class', 'chart-title')
                    .attr('x', cfg.chartWidth / 2).attr('y', -5)
                    .attr('text-anchor', 'middle')
                    .style('font-size', C2.chartTitleFontSize + 'px')
                    .style('font-weight', 'bold').style('fill', '#333')
                    .text((data.reaction_name || data.name).replace(/;\s*$/, ''));
            }

            if (cfg.showBackground) {
                group.insert('rect', ':first-child')
                    .attr('width',  cfg.chartWidth + 10)
                    .attr('height', chartH + 20 + titleH)
                    .attr('x', -10)
                    .attr('y', titleH ? -10 - titleH : -10)
                    .attr('fill', 'white').attr('opacity', 0.8).attr('rx', 5);
            }

            if (cfg.showAxes) {
                const yScale = d3.scalePoint()
                    .domain(keys)
                    .range([titleH + barHeight / 2, chartH + titleH - barHeight / 2]);

                group.append('g').attr('class', 'y-axis')
                    .call(d3.axisLeft().scale(yScale).tickSize(0))
                    .selectAll('text')
                    .style('font-size', C2.chartLabelFontSize + 'px')
                    .style('fill', '#666').style('text-anchor', 'end')
                    .style('dy', '0.35em').attr('x', -5);

                group.append('g').attr('class', 'x-axis')
                    .attr('transform', `translate(0,${chartH + titleH})`)
                    .call(d3.axisBottom().scale(xScale).ticks(2).tickFormat(d3.format('.2e')))
                    .selectAll('text')
                    .style('font-size', C2.chartLabelFontSize + 'px')
                    .style('fill', '#666').style('text-anchor', 'middle');

                if (C2.barChartYLabel) {
                    group.append('text').attr('class', 'y-axis-label')
                        .attr('transform', 'rotate(-90)')
                        .attr('y', 0 - cfg.chartWidth + 15)
                        .attr('x', 0 - (chartH + titleH) / 2)
                        .attr('dy', '1em').style('text-anchor', 'middle')
                        .style('font-size', C2.chartLabelFontSize + 'px').style('fill', '#666')
                        .text(C2.barChartYLabel);
                }
                if (C2.barChartXLabel) {
                    group.append('text').attr('class', 'x-axis-label')
                        .attr('x', cfg.chartWidth / 2)
                        .attr('y', chartH + titleH + 35)
                        .attr('text-anchor', 'middle')
                        .style('font-size', C2.chartLabelFontSize + 'px').style('fill', '#666')
                        .text(C2.barChartXLabel);
                }
            }

            group.selectAll('.bar').data(values).enter().append('rect')
                .attr('class', 'bar')
                .attr('x', 0)
                .attr('y', (d, i) => i * barHeight + titleH)
                .attr('width',  d => xScale(d))
                .attr('height', barHeight - 1)
                .attr('fill',   cfg.barColor)
                .attr('stroke', 'white').attr('stroke-width', 0.5)
                .on('mouseover', function(d, i) {
                    const C3 = EscherVisualizer.CONFIG;
                    d3.select(this).attr('fill', cfg.hoverColor).attr('height', barHeight + 2);
                    if (cfg.showValueOnHover) {
                        group.append('text').attr('class', 'value-label')
                            .attr('x', xScale(d) + 3)
                            .attr('y', i * barHeight + barHeight / 2 + titleH)
                            .attr('dy', '0.35em')
                            .style('font-size', C3.chartLabelFontSize + 'px')
                            .style('fill', 'black').style('font-weight', 'bold')
                            .text(d3.format('.2e')(d));
                    }
                })
                .on('mouseout', function() {
                    d3.select(this).attr('fill', cfg.barColor).attr('height', barHeight - 1);
                    group.selectAll('.value-label').remove();
                });
        });
    }

    static createNodeBarCharts(jsonData) {
        const C         = EscherVisualizer.CONFIG;
        const isVert    = window.CONFIG?.smallGraphLayoutVertical ?? true;
        const threshold = window.CONFIG?.nodeThresholdSmall ?? 10;
        const numNodes  = d3.select('#map_container').selectAll('.metabolite-circle').size();
        const useLeft   = isVert && numNodes < threshold;

        this.createBarCharts('.node-circle.metabolite-circle', {
            getPosition: (data, cfg) => useLeft
                ? { x: data.x + C.barChartOffsetY, y: data.y - cfg.chartHeight }
                : { x: data.x - C.barChart.width / 2, y: data.y + C.barChartOffsetY },
        });
    }

    static createSegmentBarCharts(nodeData) {
        const C = EscherVisualizer.CONFIG;
        this.createBarCharts('.segment', {
            barColor: 'red',
            getPosition: data => {
                const from = nodeData[data.to_node_id];
                return {
                    x: (from?.x || 0) - C.barChart.width,
                    y: (from?.y || 0) - C.barChart.height,
                };
            },
            getParentNode: el => d3.select(el.node().parentNode),
        });
    }

    // ─────────────────────────────────────────────────────────────────────────
    // NODE STYLING & RADII
    // ─────────────────────────────────────────────────────────────────────────

    static stylePathwayElements() {
        d3.selectAll('.segment')
            .style('stroke', 'black')
            .style('stroke-width', '3px')
            .style('stroke-dasharray', '5,5');
    }

    static equalizeNodeRadii() {
        const C    = EscherVisualizer.CONFIG;
        const metR  = C.metaboliteRadius ?? C.nodeRadius;
        const reacR = C.reactionRadius   ?? C.nodeRadius;

        const update = () => {
            d3.selectAll('circle.node-circle.metabolite-circle')
                .each(function() { this.setAttribute('r', metR); });
            d3.selectAll('circle.coproduct-circle')
                .each(function() { this.setAttribute('r', reacR); });
        };
        update();

        const container = document.getElementById('map_container');
        if (container && !this._radiusObserverActive) {
            const obs = new MutationObserver(() => {
                let needs = false;
                d3.selectAll('circle.node-circle.metabolite-circle').each(function() {
                    if (parseFloat(this.getAttribute('r')) !== metR) needs = true;
                });
                d3.selectAll('circle.coproduct-circle').each(function() {
                    if (parseFloat(this.getAttribute('r')) !== reacR) needs = true;
                });
                if (needs) update();
            });
            obs.observe(container, {
                subtree: true, attributes: true, attributeFilter: ['r']
            });
            this._radiusObserverActive = true;
            this._radiusObserver       = obs;
        }
    }

    static clearBarCharts() {
        d3.selectAll('.bar-chart').remove();
    }

    // ─────────────────────────────────────────────────────────────────────────
    // HIGHLIGHTING
    // ─────────────────────────────────────────────────────────────────────────

    static highlightPath(pathNodes) {
        d3.selectAll('.path-node').classed('highlighted-path', false);
        d3.selectAll('.path-segment').classed('highlighted-path', false);

        d3.selectAll('circle.node-circle.metabolite-circle')
            .classed('path-node', function() {
                const d = d3.select(this).data()[0];
                return d && pathNodes.includes(d.bigg_id);
            })
            .classed('highlighted-path', function() {
                return d3.select(this).classed('path-node');
            });

        d3.selectAll('path.segment')
            .classed('path-segment', function() {
                const d = d3.select(this).data()[0];
                return d && pathNodes.includes(d.to_node_id);
            })
            .classed('highlighted-path', function() {
                return d3.select(this).classed('path-segment');
            });
    }

    static highlightMultiNodes(selectedNodeIds) {
        d3.selectAll('.multi-node').classed('highlighted-multi', false);
        d3.selectAll('.multi-segment').classed('highlighted-multi', false);

        d3.selectAll('circle.node-circle.metabolite-circle')
            .classed('multi-node', function() {
                const d = d3.select(this).data()[0];
                return d && selectedNodeIds.includes(d.bigg_id);
            })
            .classed('highlighted-multi', function() {
                return d3.select(this).classed('multi-node');
            });

        d3.selectAll('path.segment')
            .classed('multi-segment', function() {
                const d = d3.select(this).data()[0];
                return d && selectedNodeIds.includes(d.to_node_id);
            })
            .classed('highlighted-multi', function() {
                return d3.select(this).classed('multi-segment');
            });
    }
}