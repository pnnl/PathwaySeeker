/**
 * EscherVisualizer - Handles Escher map visualization and D3 chart rendering.
 * Creates bar charts for node/segment data, styles pathway elements, and manages metabolite labels.
 */
class EscherVisualizer {
    static CONFIG = {
        imageSize: window.visualizationConfig?.imageSize ?? 200,
        barChart: {
            width: window.visualizationConfig?.barChartWidth ?? 100,
            height: window.visualizationConfig?.barChartHeight ?? 100,
            axisPadding: window.visualizationConfig?.barChartAxisPadding ?? 20,
            barHeight: window.visualizationConfig?.barHeight ?? 15  // Height per bar (px)
        },
        defaultToBiggId: false,  // false = show name first, true = show bigg_id first
        nodeRadius: window.visualizationConfig?.nodeRadius ?? 15,  // Node radius from backend config or fallback
        metaboliteRadius: window.visualizationConfig?.metaboliteRadius ?? 10,  // Metabolite node radius
        reactionRadius: window.visualizationConfig?.reactionRadius ?? 8,  // Reaction/coproduct node radius
        labelOffsetY: window.visualizationConfig?.labelOffsetY ?? 20,  // Label offset from backend config or fallback
        coproductLabelOffsetY: window.visualizationConfig?.coproductLabelOffsetY ?? 25,  // Coproduct label vertical offset
        barChartOffsetY: window.visualizationConfig?.barChartOffsetY ?? 60,  // Bar chart offset from backend config or fallback
        metaboliteLabelFontSize: window.visualizationConfig?.metaboliteLabelFontSize ?? 12,  // Metabolite label font size
        coproductLabelFontSize: window.visualizationConfig?.coproductLabelFontSize ?? 12,  // Coproduct label font size
        chartTitleFontSize: window.visualizationConfig?.chartTitleFontSize ?? 12,  // Chart title font size
        chartLabelFontSize: window.visualizationConfig?.chartLabelFontSize ?? 10,  // Chart axis label font size
        barChartTitle: window.visualizationConfig?.barChartTitle ?? "Data",  // Bar chart title text
        barChartXLabel: window.visualizationConfig?.barChartXLabel ?? "",  // Bar chart X-axis label
        barChartYLabel: window.visualizationConfig?.barChartYLabel ?? "Value"  // Bar chart Y-axis label
    };

    static initializeStructures(jsonData, hasImages = true) {
        this.clearBarCharts();
        this.stylePathwayElements(jsonData);
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

    static addMetaboliteLabels(jsonData) {
        d3.selectAll('g.node').each(function() {
            const circle = d3.select(this).select('.node-circle.metabolite-circle');
            if (circle.empty()) return;
            const data = circle.data()[0];
            if (!data) return;
            const parentNode = d3.select(this);
            
            // Hide the Escher-created label
            parentNode.selectAll(".label").style("display", "none");
            
            let label = parentNode.select(".node-label.metabolite-name");
            if (label.empty()) {
                label = parentNode.append("text")
                    .attr("class", "node-label metabolite-name")
                    .style("font-family", "Arial, sans-serif")
                    .style("font-size", EscherVisualizer.CONFIG.metaboliteLabelFontSize + "px")
                    .style("font-weight", "bold")
                    .style("fill", "#333")
                    .style("text-anchor", "middle")
                    .style("pointer-events", "none");
            }
            function updateLabel() {
                const transform = circle.attr("transform");
                if (transform) {
                    const match = transform.match(/translate\(([^,]+),([^)]+)\)/);
                    if (match) {
                        const x = parseFloat(match[1]);
                        const y = parseFloat(match[2]);
                        label.attr("transform", `translate(${x},${y + EscherVisualizer.CONFIG.labelOffsetY})`);
                    }
                }
                
                // Show name only, not bigg_id
                const labelText = (data.name || data.bigg_id || "Unknown").replace(/;\s*$/, "");
                
                // Always ensure font size is applied
                label.style("font-size", EscherVisualizer.CONFIG.metaboliteLabelFontSize + "px");
                label.text(labelText);
            }
            updateLabel();
            const observer = new MutationObserver(updateLabel);
            observer.observe(circle.node(), { attributes: true, attributeFilter: ["transform"] });
        });
    }

    static addCoproductLabels(jsonData) {
        d3.selectAll('.coproduct-circle').each(function() {
            const circle = d3.select(this);
            const data = circle.data()[0];
            if (!data) return;
            const parentNode = d3.select(this.parentNode);
            
            // Hide the Escher-created label
            parentNode.selectAll(".label").style("display", "none");
            
            let label = parentNode.select(".node-label.coproduct-name");
            if (label.empty()) {
                label = parentNode.append("text")
                    .attr("class", "node-label coproduct-name")
                    .style("font-family", "Arial, sans-serif")
                    .style("font-size", EscherVisualizer.CONFIG.coproductLabelFontSize + "px")
                    .style("fill", "#666")
                    .style("text-anchor", "middle")
                    .style("pointer-events", "none");
            }
            function updateLabel() {
                const transform = circle.attr("transform");
                if (transform) {
                    const match = transform.match(/translate\(([^,]+),([^)]+)\)/);
                    if (match) {
                        const x = parseFloat(match[1]);
                        const y = parseFloat(match[2]);
                        label.attr("transform", `translate(${x},${y - EscherVisualizer.CONFIG.coproductLabelOffsetY})`);
                    }
                }
                                
                // Use config boolean to determine which to show first
                const labelText = (EscherVisualizer.CONFIG.defaultToBiggId
                    ? (data.bigg_id || data.name || "Coproduct")
                    : (data.name || data.bigg_id || "Coproduct")).replace(/;\s*$/, "");
                
                // Always ensure font size is applied
                label.style("font-size", EscherVisualizer.CONFIG.coproductLabelFontSize + "px");
                label.text(labelText);
            }
            updateLabel();
            const observer = new MutationObserver(updateLabel);
            observer.observe(circle.node(), { attributes: true, attributeFilter: ["transform"] });
        });
    }

    static loadStructureImages(jsonData) {
        const config = EscherVisualizer.CONFIG;
        const maxDisplaySize = config.imageSize; // The largest molecule gets this size
        const isVerticalLayout = window.visualizationConfig?.smallGraphLayoutVertical ?? true;
        const nodeThreshold = window.visualizationConfig?.nodeThresholdSmall ?? 10;

        // Fetch image dimensions manifest, then render proportionally
        fetch('static/structure_imgs/image_dimensions.json')
            .then(resp => resp.ok ? resp.json() : {})
            .catch(() => ({}))
            .then(dimensions => {
                const circles = d3.select("#map_container").selectAll(".metabolite-circle");
                const numNodes = circles.size();

                // Use side positioning (image left, chart right) for small vertical graphs
                const useLeftPosition = isVerticalLayout && numNodes < nodeThreshold;

                // Use global max width across ALL known molecules so subgraphs
                // keep the same image scale as the full graph.
                let maxNaturalWidth = 0;
                for (const id of Object.keys(dimensions)) {
                    maxNaturalWidth = Math.max(maxNaturalWidth, dimensions[id].w);
                }
                if (maxNaturalWidth === 0) maxNaturalWidth = 1; // fallback

                const scaleFactor = maxDisplaySize / maxNaturalWidth;

                circles.each(function(data) {
                    const circle = d3.select(this);
                    const parentNode = d3.select(circle.node().parentNode);

                    if (data.highlight) circle.classed("highlighted", true);
                    circle.style("fill", "#999").style("stroke", "#666").style("fill-opacity", 0.7);

                    const imgPath = `static/structure_imgs/${data.bigg_id}.png`;
                    const dim = dimensions[data.bigg_id];

                    // Scale: proportional to natural size, largest = maxDisplaySize
                    const displayW = dim ? dim.w * scaleFactor : maxDisplaySize;
                    const displayH = dim ? dim.h * scaleFactor : maxDisplaySize;

                    // Position: left of node for small horizontal graphs, above node otherwise
                    let imgTranslateX, imgTranslateY;
                    if (useLeftPosition) {
                        // Image to the left of the node, bottom aligned with node
                        imgTranslateX = data.x - displayW - 10;
                        imgTranslateY = data.y - displayH;
                    } else {
                        // Image centered above the node (default)
                        imgTranslateX = data.x - displayW / 2;
                        imgTranslateY = data.y - displayH - 10;
                    }

                    const img = new Image();
                    img.onload = () => {
                        parentNode.insert("image", "text")
                            .attr("class", "structure-image")
                            .attr("transform", `translate(${imgTranslateX},${imgTranslateY})`)
                            .attr("width", displayW)
                            .attr("height", displayH)
                            .attr("xlink:href", imgPath);
                    };
                    img.onerror = () => { /* Silently skip missing structure images */ };
                    img.src = imgPath;

                    const updateImagePosition = () => {
                        const transform = circle.attr("transform");
                        if (transform) {
                            let offsetX, offsetY;
                            if (useLeftPosition) {
                                offsetX = -displayW - 10;
                                offsetY = -displayH;
                            } else {
                                offsetX = -displayW / 2;
                                offsetY = -displayH;
                            }
                            const updatedTransform = `${transform} translate(${offsetX},${offsetY})`;
                            parentNode.select("image").attr("transform", updatedTransform);
                        }
                    };
                    updateImagePosition();
                    const observer = new MutationObserver(updateImagePosition);
                    observer.observe(circle.node(), { attributes: true, attributeFilter: ["transform"] });
                });
            });
    }

static createBarCharts(selector, options = {}) {
    const defaults = {
        chartWidth: this.CONFIG.barChart.width,
        chartHeight: this.CONFIG.barChart.height,
        axisPadding: this.CONFIG.barChart.axisPadding,
        barColor: "#2a9d8f",
        hoverColor: "#1f7a67",
        showAxes: true,
        showBackground: true,
        showLabels: true,
        showValueOnHover: true,
        showTitle: false,
        fixedBarHeight: this.CONFIG.barChart.barHeight,  // ← Use configurable bar height
        positionOffset: { x: 0, y: 0 },
        getPosition: null,
        getParentNode: null,
        filterKeys: (key) => key !== "metabolite" && key !== "KEGG_C_number"
    };

    const config = Object.assign({}, defaults, options);

    const elements = d3.select("#map_container").selectAll(selector);
    
    elements.each(function(data, index) {
        if (!data || !data.graph_info) {
            return;
        }
        const element = d3.select(this);
        const parentNode = config.getParentNode ? config.getParentNode(element) : d3.select(element.node().parentNode);
        
        const graphInfo = data.graph_info;
        const keys = Object.keys(graphInfo).filter(config.filterKeys);
        const values = keys.map(key => {
            const average = graphInfo[key]?.average;
            return isNaN(average) ? 0 : average;
        });
        
        if (values.length === 0) {
            return;
        }

        const maxValue = Math.max(...values);
        const barHeight = config.fixedBarHeight;  // ← Use fixed height (configurable)
        const actualChartHeight = barHeight * values.length;  // ← Calculate total height
        
        const xScale = d3.scaleLinear().domain([0, maxValue]).range([0, config.chartWidth - (config.showAxes ? config.axisPadding : 0)]);

        // Show title if barChartTitle is configured (non-empty)
        const shouldShowTitle = config.showTitle || (EscherVisualizer.CONFIG.barChartTitle && EscherVisualizer.CONFIG.barChartTitle.length > 0);
        const titleHeight = shouldShowTitle && data.name ? 5 : 0;
        let position;
        
        if (config.getPosition) {
            position = config.getPosition(data, config);
            if (titleHeight > 0) position.y -= titleHeight / 2;
        } else {
            const x = (data.x || data.b1?.x || 0) - config.chartWidth + config.positionOffset.x;
            const y = (data.y || data.b1?.y || 0) - actualChartHeight / 2 + config.positionOffset.y - titleHeight;
            position = { x, y };
        }

        const chartGroup = parentNode.append("g")
            .attr("class", "bar-chart")
            .attr("transform", `translate(${position.x}, ${position.y})`);

        // Add observer to update bar chart position when node moves
        const updateChartPosition = () => {
            let newPosition;
            if (config.getPosition) {
                newPosition = config.getPosition(data, config);
                if (titleHeight > 0) newPosition.y -= titleHeight / 2;
            } else {
                const x = (data.x || data.b1?.x || 0) - config.chartWidth + config.positionOffset.x;
                const y = (data.y || data.b1?.y || 0) - actualChartHeight / 2 + config.positionOffset.y - titleHeight;
                newPosition = { x, y };
            }
            chartGroup.attr("transform", `translate(${newPosition.x}, ${newPosition.y})`);
        };

        // Observe the element for transform changes
        const observer = new MutationObserver(updateChartPosition);
        observer.observe(element.node(), { attributes: true, attributeFilter: ["transform"] });

        if (shouldShowTitle && data.name) {
            chartGroup.append("text")
                .attr("class", "chart-title")
                .attr("x", config.chartWidth / 2)
                .attr("y", -5)
                .attr("text-anchor", "middle")
                .style("font-size", EscherVisualizer.CONFIG.chartTitleFontSize + "px")
                .style("font-weight", "bold")
                .style("fill", "#333")
                .text((data.reaction_name || data.name).replace(/;\s*$/, ""));
        }

        if (config.showBackground) {
            chartGroup.insert("rect", ":first-child")
                .attr("width", config.chartWidth + 10)
                .attr("height", actualChartHeight + 20 + titleHeight)
                .attr("x", -10)
                .attr("y", titleHeight > 0 ? -10 - titleHeight : -10)
                .attr("fill", "white")
                .attr("opacity", 0.8)
                .attr("rx", 5);
        }

        if (config.showAxes) {
            // Y-AXIS: Center labels vertically within each bar
            const yScale = d3.scalePoint()
                .domain(keys)
                .range([titleHeight + barHeight / 2, actualChartHeight + titleHeight - barHeight / 2]);
            
            const yAxis = d3.axisLeft()
                .scale(yScale)
                .tickSize(0);

            chartGroup.append("g")
                .attr("class", "y-axis")
                .call(yAxis)
                .selectAll("text")
                .style("font-size", EscherVisualizer.CONFIG.chartLabelFontSize + "px")
                .style("fill", "#666")
                .style("text-anchor", "end")
                .style("dy", "0.35em")
                .attr("x", -5);

            // X-AXIS: Center ticks and labels on bar
            const xAxis = d3.axisBottom()
                .scale(xScale)
                .ticks(2)
                .tickFormat(d3.format(".2e"));

            chartGroup.append("g")
                .attr("class", "x-axis")
                .attr("transform", `translate(0, ${actualChartHeight + titleHeight})`)
                .call(xAxis)
                .selectAll("text")
                .style("font-size", EscherVisualizer.CONFIG.chartLabelFontSize + "px")
                .style("fill", "#666")
                .style("text-anchor", "middle");

            // Add Y-axis label if configured
            if (EscherVisualizer.CONFIG.barChartYLabel) {
                chartGroup.append("text")
                    .attr("class", "y-axis-label")
                    .attr("transform", "rotate(-90)")
                    .attr("y", 0 - config.chartWidth + 15)
                    .attr("x", 0 - (actualChartHeight + titleHeight) / 2)
                    .attr("dy", "1em")
                    .style("text-anchor", "middle")
                    .style("font-size", EscherVisualizer.CONFIG.chartLabelFontSize + "px")
                    .style("fill", "#666")
                    .text(EscherVisualizer.CONFIG.barChartYLabel);
            }

            // Add X-axis label if configured
            if (EscherVisualizer.CONFIG.barChartXLabel) {
                chartGroup.append("text")
                    .attr("class", "x-axis-label")
                    .attr("x", config.chartWidth / 2)
                    .attr("y", actualChartHeight + titleHeight + 35)
                    .attr("text-anchor", "middle")
                    .style("font-size", EscherVisualizer.CONFIG.chartLabelFontSize + "px")
                    .style("fill", "#666")
                    .text(EscherVisualizer.CONFIG.barChartXLabel);
            }
        }

        const barsSelection = chartGroup.selectAll(".bar")
            .data(values);

        barsSelection.enter()
            .append("rect")
            .attr("class", "bar")
            .attr("x", 0)
            .attr("y", (d, i) => i * barHeight + titleHeight)
            .attr("width", d => xScale(d))
            .attr("height", barHeight - 1)
            .attr("fill", config.barColor)
            .attr("stroke", "white")
            .attr("stroke-width", 0.5)
            .on("mouseover", function(d, i) {
                d3.select(this).attr("fill", config.hoverColor).attr("height", barHeight + 2);
                if (config.showValueOnHover) {
                    chartGroup.append("text")
                        .attr("class", "value-label")
                        .attr("x", xScale(d) + 3)
                        .attr("y", i * barHeight + barHeight / 2 + titleHeight)
                        .attr("dy", "0.35em")
                        .style("font-size", EscherVisualizer.CONFIG.chartLabelFontSize + "px")
                        .style("fill", "black")
                        .style("font-weight", "bold")
                        .text(d3.format(".2e")(d));
                }
            })
            .on("mouseout", function() {
                d3.select(this).attr("fill", config.barColor).attr("height", barHeight - 1);
                chartGroup.selectAll(".value-label").remove();
            });
    });
}
    static createNodeBarCharts(jsonData) {
        const isVerticalLayout = window.visualizationConfig?.smallGraphLayoutVertical ?? true;
        const nodeThreshold = window.visualizationConfig?.nodeThresholdSmall ?? 10;
        const numNodes = d3.select("#map_container").selectAll(".metabolite-circle").size();
        const useLeftPosition = isVerticalLayout && numNodes < nodeThreshold;

        this.createBarCharts(".node-circle.metabolite-circle", {
            chartWidth: this.CONFIG.barChart.width,
            chartHeight: this.CONFIG.barChart.height,
            barColor: "#2a9d8f",
            getPosition: (data, config) => {
                if (useLeftPosition) {
                    // Bar chart to the right of the node, bottom aligned with node
                    return { x: data.x + this.CONFIG.barChartOffsetY, y: data.y - config.chartHeight };
                }
                // Default: bar chart below the node, centered
                return { x: data.x - this.CONFIG.barChart.width / 2, y: data.y + this.CONFIG.barChartOffsetY };
            }
        });
        
    }

    static createSegmentBarCharts(nodeData) {
        this.createBarCharts(".segment", {
            chartWidth: this.CONFIG.barChart.width,
            chartHeight: this.CONFIG.barChart.height,
            barColor: "red",
            getPosition: (data) => {
                const fromNode = nodeData[data.to_node_id];
                return { x: (fromNode?.x || 0) - this.CONFIG.barChart.width, y: (fromNode?.y || 0) - this.CONFIG.barChart.height };
            },
            getParentNode: (element) => d3.select(element.node().parentNode)
        });
    }

    static stylePathwayElements(jsonData) {
        d3.selectAll(".segment").style("stroke", "black").style("stroke-width", "3px").style("stroke-dasharray", "5,5");
    }

    static equalizeNodeRadii() {
        // Use metabolite and reaction specific radii if available, otherwise fall back to nodeRadius
        const metaboliteRadius = this.CONFIG.metaboliteRadius ?? this.CONFIG.nodeRadius;
        const reactionRadius = this.CONFIG.reactionRadius ?? this.CONFIG.nodeRadius;
        
        const updateRadii = () => {
            // Update metabolite circles using attribute (SVG standard)
            d3.selectAll('circle.node-circle.metabolite-circle').each(function() {
                this.setAttribute('r', metaboliteRadius);
            });
            
            // Update coproduct circles using attribute
            d3.selectAll('circle.coproduct-circle').each(function() {
                this.setAttribute('r', reactionRadius);
            });
        };
        
        // Apply radius updates
        updateRadii();
        
        // Also set up a mutation observer to keep radii locked if Escher tries to change them
        const mapContainer = document.getElementById('map_container');
        if (mapContainer && !this.radiusObserverActive) {
            const radiusObserver = new MutationObserver(() => {
                // Check if radii changed, and restore them
                const metaboliteCircles = d3.selectAll('circle.node-circle.metabolite-circle');
                const coproductCircles = d3.selectAll('circle.coproduct-circle');
                
                let needsUpdate = false;
                metaboliteCircles.each(function() {
                    if (parseFloat(this.getAttribute('r')) !== metaboliteRadius) {
                        needsUpdate = true;
                    }
                });
                coproductCircles.each(function() {
                    if (parseFloat(this.getAttribute('r')) !== reactionRadius) {
                        needsUpdate = true;
                    }
                });
                
                if (needsUpdate) {
                    updateRadii();
                }
            });
            
            radiusObserver.observe(mapContainer, {
                subtree: true,
                attributes: true,
                attributeFilter: ['r']
            });
            
            this.radiusObserverActive = true;
            this.radiusObserver = radiusObserver;
        }
    }

    static initializeLabels(jsonData) {
        this.removeStoichiometryLabels();
        this.addMetaboliteLabels(jsonData);
        this.addCoproductLabels(jsonData);
    }

    static clearBarCharts() {
        d3.selectAll(".bar-chart").remove();
    }

    static highlightPath(pathNodes) {
        d3.selectAll('.path-node').classed('highlighted-path', false);
        d3.selectAll('.path-segment').classed('highlighted-path', false);

        d3.selectAll('circle.node-circle.metabolite-circle').classed('path-node', function() {
            const data = d3.select(this).data()[0];
            return data && pathNodes.includes(data.bigg_id);
        }).classed('highlighted-path', function() {
            return d3.select(this).classed('path-node');
        });

        d3.selectAll('path.segment').classed('path-segment', function() {
            const data = d3.select(this).data()[0];
            if (!data) return false;
            const fromNode = data.to_node_id;
            return fromNode && pathNodes.includes(fromNode);
        }).classed('highlighted-path', function() {
            return d3.select(this).classed('path-segment');
        });
    }

    static highlightMultiNodes(selectedNodeIds) {
        d3.selectAll('.multi-node').classed('highlighted-multi', false);
        d3.selectAll('.multi-segment').classed('highlighted-multi', false);

        d3.selectAll('circle.node-circle.metabolite-circle').classed('multi-node', function() {
            const data = d3.select(this).data()[0];
            return data && selectedNodeIds.includes(data.bigg_id);
        }).classed('highlighted-multi', function() {
            return d3.select(this).classed('multi-node');
        });

        d3.selectAll('path.segment').classed('multi-segment', function() {
            const data = d3.select(this).data()[0];
            if (!data) return false;
            const fromNode = data.to_node_id;
            return fromNode && selectedNodeIds.includes(fromNode);
        }).classed('highlighted-multi', function() {
            return d3.select(this).classed('multi-segment');
        });
    }
}
