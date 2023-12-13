"use strict";

// normal random variable
function randn_bm(sigma2) {
    let u = 0, v = 0;
    while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
    while(v === 0) v = Math.random();
    return Math.sqrt( -2.0 * Math.log( u ) * sigma2 ) * Math.cos( 2.0 * Math.PI * v );
}

let graph_ctx = document.getElementById("graph").getContext('2d');
let myChart;
let graph_data = new TreeMap();
let ND = 1024;
let sigma2 = 1.0;

graph_data.set(0.0, 0.0);
//graph_data.set(1.0, randn_bm(sigma2));
graph_data.set(1.0, 1.0);

var downSamplePlugin = {
    beforeUpdate: function (chart, options) {
        let tmin = chart.scales["x-axis-0"].min;
        let tmax = chart.scales["x-axis-0"].max;

        if (!tmin) {
            tmin = 0.0;
        }
        if (!tmax) {
            tmax = 1.0;
        }

        let data = [];

        let H = 1;
        while (H * (tmax - tmin) < ND) {
            H *= 2;
        }

        let x0 = Math.floor(H * tmin) / H;

        for (let x = x0; x < tmax + 1/H; x += 1/H) {
            let it = graph_data.lowerBound(x);
            if (it.key == x) {
                data.push({x: x, y: graph_data.get(x)});
            } else {
                let s = {};
                if (it.equals(graph_data.end())) {
                    it.prev();
                    s = {x: x, y: it.value + randn_bm(x - it.key)};
                } else if (it.equals(graph_data.begin())) {
                    s = {x: x, y: it.value + randn_bm(it.key - x)};
                } else {
                    let right = {x: it.key, y: graph_data.get(it.key)};
                    it.prev();
                    let left = {x: it.key, y: graph_data.get(it.key)};
                    let mean = left.y + (x - left.x)*(right.y - left.y)/(right.x - left.x);
                    let varr = sigma2 * (right.x - x) * (x - left.x) / (right.x - left.x);
                    s = {x: x, y: mean + randn_bm(varr)};
                }
                data.push(s);
                graph_data.set(s.x, s.y);
            }
        }
        
        chart.data.datasets[0].data = data;
    }
};

function plot() {

    if (myChart == null) {
        myChart = new Chart(graph_ctx, {
            type: 'line',
            data: {
                datasets: [{
                    label: 'W(t)',
                    lineTension: 0,
                    steppedLine: 'before',
                    borderColor: 'rgba(255, 142, 0, 0.842)',
                    pointRadius: 0,
                    borderWidth: 2,
                    fill: false,
                    data: [{x:0.0, y:0.0}, {x:0.5,y:1.0}, {x:1.0, y:0.0}]
                }]
            },
            options: {
                scales: {
                    xAxes: [{
                        type: 'linear',
                        position: 'bottom',
                        ticks: {
                            maxRotation: 0
                        }
                    }]
                },
                pan: {
                    enabled: true,
                    mode: 'x'
                },
                zoom: {
                    enabled: true,
                    mode: 'x'
                }
            },
            plugins: [downSamplePlugin]
        });
    } else { // myChart is not null
        myChart.update();
        myChart.resetZoom();
    }
}

plot();
