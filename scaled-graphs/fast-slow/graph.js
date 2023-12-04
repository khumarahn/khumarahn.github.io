"use strict";

function set_defaults() {
    let sel = document.getElementById("preset").value;

    let f = document.getElementById("f");
    let x = document.getElementById("x0");
    let o = document.getElementById("observable");
    let c = document.getElementById("cn");
    let s0 = document.getElementById("s0");
    let n = document.getElementById("n");
    let pvar_p =  document.getElementById("pvar_p");
    if (sel == "SLSV:Levy") {
        f.value = "SLSV(x, 0.6)";
        x.value = "0.2";
        s0.value = "{x: 1.0, y: 0.0}";
        o.value = "{x: s.x - eps * 4 * s.y * (x-0.5), y: s.y + eps * 4 * s.x * (x-0.5)}";
        c.value = "Math.pow(n, -0.6)";
    } else if (sel == "CUSP:CKM23") {
        f.value = "CUSP(x, 1/3, 3)";
        x.value = "{x: 0.5, y: 0.0, theta: 0.23, phi: 0}";
        s0.value = "{x: 0.0, y: 0.0}";
        o.value = "{x: s.x + eps * Math.cos(3*x.phi), y: s.y + eps * Math.cos(5*x.phi)}";
        c.value = "Math.pow(n, -2/3)";
        n.value = "5000";
        pvar_p.value = "3/2 + 1/4";
    }
}

let input = {}; //f,phi,x0,n,cn,pvar_p,pvar_compute;

let variation = {}; //pvar,onevar,qvar;

let graph_ctx = document.getElementById("graph").getContext('2d');
let myChart;
let graph_data = null;

var downSamplePlugin = {
    beforeUpdate: function (chart, options) {

        if (graph_data == null) { // no init yet
            return;
        }

        function straighten(d, dx, dy) {
            if (d.length > 2 && d[d.length - 1] != null && d[d.length - 2] != null && d[d.length - 3] != null) {
                let x1 = d[d.length - 3].x;
                let y1 = d[d.length - 3].y;
                let x2 = d[d.length - 2].x;
                let y2 = d[d.length - 2].y;
                let x3 = d[d.length - 1].x;
                let y3 = d[d.length - 1].y;
                let ddx = 0.0, ddy = 0.0;
                let a = (x2 - x1) * (x3 - x1) + (y2 - y1) * (y3 - y1);
                let b = (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1);
                if (y3 != y1 && x3 != x1 && a >= 0 && a <= b) {
                    ddx = Math.abs(x1 - x2 + (x3 - x1) * (y2 - y1) / (y3 - y1));
                    ddy = Math.abs(y1 - y2 + (y3 - y1) * (x2 - x1) / (x3 - x1));
                    if (ddx < dx && ddy < dy) {
                        let dl = d.pop();
                        d.pop();
                        d.push(dl);
                    }
                }
            }
        }
        
        let xmin = chart.scales["x-axis-0"].min;
        let xmax = chart.scales["x-axis-0"].max;
        if (!xmin) {
            xmin = 0.0;
        }
        if (!xmax) {
            xmax = 0.0;
        }

        let n = input.n;
        
        let ymin = chart.scales["y-axis-0"].min;
        let ymax = chart.scales["y-axis-0"].max;
        if (!ymin) {
            ymin = 0.0;
        }
        if (!ymax) {
            ymax = 0.0;
        }

        let dx = (xmax - xmin) / 512;
        let dy = (ymax - ymin) / 512;
        
        let data_x = [null];
        let data_y = [null];
        let data_xy = [null];
        
        for (let k = 0; k < graph_data.jxy.length; k++) {
            let j = graph_data.jxy[k].j;
            let x = graph_data.jxy[k].x;
            let y = graph_data.jxy[k].y;
            let t = j / n;

            let data_x_last = data_x[data_x.length - 1];
            let data_y_last = data_y[data_y.length - 1];
            let data_xy_last = data_xy[data_xy.length - 1];
            if (t >= xmin && t <= xmax) {
                if (data_x_last == null || Math.abs(t - data_x_last.x) > dx || Math.abs(x - data_x_last.y) > dy) {
                    data_x.push({x:t, y:x});
                }
                if (data_y_last == null || Math.abs(t - data_y_last.x) > dx || Math.abs(y - data_y_last.y) > dy) {
                    data_y.push({x:t, y:y});
                }
            } else {
                if (data_x_last != null) {
                    data_x.push(null);
                }
                if (data_y_last != null) {
                    data_y.push(null);
                }
            }

            if ((x >= xmin && x <= xmax && y >= ymin && y <= ymax)
                || (data_xy_last != null && data_xy_last.x >= xmin && data_xy_last.x <= xmax && data_xy_last.y >= ymin && data_xy_last.y <= ymax)) {
                if (data_xy_last == null || Math.abs(x - data_xy_last.x) > dx || Math.abs(y - data_xy_last.y) > dy) {
                    data_xy.push({t:t, x:x, y:y});
                }
            } else {
                if (data_xy_last != null) {
                    data_xy.push(null);
                }
            }

            //straighten(data_x, dx, dy);
            //straighten(data_y, dx, dy);
            //straighten(data_xy, dx, dy);
        }

        //chart.data.datasets[0].pointRadius = 0;
        chart.data.datasets[0].data = data_x;
        chart.data.datasets[1].data = data_y;
        chart.data.datasets[2].data = data_xy;
    }
};

function update_graph_data() {
    if (graph_data == null) {
        graph_data = {};
    }
    let f = input.f;
    let n = input.n;
    let phi = input.phi;
    let x0 = input.x0;
    let s0 = input.s0;
    let pvar_p = input.pvar_p;

    // normalization
    let c = input.cn(n);
    graph_data.c = c;
    
    let euc = (a, b) => Math.sqrt(a*a + b*b);

    let s = s0;
    let x = x0;

    let qvar = 0.0;
    let onevar = 0.0;

    graph_data.jxy = [];

    for (let j=0; j<n; j++) {
        graph_data.jxy.push({j: j, x: s.x, y: s.y});
        if (graph_data.jxy.length > 1) {
            let p1 = graph_data.jxy[graph_data.jxy.length - 1];
            let p2 = graph_data.jxy[graph_data.jxy.length - 2];
            let d12 = euc(p1.x - p2.x, p1.y - p2.y);
            qvar += d12*d12;
            onevar += d12;
        }

        s = phi(s, x, c);
        x = f(x);
    }

    variation.onevar = onevar;
    variation.qvar = Math.sqrt(qvar);

    let pvar_dist = (a, b) => euc(graph_data.jxy[b].x - graph_data.jxy[a].x, graph_data.jxy[b].y - graph_data.jxy[a].y);
    let pvar = p_var_backbone(graph_data.jxy.length, pvar_p, pvar_dist);
    variation.pvar = Math.pow(pvar.p_var, 1./pvar_p);
    variation.pvar_points = [];
    variation.pvar_all = pvar;
    for (let q=0; q<pvar.points.length; q++) {
        let qq = graph_data.jxy[pvar.points[q]];
        variation.pvar_points.push({x: qq.x, y: qq.y});
    }
}

function init_plot() {

    if (myChart == null) {
        myChart = new Chart(graph_ctx, {
            type: 'line',
            data: {
                datasets: [
                    {
                        label: 'X(t)',
                        lineTension: 0,
                        steppedLine: 'before',
                        borderColor: 'rgba(142, 0, 0, 0.842)',
                        pointRadius: 0,
                        borderWidth: 2,
                        fill: false,
                        data: [{x:0.0, y:0.0}, {x:0.5,y:1.0}, {x:1.0, y:0.0}],
                        hidden: true
                    },
                    {
                        label: 'Y(t)',
                        lineTension: 0,
                        steppedLine: 'before',
                        borderColor: 'rgba(0, 42, 255, 0.842)',
                        pointRadius: 0,
                        borderWidth: 2,
                        fill: false,
                        data: [{x:0.0, y:0.0}, {x:0.5,y:1.0}, {x:1.0, y:0.0}],
                        hidden: true
                    },
                    {
                        label: 'XY',
                        lineTension: 0,
                        borderColor: 'rgba(0, 142, 0, 0.842)',
                        pointRadius: 0,
                        borderWidth: 2,
                        fill: false,
                        data: [{x:0.0, y:0.0}, {x:0.5,y:0.5}, {x:0.0, y:1.0}]
                    },
                    {
                        label: 'p-variation points',
                        backgroundColor: 'rgba(0, 0, 0)',
                        pointRadius: 0,
                        borderWidth: 0,
                        fill: false,
                        showLine: false,
                        pointRadius: 3,
                        pointHoverRadius: 5,
                        data: [{x:0.0, y:0.0}, {x:0.5,y:0.5}, {x:0.0, y:1.0}]
                    }
                ]
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
                    mode: 'xy'
                },
                zoom: {
                    enabled: true,
                    mode: 'xy'
                }
            },
            plugins: [downSamplePlugin]
        });
    } else { // myChart is not null
        myChart.resetZoom();
    }

}

function LSV_(x, alpha) {
  return x * (1 + Math.pow(2*x, alpha));
}

function LSV(x, alpha) {
  if (x <= 0) {
    return 0;
  }
  if (x >= 1) {
    return 1;
  }
  if (x <= 0.5) {
    return LSV_(x, alpha);
  } else {
    return 2*x - 1;
  }
}

function SLSV(x, alpha) {
  if (x <= 0) {
    return 0;
  }
  if (x >= 1) {
    return 1;
  }
function LSV_(x, alpha) {
  return x * (1 + Math.pow(2*x, alpha));
}
  if (x <= 0.5) {
    return LSV_(x, alpha);
  } else {
    return 1 - LSV_(1 - x, alpha);
  }
}

function LSV3_(x, alpha) {
  return x * (1 + Math.pow(3*x, alpha));
}

function SLSV3(x, alpha) {
  if (x <= 0) {
    return 0;
  }
  if (x >= 1) {
    return 1;
  }
  if (x <= 1./3.) {
    return LSV3_(x, alpha);
  } else if (x <= 2./3.) {
    return 3.*x - 1.;
  } else {
    return 1 - LSV3_(1 - x, alpha);
  }
}

function downloadText(filename,text){
    // Set up the link
    let link = document.createElement("a");
    link.setAttribute("target","_blank");
    if(Blob !== undefined) {
        let blob = new Blob([text], {type: "text/plain"});
        link.setAttribute("href", URL.createObjectURL(blob));
    } else {
        link.setAttribute("href","data:text/plain," + encodeURIComponent(text));
    }
    link.setAttribute("download",filename);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

function get_tikz_data() {
    let s = "t\tx\ty\n";
    let data_xy = myChart.data.datasets[2].data;
    for (let k=0; k< data_xy.length; k++) {
        if (data_xy[k] == null || data_xy[k].x == null || data_xy[k].y == null) {
            s += "\n";
        } else {
            s += data_xy[k].t.toExponential(6) + "\t" +  data_xy[k].x.toExponential(6) + "\t" + data_xy[k].y.toExponential(6) + "\n";
        }
    }
    downloadText("data.txt", s);
}

function form_click() {
    let i_f =  document.getElementById("f").value;
    let i_n = document.getElementById("n").value;
    let i_cn =  document.getElementById("cn").value;
    let i_x0 = document.getElementById("x0").value;
    let i_s0 = document.getElementById("s0").value;
    let i_phi =  document.getElementById("observable").value;
    let i_pvar_p =  document.getElementById("pvar_p").value;

    input.f = new Function('x', 'return ' + i_f);
    input.phi = new Function('s', 'x', 'eps', 'return ' + i_phi);
    input.cn = new Function('n', 'return ' + i_cn);
    input.x0 = eval("(" + i_x0 + ")");
    input.s0 = eval("(" + i_s0 + ")");
    input.n = eval(i_n);
    input.pvar_p = eval(i_pvar_p);

    update_graph_data();

    myChart.data.datasets[0].data = graph_data.x;
    myChart.data.datasets[1].data = graph_data.y;
    myChart.data.datasets[2].data = graph_data.xy;
    myChart.data.datasets[3].data = variation.pvar_points;
    myChart.update();
    myChart.resetZoom();
    myChart.update();
    let info = document.getElementById("info");
    info.innerHTML =
        input.pvar_p.toString() + '-variation: ' + variation.pvar.toFixed(4) +  ';&emsp;&emsp;&emsp;&emsp;' +
        '1-variation: ' + variation.onevar.toFixed(4) + ';&emsp;&emsp;&emsp;&emsp;' +
        'quadratic variation: ' + variation.qvar.toFixed(4);

}

init_plot();
form_click();
