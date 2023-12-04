"use strict";

function set_defaults() {
    let sel = document.getElementById("preset").value;

    let f = document.getElementById("f");
    let x = document.getElementById("x0");
    let o = document.getElementById("observable");
    let c = document.getElementById("cn");
    let a = document.getElementById("observable_bridge");
    if (sel == "Gauss:Brownian") {
        f.value = "(1/x) % 1";
        x.value = "0.18";
        o.value = "x + 1 - 1 / Math.log(2)";
        c.value = "Math.pow(n, -0.5)";
        a.checked = false;
    }
    if (sel == "LSV:Levy") {
        f.value = "LSV(x, 0.6)";
        x.value = "0.2";
        o.value = "x";
        c.value = "Math.pow(n, -0.6)";
        a.checked = true;
    }
    if (sel == "LSV:Brownian") {
        f.value = "LSV(x, 0.25)";
        x.value = "0.2";
        o.value = "x";
        c.value = "Math.pow(n, -0.5)";
        a.checked = true;
    }
    if (sel == "SLSV:Levy") {
        f.value = "SLSV(x, 0.6)";
        x.value = "0.2";
        o.value = "x-0.5";
        c.value = "Math.pow(n, -0.6)";
        a.checked = false;
    }
    if (sel == "SLSV:Brownian") {
        f.value = "SLSV(x, 0.25)";
        x.value = "0.2";
        o.value = "x-0.5";
        c.value = "Math.pow(n, -0.5)";
        a.checked = false;
    }
    if (sel == "SLSV:uBrownian") {
        f.value = "SLSV(x, 0.6)";
        x.value = "0.2";
        o.value = "Math.abs(x-0.5) < 0.4 ? (x-0.5) : 0";
        c.value = "Math.pow(n, -0.5)";
        a.checked = false;
    }
    if (sel == "CUSP:Levy") {
        f.value = "CUSP(x, 1/3, 3)";
        x.value = "({x: 0.5, y: 0.0, theta: 0.23, phi: 0})";
        o.value = "Math.cos(3*x.phi)";
        c.value = "Math.pow(n, -2/3)";
        a.checked = false;
    }
    if (sel == "CUSP:LevyM") {
        f.value = "CUSP(x, 1/3, 3)";
        x.value = "({x: 0.5, y: 0.0, theta: 0.23, phi: 0})";
        o.value = "((Math.abs(x.phi) > Math.PI/2.1) ? -1.0 : 0.0028038)";
        c.value = "Math.pow(n, -2/3)";
        a.checked = false;
    }
    if (sel == "CUSP:Brownian") {
        f.value = "CUSP(x, 1/3, 3)";
        x.value = "({x: 0.5, y: 0.0, theta: 0.23, phi: 0})";
        o.value = "x.x < 0.25 ? 0 : Math.cos(3*x.phi)";
        c.value = "Math.pow(n, -0.5)";
        a.checked = false;
    }
    if (sel == "3x:Brownian") {
        f.value = "(3*x) % 1";
        x.value = "0.2";
        o.value = "Math.cos(2*Math.PI*x)";
        c.value = "Math.pow(n, -0.5)";
        a.checked = false;
    }
    if (sel == "3x:Levy") {
        f.value = "(3*x) % 1";
        x.value = "0.2";
        o.value = "Math.pow(x, -0.6) - 1/(1-0.6)";
        c.value = "Math.pow(n, -0.6)";
        a.checked = false;
    }
}

function getJsonFromUrl(url) {
    if(!url) url = location.search;
    var query = url.substr(1);
    var result = {};
    query.split("&").forEach(function(part) {
        var item = part.split("=");
        result[item[0]] = decodeURIComponent(item[1]);
    });
    return result;
}

document.getElementById("preset").value = getJsonFromUrl().preset;
set_defaults();


let input = {}; //f,phi,phi_bridge,x0,n,cn,pvar_p,pvar_compute;

let variation = {}; //pvar,onevar,qvar;

let graph_ctx = document.getElementById("graph").getContext('2d');
let myChart;
let graph_data = {};

var downSamplePlugin = {
    beforeUpdate: function (chart, options) {
        let tmin = chart.scales["x-axis-0"].min;
        let tmax = chart.scales["x-axis-0"].max;

        let n = input.n;
        let f = input.f;
        let phi = input.phi;
        let c = graph_data.c;
        let a = graph_data.a;

        if (!tmin) {
            tmin = 0.0;
        } else if (tmin < 0) {
            tmin = 0;
        }
        if (!tmax) {
            tmax = 1.0;
        } else if (tmax > 1.0) {
            tmax = 1.0;
        }

        let data = [];

        let kmin = 0;
        let kmax = graph_data.data.length - 1;
        while (kmax - kmin > 1) {
            let k = Math.floor((kmin + kmax)/2);
            if (graph_data.data[k].j < n*tmin) {
                kmin = k;
            } else {
                kmax = k;
            }
        }
        let k = kmin;

        let s = graph_data.data[k].s;
        let x = graph_data.data[k].x;
        let j = graph_data.data[k].j;
        let dt = (tmax - tmin) / 1024;

        data.push({x:j/n,y:s});

        {
            chart.data.datasets[0].pointRadius = 0;
            let local_min = null, 
                local_max = null;
            while (j <= n*tmax) {
                let t = j/n;
                if (local_min == null || s < local_min.s) {
                    local_min = {t:t, s:s};
                }
                if (local_max == null || s > local_max.s) {
                    local_max = {t:t, s:s};
                }
                if (t < tmin) {
                    data[0] = {x:t, y:s};
                    local_min = null;
                    local_max = null;
                } else if (t - data[data.length - 1].x < dt) {
                } else {
                    if (local_min.t == local_max.t) {
                        data.push({x:local_min.t, y:local_min.s});
                    } else if (local_min.t < local_max.t) {
                        data.push({x:local_min.t, y:local_min.s});
                        data.push({x:local_max.t, y:local_max.s});
                    } else {
                        data.push({x:local_max.t, y:local_max.s});
                        data.push({x:local_min.t, y:local_min.s});
                    }
                    local_min = null;
                    local_max = null;
                }

                if (k + 1 < graph_data.data.length && graph_data.data[k + 1].j - j < n*dt) {
                    k++;
                    s = graph_data.data[k].s;
                    x = graph_data.data[k].x;
                    j = graph_data.data[k].j;
                } else {
                    s += c*(phi(x) - a);
                    x = f(x);
                    j++;
                }
            }
        }
        chart.data.datasets[0].data = data;
        //chart.options.scales.yAxes[0].ticks.suggestedMin = - graph_data.max_abs_s;
        //chart.options.scales.yAxes[0].ticks.suggestedMax = + graph_data.max_abs_s;
    }
};

function update_graph_data() {
    let f = input.f;
    let n = input.n;
    let phi = input.phi;
    let x0 = input.x0;

    // average?
    let a = 0;
    if (input.phi_bridge) {
        let x = x0;
        for (let j=0; j<2*n; j++) {
            x = f(x);
        }
        let s = 0;
        for (let j=0; j<2*n; j++) {
            s += phi(x);
            x = f(x);
        }
        a = s / (2*n);
    }
    graph_data.a = a;

    // normalization
    let c = input.cn(n);
    graph_data.c = c;

    let s = 0;
    let max_abs_s = 0;
    let x = x0;

    let qvar = 0.0;
    let onevar = 0.0;
    let points_pv = [];

    graph_data.data = [];
    let graph_data_every = Math.floor(Math.sqrt(n));

    let local_min, local_max;
    for (let j=0; j<n; j++) {
        
        if (j == 0 || s < local_min.s) {
            local_min = {j:j, x:x, s:s};
        }
        if (j == 0 || s > local_max.s) {
            local_max = {j:j, x:x, s:s};
        }
        
        if (j % graph_data_every == 1) {
            local_min = {j:j, x:x, s:s};
            local_max = {j:j, x:x, s:s};
        }

        if (j % graph_data_every == 0) {
            if (local_min.j == local_max.j) {
                graph_data.data.push({j: local_min.j, x: local_min.x, s: local_min.s});
            } else if (local_min.j < local_max.j) {
                graph_data.data.push({j: local_min.j, x: local_min.x, s: local_min.s});
                graph_data.data.push({j: local_max.j, x: local_max.x, s: local_max.s});
            } else {
                graph_data.data.push({j: local_max.j, x: local_max.x, s: local_max.s});
                graph_data.data.push({j: local_min.j, x: local_min.x, s: local_min.s});
            }
        }

        if (Math.abs(s) > max_abs_s) {
            max_abs_s = Math.abs(s);
        }

        let inc = c*(phi(x) - a);
        s += inc;
        x = f(x);

        qvar += inc*inc;
        onevar += Math.abs(inc);
        if (pvar_compute) {
            // add s to points_pv; only keep endpoints and local extrema
            let p_size = points_pv.length;
            if (p_size < 2) {
                points_pv.push({j: j, s: s});
            } else {
                let p_last1 = points_pv[p_size - 1].s;
                let p_last2 = points_pv[p_size - 2].s;
                if ( (p_last2 <= p_last1 && p_last1 <= s) ||
                    (p_last2 >= p_last1 && p_last1 >= s)) {
                    points_pv[p_size - 1] = {j: j, s: s};
                } else {
                    points_pv.push({j: j, s: s});
                }
            }
        }
    }

    graph_data.max_abs_s = max_abs_s;

    variation.onevar = onevar;
    variation.qvar = Math.sqrt(qvar);
    if (input.pvar_compute) {
        let dist = (a, b) => Math.abs(points_pv[b].s - points_pv[a].s);
        let p_var = p_var_backbone(points_pv.length, input.pvar_p, dist);
        variation.pvar = Math.pow(p_var.p_var, 1./input.pvar_p);
        variation.pvar_points = [];
        for (let q = 0; q < p_var.points.length; q++) {
            let qq = points_pv[p_var.points[q]];
            variation.pvar_points.push({x: qq.j / n, y: qq.s});
        }
    }
}

function plot() {

    if (myChart == null) {
        myChart = new Chart(graph_ctx, {
            type: 'line',
            data: {
                datasets: [{
                    label: 'S(t)',
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
        myChart.update();
        myChart.resetZoom();
    }
    if (input.pvar_compute) {
        let pvar_data = {
            label: 'Points for p-variation',
            backgroundColor: 'rgba(0, 0, 0, 1.0)',
            borderWidth: 0,
            fill: false,
            showLine: false,
            pointRadius: 3,
            pointHoverRadius: 5,
            data: variation.pvar_points
        };
        if (myChart.data.datasets.length < 2) {
            myChart.data.datasets.push(pvar_data);
        } else {
            myChart.data.datasets[1] = pvar_data;
        }
        myChart.update();
    } else {
        if (myChart.data.datasets.length > 1) {
            myChart.data.datasets.pop();
            myChart.update();
        }
    }

    let info = document.getElementById("info");
    info.innerHTML = '';
    if (input.pvar_compute) {
        info.innerHTML = input.pvar_p.toFixed(2) + '-variation: ' + variation.pvar.toFixed(4) + ';&emsp;&emsp;&emsp;&emsp;';
    }
    info.innerHTML += 
        '1-variation: ' + variation.onevar.toFixed(4) + ';&emsp;&emsp;&emsp;&emsp;' +
        'quadratic variation: ' + variation.qvar.toFixed(4);

}

function download_graph_data() {

    let data_text = "x\ty\n";

    if (myChart.data.datasets.length > 0) {
        for (let i = 0; i < myChart.data.datasets[0].data.length; i++) {
            let x = myChart.data.datasets[0].data[i].x;
            let y = myChart.data.datasets[0].data[i].y;
            data_text += x.toExponential(6) + "\t" + y.toExponential(6) + "\n";
        }
    }

    // https://stackoverflow.com/questions/3665115
    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(data_text));
    element.setAttribute('download', 'data.txt');

    element.style.display = 'none';
    document.body.appendChild(element);

    element.click();

    document.body.removeChild(element);
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
  if (x <= 0.5) {
    return LSV_(x, alpha);
  } else {
    return 1 - LSV_(1 - x, alpha);
  }
}

function form_click() {
  let i_f =  document.getElementById("f").value;
  let i_n = document.getElementById("n").value;
  let i_cn =  document.getElementById("cn").value;
  let i_x0 = document.getElementById("x0").value;
  let i_phi =  document.getElementById("observable").value;
  let i_pvar_p =  document.getElementById("pvar_p").value;

  input.phi_bridge = document.getElementById("observable_bridge").checked;

  input.f = new Function('x', 'return ' + i_f);
  input.phi = new Function('x', 'return ' + i_phi);
  input.cn = new Function('n', 'return ' + i_cn);
  input.x0 = eval(i_x0);
  input.n = eval(i_n);
  
  input.pvar_p = eval(i_pvar_p);
  input.pvar_compute = document.getElementById("pvar_compute").checked;

  update_graph_data();
  plot();
}

form_click();
