"use strict";

let alpha = 0.75;

function lsv_trace() {
    let trace = {
        x: [],
        y: [],
        type: 'scatter'
    };

    for (let x = 0; x <= 0.5; x+=1./128.) {
        trace.x.push(x);
        trace.y.push(LSV_left(x, alpha));
    }
    trace.x.push(0.5);
    trace.y.push(NaN);
    for (let x = 0.5; x <= 1.0; x+=1./4.) {
        trace.x.push(x);
        trace.y.push(LSV_right(x, alpha));
    }

    return trace;
}

function alphaChange() {
    alpha = alphaInput();
    document.getElementById("alphaValue").innerHTML = alpha.toFixed(4);
    Plotly.newPlot('lsv', [lsv_trace()]);

    let v = new Function('x', 'return ' + document.getElementById("v").value);

    Plotly.newPlot('Ln', computeLv(v,10), { yaxis: {range: [0.0, 2.0]}, xaxis: {dtick: 0.125}});

    let h = compute_h(14);
    Plotly.newPlot('hn', h[0], { yaxis: {range: [-8.0, 8.0]}, xaxis: {dtick: 0.125}});
    Plotly.newPlot('hph', h[1], { yaxis: {range: [0, 4.0]}, xaxis: {dtick: 0.125}});

    Plotly.newPlot('ccc', three_conditions(14), { yaxis: {range: [-1.0, 8.0]}, xaxis: {dtick: 0.125}});
}
function alphaInput() {
    let a = document.getElementById("alphaSlider").value / 64;
    document.getElementById("alphaValue").innerHTML = a.toFixed(4) + " (wait for calculations!)";
    return a;
}

function computeLv(v, n) {
    let traces = [];
    for (let k = 0; k <= n; k++) {
        let trace = {
            x: [],
            y: [],
            type: 'scatter',
            name: 'log(L^'+k.toString()+' v)',
            visible: 'legendonly'
        };
        for (let x = 0.0; x <= 1.0; x += 1./128) {
            trace.x.push(x);
            trace.y.push(Math.log(LSV_Ln(v, x, alpha, k))); // / LSV_Ln(v, x, alpha, 12));
        }
        traces.push(trace);
    }
    return traces;
}

function compute_h(N) {
    function vvv(x) {
        let gamma = 1. / alpha;
        let f = Math.pow(x, - gamma - 2) / 4;
        return [x * x * f, - gamma * f * x, gamma * (gamma + 1) * f];
    }
    let s = 'L^' + N.toString() + ' 1';
    let h = {
        x: [],
        y: [],
        name: s
    };
    let hp = {
        x: [],
        y: [],
        name: '(' + s + ")'"
    };
    let hpp = {
        x: [],
        y: [],
        name: '(' + s + ")''"
    };
    let hph = {
        x: [],
        y: [],
        name: '- x (' + s + ")'/" + s
    };
    let hpph = {
        x: [],
        y: [],
        name: 'x^2 (' + s + ")''/" + s
    };
    for (let x = 0.0; x <= 1.0; x += 1./128) {
        let g = LSV_Ln_pp(vvv, x, alpha, N);
        h.x.push(x);
        h.y.push(Math.log(g[0]));

        //hp.x.push(x);
        //hp.y.push(g[1]);

        //hpp.x.push(x);
        //hpp.y.push(g[2]);

        hph.x.push(x);
        hph.y.push(- x * g[1] / g[0]);

        hpph.x.push(x);
        hpph.y.push(x * x * g[2] / g[0]);
    }
    return [[h, hp, hpp], [hph, hpph]];
}

function three_conditions(N) {
    function vvv(x) {
        return [1, 0, 0];
    }
    function h(x) {
        return LSV_Ln_pp(vvv, x, alpha, N);
    }
    function COND(x) {
        let y1 = LSV_left_i(x, alpha),
            y2 = LSV_right_i(x, alpha),
            h1 = h(y1),
            h2 = h(y2),
            w1 = LSV_left_w(y1, alpha),
            w2 = LSV_right_w(y2, alpha),
            wp1 = LSV_left_wp(y1, alpha),
            wp2 = LSV_right_wp(y2, alpha),
            wpp1 = LSV_left_wpp(y1, alpha),
            wpp2 = LSV_right_wpp(y2, alpha);

        let A = h1[1] / h1[0] * w1 + wp1 - (h2[1] / h2[0] * w2 + wp2);
        A = -A;

        let B = h1[2] / h1[0] * w1 * w1 + 3 * h1[1] / h1[0] * wp1 * w1 + wpp1 * w1 + wp1 * wp1 -
            (   h2[2] / h2[0] * w2 * w2 + 3 * h2[1] / h2[0] * wp2 * w2 + wpp2 * w2 + wp2 * wp2 );

        let C = y1 / (1 + h2[0] * w2 / (h1[0] * w1))  +  y2 / (1 + h1[0] * w1 / (h2[0] * w2));
        C = x - C;

        return [A, B, C];
    }

    let c1 = {
        x: [],
        y: [],
        name: 'C1'
    };
    let c2 = {
        x: [],
        y: [],
        name: 'C2'
    };
    let c3 = {
        x: [],
        y: [],
        name: 'C3'
    };

    for (let x = 0.0; x <= 1.0; x += 1./128) {
        let c = COND(x);

        c1.x.push(x);
        c1.y.push(c[0]);

        c2.x.push(x);
        c2.y.push(c[1]);

        if (x >= 0.5) {
            c3.x.push(x);
            c3.y.push(c[2]);
        }
    }
    return [c1, c2, c3];
}

alphaChange();

let Lv_loop_n = 0;
function Lv_loop() {
    let N = document.getElementById('Ln').data.length;
    Plotly.restyle('Ln', {visible: 'legendonly', opacity: 1.0});
    if (Lv_loop_n < N - 1) {
        Plotly.restyle('Ln', {visible: 'true', opacity: 1.0}, Lv_loop_n);
    }
    Lv_loop_n = (Lv_loop_n + 1) % N;
    Plotly.restyle('Ln', {visible: true}, Lv_loop_n);
}
setInterval(Lv_loop, 1500);
