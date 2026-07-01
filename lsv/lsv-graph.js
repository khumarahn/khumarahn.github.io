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
    Plotly.newPlot('lsv', [lsv_trace()],{ yaxis: {range: [0,1.02], dtick: 0.25}, xaxis: {range: [0,1.02], dtick: 0.25}});

    // lsv_cpp.set_gamma(1. / alpha);

    let v = new Function('x', 'return ' + document.getElementById("v").value);

    Plotly.newPlot('Ln', computeLv(v,10), { yaxis: {range: [0.0, 2.0]}, xaxis: {range: [0,1], dtick: 0.125}});

    let h = compute_h();
    Plotly.newPlot('hn', h[0], { yaxis: {range: [-24.0, 24.0]}, xaxis: {range: [0,1], dtick: 0.125}});
    Plotly.newPlot('hph', h[1], { yaxis: {autoscale: true}, xaxis: {range: [0,1], dtick: 0.125}});

    Plotly.newPlot('ccc', three_conditions(14), { yaxis: {type: 'log', autorange: true}, xaxis: {range: [0,1], dtick: 0.125}});

    if (typeof MathJax !== 'undefined') {
        MathJax.typeset();
    }
}

function alphaInput() {
    let a = document.getElementById("alphaSlider").value / 64;
    document.getElementById("alphaValue").innerHTML = a.toFixed(4) + " (wait for calculations!)";
    return a;
}

function computeLv(v, n) {
    let traces = [];
    for (let k = 0; k <= n; k++) {if (typeof variable !== 'undefined') {
    // the variable is defined
}
        let trace = {
            x: [],
            y: [],
            type: 'scatter',
            name: 'L^'+k.toString()+' v',
            visible: 'legendonly'
        };
        for (let x = 0.0; x <= 1.0; x += 1./128) {
            trace.x.push(x);
            trace.y.push(LSV_Ln(v, x, alpha, k));
        }
        traces.push(trace);
    }
    return traces;
}

function compute_h() {
    let s = 'h';
    let h = {
        x: [],
        y: [],
        name: s
    };
    let hp = {
        x: [],
        y: [],
        name: "h'"
    };
    let hpp = {
        x: [],
        y: [],
        name: "h''"
    };
    for (let x = 1./16; x <= 1.0; x += 1./128) {
        h.x.push(x);
        h.y.push(x); //lsv_cpp.h(x));

        hp.x.push(x);
        hp.y.push(1.2 * x); //lsv_cpp.h_p(x));

        hpp.x.push(x);
        hpp.y.push(1.3 * x); //lsv_cpp.h_pp(x));
    }
    return [[h, hp, hpp]];
}

function three_conditions(N) {
    function h(x) {
        //return [lsv_cpp.h(x), lsv_cpp.h_p(x), lsv_cpp.h_pp(x)];
        return [ 1 / x, - 1 / (x * x), 2 / (x * x * x)];
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

    for (let x = 1./16; x <= 1.0; x += 1./128) {
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
    //return [c1, c2, c3];
    return [c1, c2];
}

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

let worker = new Worker("lsv-worker.js");

worker.onmessage = function(event) {
    let message = event.data;

    console.log("gg gamma:", message.gamma);

    alphaChange();
};

function startup() {
    if (!(window.MathJax && window.Plotly)) {
        setTimeout(startup, 1000);
        return;
    }

    alphaChange();
    setInterval(Lv_loop, 1500);

    let worker = new Worker("lsv-worker.js");

    worker.postMessage({
        gamma: 1.2
    });
}

startup();
