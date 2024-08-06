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

    lsv_cpp.set_gamma(1. / alpha);

    let v = new Function('x', 'return ' + document.getElementById("v").value);

    Plotly.newPlot('Ln', computeLv(v,10), { yaxis: {range: [0.0, 2.0]}, xaxis: {range: [0,1], dtick: 0.125}});

    let h = compute_h();
    Plotly.newPlot('hn', h[0], { yaxis: {range: [-24.0, 24.0]}, xaxis: {range: [0,1], dtick: 0.125}});
    Plotly.newPlot('hph', h[1], { yaxis: {autoscale: true}, xaxis: {range: [0,1], dtick: 0.125}});

    Plotly.newPlot('ccc', three_conditions(14), { yaxis: {type: 'log', autorange: true}, xaxis: {range: [0,1], dtick: 0.125}});


    let R_size = lsv_cpp.R_size();
    document.getElementById("R_size").innerHTML = R_size.toString() + "x" + R_size.toString();
    let R_div = document.getElementById("R_matrix");
    let R_html = "\\[ R = \\begin{pmatrix} ";
    for (let i = 0; i < 5; i++) {
        for (let j = 0; j < 5; j++) {
            R_html += lsv_cpp.R_coef(i,j).toFixed(6) + " & ";
        }
        if (i == 0) {
            R_html += "\\ldots";
        }
        R_html += " & \\text{" + lsv_cpp.R_coef(i, R_size - 1).toExponential(3) + "}";
        R_html += "\\\\";
    }
    R_html += "\\vdots & & & & & \\ddots & \\\\";
    for (let j = 0; j < 5; j++) {
        R_html += "\\text{" + lsv_cpp.R_coef(R_size - 1, j).toExponential(3) + "}" + " & ";
    }
    R_html += " & \\text{" + lsv_cpp.R_coef(R_size - 1, R_size - 1).toExponential(3) + "}";
    R_html += "\\end{pmatrix} \\]";
    R_div.innerHTML = R_html;

    let R_evalues_div = document.getElementById("R_evalues");
    for (let j = 0; j < 5; j++) {
        R_evalues_div.innerHTML += " \\(" + lsv_cpp.R_evalues_real(j).toFixed(4) + "+" + lsv_cpp.R_evalues_imag(j).toFixed(4) + "i \\)";
            R_evalues_div.innerHTML += ",";
    }
    R_evalues_div.innerHTML += "...";

    let R_evector_div = document.getElementById("R_evector");
    R_evector_div.innerHTML = "\\[ \\begin{pmatrix} ";
    for (let j = 0; j < 5; j++) {
        R_evector_div.innerHTML += lsv_cpp.h_coef(j).toFixed(6) + " \\\\ ";
    }
    R_evector_div.innerHTML += "\\vdots \\end{pmatrix} \\]";

    MathJax.typeset();
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
    let hph = {
        x: [],
        y: [],
        name: "- x h'(x) / h(x)"
    };
    let hpph = {
        x: [],
        y: [],
        name: "x^2 h''(x) / h(x)"
    };
    for (let x = 1./16; x <= 1.0; x += 1./128) {
        h.x.push(x);
        h.y.push(lsv_cpp.h(x));

        hp.x.push(x);
        hp.y.push(lsv_cpp.h_p(x));

        hpp.x.push(x);
        hpp.y.push(lsv_cpp.h_pp(x));

        hph.x.push(x);
        hph.y.push(-x * lsv_cpp.h_p(x) / lsv_cpp.h(x));

        hpph.x.push(x);
        hpph.y.push(x * x * lsv_cpp.h_pp(x) / lsv_cpp.h(x));
    }
    return [[h, hp, hpp], [hph, hpph]];
}

function three_conditions(N) {
    function h(x) {
        return [lsv_cpp.h(x), lsv_cpp.h_p(x), lsv_cpp.h_pp(x)];
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
    return [c1, c2, c3];
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

var lsv_cpp;

var Module = {
    'onRuntimeInitialized': function() {
        lsv_cpp = new Module.LSV();

        alphaChange();
        setInterval(Lv_loop, 1500);
    }
};
