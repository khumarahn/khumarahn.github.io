"use strict";

let alpha = 0.75;

function SQR(x) {
    return x * x;
}

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

    Plotly.newPlot('Ln', computeLv(v,10), { yaxis: {range: [0, 1.25]}, xaxis: {range: [0,6], rangemode: 'tozero'}});

    Plotly.newPlot('hn', compute_abel_h(), { yaxis: {autorange: true, rangemode: 'tozero'}, xaxis: {autorange: true, rangemode: 'tozero'}});
    Plotly.newPlot('ccc', four_conditions(), { yaxis: {type: 'log', autorange: true}, xaxis: {range: [0,1], dtick: 0.125}});
    Plotly.newPlot('qqq', q1q2(), { yaxis: {autorange: true}, xaxis: {range: [0,1], dtick: 0.125}});

    let R_cols = lsv_cpp.R_cols();
    document.getElementById("R_cols").innerHTML = R_cols.toString() + "x" + R_cols.toString();
    let R_div = document.getElementById("R_matrix");
    let R_html = "\\[ R = \\begin{pmatrix} ";
    for (let i = 0; i < 5; i++) {
        for (let j = 0; j < 5; j++) {
            R_html += lsv_cpp.R_coef(i,j).toFixed(6) + " & ";
        }
        if (i == 0) {
            R_html += "\\ldots";
        }
        R_html += " & \\text{" + lsv_cpp.R_coef(i, R_cols - 1).toExponential(3) + "}";
        R_html += "\\\\";
    }
    R_html += "\\vdots & & & & & \\ddots & \\\\";
    for (let j = 0; j < 5; j++) {
        R_html += "\\text{" + lsv_cpp.R_coef(R_cols - 1, j).toExponential(3) + "}" + " & ";
    }
    R_html += " & \\text{" + lsv_cpp.R_coef(R_cols - 1, R_cols - 1).toExponential(3) + "}";
    R_html += "\\end{pmatrix} \\]";
    R_div.innerHTML = R_html;

    let R_evalues_div = document.getElementById("R_evalues");
    R_evalues_div.innerHTML = "";
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

    document.getElementById("lambda_error").innerHTML = "\\( \\varepsilon = \\text{"
        + (lsv_cpp.R_evalues_real(0) - 1).toExponential(2)
        + "} + \\text{" + lsv_cpp.R_evalues_imag(0).toExponential(2) + "} i \\)";
    document.getElementById("h12_error").innerHTML = "\\( \\text{ "
        + (lsv_cpp.h(0.5) / lsv_cpp.h(1) - (2 + lsv_cpp.gamma()) / 2).toExponential(2) + "} \\)";

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
        for (let u = 1./16; u <= 12.0; u += 1./128) {
            let x = lsv_cpp.abel_inv(u),
                p = - lsv_cpp.abel_inv_p(u);
            trace.x.push(u);
            trace.y.push(LSV_Ln(v, x, alpha, k) * p / abel_h(u));
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

function abel_h(u) {
    let x = lsv_cpp.abel_inv(u),
        p = - lsv_cpp.abel_inv_p(u);

    return lsv_cpp.h(x) * p;
}

function compute_abel_h() {
    let h_A = {
        x: [],
        y: [],
        name: 'h_A'
    };
    for (let u = 1./16; u <= 12.0; u += 1./32) {
        h_A.x.push(u);
        h_A.y.push(abel_h(u));
    }
    return [h_A];

}

function four_conditions() {
    function h(x) {
        return [lsv_cpp.h(x), lsv_cpp.h_p(x), lsv_cpp.h_pp(x)];
    }
    function COND(x) {
        let y1 = LSV_left_i(x, alpha),
            y2 = LSV_right_i(x, alpha),
            //
            hx = lsv_cpp.h(x),
            hpx = lsv_cpp.h_p(x),
            hppx = lsv_cpp.h_pp(x),
            //
            h1 = lsv_cpp.h(y1),
            hp1 = lsv_cpp.h_p(y1),
            hpp1 = lsv_cpp.h_pp(y1),
            //
            h2 = lsv_cpp.h(y2),
            hp2 = lsv_cpp.h_p(y2),
            hpp2 = lsv_cpp.h_pp(y2),
            //
            w1 = LSV_left_w(y1, alpha),
            wp1 = LSV_left_wp(y1, alpha),
            wpp1 = LSV_left_wpp(y1, alpha),
            //
            w2 = LSV_right_w(y2, alpha),
            wp2 = LSV_right_wp(y2, alpha),
            wpp2 = LSV_right_wpp(y2, alpha),
            //
            q1 = w1 * h1 / hx,
            q2 = w2 * h2 / hx,
            //
            qp1 = (
                (hp1 * w1 + h1 * wp1) / hx
                - h1 * hpx / (hx * hx)
            ),
            qp2 = (
                (hp2 * w2 + h2 * wp2) / hx
                - h2 * hpx / (hx * hx)
            ),
            qpp1 = (
                (hpp1 * w1 + 2 * hp1 * wp1 + h1 * wpp1) / hx 
                - ( 2 * hp1 * w1 * hpx + h1 * wp1 * hpx + h1 * hppx) / (hx * hx * w1)
                + 2 * h1 * hpx * hpx / (hx * hx * hx * w1)
            ),
            qpp2 = (
                (hpp2 * w2 + 2 * hp2 * wp2 + h2 * wpp2) / hx 
                - ( 2 * hp2 * w2 * hpx + h2 * wp2 * hpx + h2 * hppx) / (hx * hx * w2)
                + 2 * h2 * hpx * hpx / (hx * hx * hx * w2)
            );

        let A = y1 / (1 + h2 * w2 / (h1 * w1))  +  y2 / (1 + h1 * w1 / (h2 * w2));
        A = x - A;

        let B = qp2;

        let C = - qpp2;

        let D = - (
            + (2 * qp1 * w1 + q1 * wp1) * w1
            + (2 * qp2 * w2 + q2 * wp2) * w2
        );

        let E1 = qp1 * w1 + qp2 * w2;

        let E2 = (
            + (qpp1 * w1 + qp1 * wp1) * w1
            + (qpp2 * w2 + qp2 * wp2) * w2
        );

        return [A, B, C, D];
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
    let c4 = {
        x: [],
        y: [],
        name: 'C4'
    };

    for (let x = 1./16; x <= 1.0; x += 1./128) {
        let c = COND(x);

        if (x >= 0.5) {
            c1.x.push(x);
            c1.y.push(c[0]);
        }

        c2.x.push(x);
        c2.y.push(c[1]);

        c3.x.push(x);
        c3.y.push(c[2]);

        c4.x.push(x);
        c4.y.push(c[3]);
    }
    return [c1, c2, c3, c4];
}

function q1q2() {
    function QQ(x) {
        let y1 = LSV_left_i(x, alpha),
            y2 = LSV_right_i(x, alpha),
            //
            hx = lsv_cpp.h(x),
            //
            h1 = lsv_cpp.h(y1),
            //
            h2 = lsv_cpp.h(y2),
            //
            w1 = LSV_left_w(y1, alpha),
            //
            w2 = LSV_right_w(y2, alpha),
            //
            q1 = w1 * h1 / hx,
            q2 = w2 * h2 / hx;

        return [y1, q1, y2, q2];
    }

    let c1 = {
        x: [],
        y: [],
        name: 'q(x) (left)'
    };
    let c2 = {
        x: [],
        y: [],
        name: 'q(x) (right)'
    };

    for (let x = 1./16; x <= 1.0; x += 1./128) {
        let c = QQ(x);

        c1.x.push(c[0]);
        c1.y.push(c[1]);

        c2.x.push(c[2]);
        c2.y.push(c[3]);
    }
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

var lsv_cpp;

var Module = {
    'onRuntimeInitialized': function() {
        lsv_cpp = new Module.LSV();

        alphaChange();
        setInterval(Lv_loop, 1500);
    }
};
