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

    document.getElementById('thm').innerHTML += '<br>... formulating for alpha ≈ ' + alpha.toFixed(6) + ' ...';

    worker.postMessage({
        gamma: 1. / alpha
    });

    MathJax.typeset();
}

function alphaInput() {
    let a = 0.5 + document.getElementById("alphaSlider").value / 256;
    document.getElementById("alphaValue").innerHTML = a.toFixed(4) + " (wait for calculations!)";
    return a;
}

let worker = new Worker('lsv-worker.js');

worker.onmessage = function(event) {
    let message = event.data;

    let thm = String.raw`
        For \( \alpha \approx ::alpha:: \), the invariant density \(h(x)\)
        for \(x \in (0, 1]\) satisfies:
        \[
            \Bigl( \frac{h'}{h} \Bigr)' \geq ::hp::
            \quad \text{and} \quad
            \Bigl( \frac{h''}{h} \Bigr)' \leq ::hpp::
            .
        \]
    `;

    thm = thm.replaceAll('::alpha::', message.alpha.toFixed(6))
        .replaceAll('::hp::', message.min_hp_h_prime.toFixed(6))
        .replaceAll('::hpp::', message.max_hpp_h_prime.toFixed(6));

    document.getElementById("thm").innerHTML = thm;

    MathJax.typeset();

    console.log("gg gamma:", message.gamma);

};

function startup() {
    if (!(window.MathJax && window.Plotly)) {
        setTimeout(startup, 1000);
        return;
    }

    alphaChange();

    worker.postMessage({
        gamma: 1. / alpha
    });
}

startup();
