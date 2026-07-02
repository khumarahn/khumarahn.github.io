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
    Plotly.newPlot('lsv', [lsv_trace()],{ yaxis: {range: [0,1.02], dtick: 0.25}, xaxis: {range: [0,1.02], dtick: 0.25}});

    if (!workerBusy) {
        document.getElementById('theorem-control').innerHTML = `
            <br>
                <button id="theorem-btn" style=" background-color: #007bff; color: white; border: none; 
                    border-radius: 0.3rem; padding: 0.5rem 1rem; font-size: 1rem; cursor: pointer;  margin-top: 0.5rem;">
                    Formulate the result for \\(\\alpha=${alpha}\\)
                </button>
            `;

        document.getElementById('theorem-btn').addEventListener('click', function(e) {
            e.preventDefault(); // no page jump

            document.getElementById('theorem-control').innerHTML =
                `<span>Computing bounds for \\(\\alpha \\approx ${alpha}\\)...</span>`;

            document.getElementById('reasoning-log').textContent += '\n\n';
            document.getElementById('reasoning-summary').textContent = 'Reasoning...';

            worker.postMessage({ type: 'compute-bounds', gamma: 1 / alpha });
            workerBusy = true;

            MathJax.typeset();
        });
    }

    MathJax.typeset();
}

function alphaInput() {
    let a = 0.5 + document.getElementById("alphaSlider").value / 256;
    document.getElementById("alphaValue").innerHTML = `${a}`;
    return a;
}

let worker = new Worker('lsv-worker.js');
let workerBusy = false;

worker.onmessage = function(e) {
    let m = e.data;

    if (m.type === 'stdout') {
        let log = document.getElementById('reasoning-log');

        log.textContent += e.data.text + '\n';
        log.scrollTop = log.scrollHeight; // auto-scroll to bottom
    } else if (m.type === 'bounds') {

        let thm = String.raw`
        <div>
        For
        \begin{align*}
            \alpha \in [
            & ::alpha_minus:: , \\
            & ::alpha_plus:: ],
        \end{align*}
        \(h(x)\) on \((0, 1]\) satisfies:
        \[
            \Bigl( \frac{h'}{h} \Bigr)' \geq ::hp::
            \quad \text{and} \quad
            \Bigl( \frac{h''}{h} \Bigr)' \leq ::hpp::
            .
        \]
        </div>
       `;

        console.log("huj", m.alpha_minus);
        thm = thm.replaceAll('::alpha_minus::', m.alpha_minus)
            .replaceAll('::alpha_plus::', m.alpha_plus)
            .replaceAll('::hp::', m.min_hp_h_prime)
            .replaceAll('::hpp::', m.max_hpp_h_prime);

        document.getElementById("theorem-text").innerHTML += thm;

        document.getElementById('theorem-control').innerHTML = '';
        document.getElementById('reasoning-summary').innerHTML = 'Done reasoning';

        workerBusy = false;

        MathJax.typeset();
    }
};

function startup() {
    if (!(window.MathJax && window.Plotly)) {
        setTimeout(startup, 1000);
        return;
    }

    // check wasm64 support
    try {
        new WebAssembly.Memory({ address: 'i64', initial: 1n });
    } catch (error) {
        alert("Wasm64 is not supported. Please use a different browser...");
    }

    alphaChange();
}

startup();
