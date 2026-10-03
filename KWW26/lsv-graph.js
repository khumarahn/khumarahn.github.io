"use strict";

// alpha as a BigInt
const ALPHA_DIGITS = 54,
      ALPHA_UNIT = 10n ** BigInt(ALPHA_DIGITS),
      ALPHA_MIN = ALPHA_UNIT / 2n,          // 0.5
      ALPHA_MAX = ALPHA_UNIT * 3n / 2n,     // 1.5
      ALPHA_SLIDER = ALPHA_UNIT / 256n,     // one tick of the slider
      ALPHA_WHEEL = ALPHA_UNIT / 10n ** BigInt(ALPHA_DIGITS - 3);

let alpha = ALPHA_MIN,
    alphaTimer;

function alphaString(a) {
    return `${a / ALPHA_UNIT}.${(a % ALPHA_UNIT).toString().padStart(ALPHA_DIGITS, "0")}`
        .replace(/\.?0+$/, "");
}

function alphaParse(s) {
    let m = /^([01])(?:\.(\d*))?$/.exec(s.trim());
    if (!m) {
        return null;
    }
    let digits = (m[2] ?? "").padEnd(ALPHA_DIGITS + 1, "0"),   // past the grid, to round on
        a = BigInt(m[1]) * ALPHA_UNIT + BigInt(digits.slice(0, ALPHA_DIGITS))
            + (digits[ALPHA_DIGITS] > "4" ? 1n : 0n);
    return (a < ALPHA_MIN || a > ALPHA_MAX) ? null : a;
}

function lsv_trace() {
    let a = Number(alpha) / Number(ALPHA_UNIT); // the picture needs no exactness
    let trace = {
        x: [],
        y: [],
        type: 'scatter'
    };

    for (let x = 0; x <= 0.5; x+=1./128.) {
        trace.x.push(x);
        trace.y.push(LSV_left(x, a));
    }
    trace.x.push(0.5);
    trace.y.push(NaN);
    for (let x = 0.5; x <= 1.0; x+=1./4.) {
        trace.x.push(x);
        trace.y.push(LSV_right(x, a));
    }

    return trace;
}

function alphaSet(a) {
    alpha = a < ALPHA_MIN ? ALPHA_MIN : (a > ALPHA_MAX ? ALPHA_MAX : a);

    document.getElementById("alphaValue").innerHTML = alphaString(alpha);
    document.getElementById("alphaSlider").value = Number((alpha - ALPHA_MIN) / ALPHA_SLIDER);

    clearTimeout(alphaTimer);
    alphaTimer = setTimeout(alphaChange, 300);
}

function alphaInput() {
    alphaSet(ALPHA_MIN + BigInt(document.getElementById("alphaSlider").value) * ALPHA_SLIDER);
}

function alphaChange() {
    clearTimeout(alphaTimer);
    Plotly.newPlot('lsv', [lsv_trace()], {yaxis: {range: [0,1.02], dtick: 0.25}, xaxis: {range: [0,1.02], dtick: 0.25}});

    if (!workerBusy) {
        prepareWorker();
    }

    MathJax.typeset();
}

function alphaEdit() {
    clearTimeout(alphaTimer);
    let el = document.getElementById("alphaValue");
    el.innerHTML = `<input size="${ALPHA_DIGITS + 2}" value="${alphaString(alpha)}">`;

    let inp = el.firstChild;
    inp.focus();
    inp.select();

    inp.onblur = function() {
        inp.onblur = null;
        alphaSet(alphaParse(inp.value) ?? alpha);
    };
    inp.onkeydown = function(e) {
        if (e.key === "Enter") {
            inp.blur();
        } else if (e.key === "Escape") {
            inp.onblur = null;
            alphaSet(alpha);
        }
    };
}

document.getElementById("alphaValue").addEventListener("wheel", function(e) {
    e.preventDefault();
    let notch = (alpha + ALPHA_WHEEL / 2n) / ALPHA_WHEEL + (e.deltaY < 0 ? 1n : -1n);
    alphaSet(notch * ALPHA_WHEEL);
}, { passive: false });

document.getElementById("alphaValue").addEventListener("click", function(e) {
    if (e.target === e.currentTarget) {
        alphaEdit();
    }
});

function prepareWorker() {
    document.getElementById('theorem-control').innerHTML = `
            <br>
                <button id="theorem-btn" style=" background-color: #007bff; color: white; border: none; 
                    border-radius: 0.3rem; padding: 0.5rem 1rem; font-size: 1rem; cursor: pointer;  margin-top: 0.5rem;">
                    Formulate the result for \\(\\alpha=${alphaString(alpha)}\\)
                </button>
            `;

    document.getElementById('theorem-btn').addEventListener('click', function(e) {
        e.preventDefault(); // no page jump

        document.getElementById('theorem-control').innerHTML =
            `<span>Computing bounds for \\(\\alpha \\approx ${alphaString(alpha)}\\)...</span>`;

        document.getElementById('reasoning-log').textContent += '\n\n';
        document.getElementById('reasoning-summary').textContent = 'Reasoning...';

        worker.postMessage({
            type: 'compute-bounds',
            alpha_num: alpha.toString(),
            alpha_den: ALPHA_UNIT.toString()
        });
        workerBusy = true;

        MathJax.typeset();
    });
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
        For all
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

        thm = thm.replaceAll('::alpha_minus::', m.alpha_minus)
            .replaceAll('::alpha_plus::', m.alpha_plus)
            .replaceAll('::hp::', m.min_hp_h_prime)
            .replaceAll('::hpp::', m.max_hpp_h_prime);

        if (parseFloat(m.gamma_plus) < 1) {
            let extra = String.raw`
            <div style="margin-top: 1em;">
            Average of the first return time \(\tau (x) = \inf \{k \geq 1 : T^k(x) \in [1/2,1]\) is:
            \[
                \frac{\int_{1/2}^1 \tau(x) h(x) \, dx}{\int_{1/2}^1 h(x) \, dx}
                \in [::tau_minus::, ::tau_plus::]
                ,
            \]
            and the Lyapunov exponent is:
            \[
                \frac{\int_0^1 \log (T'(x)) h(x) \, dx}{\int_0^1 h(x) \, dx}
                \in [::lambda_minus::, ::lambda_plus::]
                .
            \]
            </div>
            `;
            extra = extra.replaceAll('::tau_minus::', m.tau_minus)
                         .replaceAll('::tau_plus::', m.tau_plus)
                         .replaceAll('::lambda_minus::', m.lambda_minus)
                         .replaceAll('::lambda_plus::', m.lambda_plus);
            thm += extra;
        }

        document.getElementById("theorem-text").innerHTML += thm;

        document.getElementById('theorem-control').innerHTML = '';
        document.getElementById('reasoning-summary').innerHTML = 'Done reasoning';

        workerBusy = false;
        prepareWorker();

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

    alphaInput();
}

startup();
