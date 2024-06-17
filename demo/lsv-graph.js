let alpha = 0.75;

function lsv_trace() {
    trace = {
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
    alpha = document.getElementById("alphaSlider").value / 64;
    document.getElementById("alphaValue").innerHTML = alpha.toFixed(4);
    Plotly.newPlot('lsv', [lsv_trace()]);
    
    let v = new Function('x', 'return ' + document.getElementById("v").value);

    Plotly.newPlot('Ln', computeLv(v,10), { yaxis: {range: [0.0, 2.0]}, xaxis: {dtick: 0.125}});
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
        for (x=0.0; x <=1.0; x+=1./128) {
            trace.x.push(x);
            trace.y.push(LSV_Ln(v, x, alpha, k));
        }
        traces.push(trace);
    }
    return traces;
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
