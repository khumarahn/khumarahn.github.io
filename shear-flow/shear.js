"use strict";

let f_table = document.querySelector("#table");
let f_table2d = f_table.getContext("2d");

let graph_ctx = document.getElementById('graph');
let graph_chart = new Chart(graph_ctx, {
    type: 'line',
    data: {
        datasets: [
            {
                label: 'VX',
                lineTension: 0,
                borderColor: 'rgba(255, 142, 0, 0.842)',
                fill: false,
                data: [{x:0.0, y:0.0}, {x:0.5,y:1.0}, {x:1.0, y:0.0}]
            },
            {
                label: 'TX+TY',
                lineTension: 0,
                borderColor: 'rgba(255, 15, 15, 0.842)',
                fill: false,
                data: [{x:0.0, y:0.0}, {x:0.5,y:1.0}, {x:1.0, y:0.0}]
            },
            {
                label: 'DENSITY',
                lineTension: 0,
                borderColor: 'rgba(15, 120, 255, 0.842)',
                fill: false,
                data: [{x:0.0, y:0.0}, {x:0.5,y:1.0}, {x:1.0, y:0.0}]
            },
            {
                label: 'TX',
                lineTension: 0,
                borderColor: 'rgba(125, 15, 15, 0.842)',
                fill: false,
                data: [{x:0.0, y:0.0}, {x:0.5,y:1.0}, {x:1.0, y:0.0}]
            },
            {
                label: 'TY',
                lineTension: 0,
                borderColor: 'rgba(15, 125, 15, 0.842)',
                fill: false,
                data: [{x:0.0, y:0.0}, {x:0.5,y:1.0}, {x:1.0, y:0.0}]
            },
            {
                label: 'VX_SIGMA2',
                lineTension: 0,
                borderColor: 'rgba(0, 0, 0, 0.842)',
                fill: false,
                data: [{x:0.0, y:0.0}, {x:0.5,y:1.0}, {x:1.0, y:0.0}]
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
        }
    }
});
let graph_last_update = 0;

let W = f_table.width;
let H = f_table.height;

let table = new OffscreenCanvas(W, H);
let table2d = table.getContext("2d");

let recalibrate = true;

let zN;
let realN;
let zNSlider = document.getElementById("zNSlider");
let zNValue = document.getElementById("zNValue");
let r;
let rSlider = document.getElementById("rSlider");
let rValue = document.getElementById("rValue");
function r_zN_oninput() {
    zN = Number(zNSlider.value);
    realN = Math.pow(2, (zN % 4)) * Math.pow(10, Math.floor(zN / 4));
    zNValue.innerHTML = "Number of particles: " + realN;

    let volumeFraction = Number(rSlider.value) / 100.;
    r = Math.sqrt(0.5 * volumeFraction / (Math.PI * realN));
    rValue.innerHTML = "Particle volume density: " + volumeFraction.toFixed(4) + " (radius = " + r.toExponential(4) +
        ", L = " + (0.5 / (2*r)).toFixed(4) + ")";
    recalibrate = true;
}
zNSlider.oninput = r_zN_oninput;
rSlider.oninput =  r_zN_oninput;

r_zN_oninput();

let speed;
let speedSlider = document.getElementById("speedSlider");
let speedValue = document.getElementById("speedValue");
speedSlider.oninput = function() {
    speed = Math.exp((this.value-256)/256*Math.log(100));
    speedValue.innerHTML = "Maximum simulation speed: " + speed.toFixed(2);
}
speedSlider.oninput();

let gravityX, gravityY;
let gravitySelect = document.getElementById("gravitySelect");
let gravityDiv = document.getElementById("gravityDiv");
let gravitySpan = document.getElementById("gravitySpan");
let gravityXSlider = document.getElementById("gravityXSlider");
let gravityYSlider = document.getElementById("gravityYSlider");

function setGravity() {
    if (gravitySelect.value == "off") {
        gravityX = 0.0;
        gravityY = 0.0;
        gravityDiv.style.display = "none";
        gravitySpan.style.display = "none";
    } else {
        gravityX = gravityXSlider.value / 32;
        gravityY = gravityYSlider.value / 32;
        gravityDiv.style.display = "";
        gravitySpan.style.display = "";
    }
    gravitySpan.innerHTML = "value: (" + gravityX.toFixed(2) + ", " + gravityY.toFixed(2) + ")";

    recalibrate = true;
}
gravitySelect.oninput = setGravity;
gravityXSlider.oninput = setGravity;
gravityYSlider.oninput = setGravity;
setGravity();


// top wall
let topColRule, topWallAlpha, topWallBeta;
let topWallAlphaSlider  = document.getElementById("topWallAlphaSlider");
let topWallAlphaDiv     = document.getElementById("topWallAlphaDiv");
let topWallAlphaSpan    = document.getElementById("topWallAlphaSpan");
let topWallBetaSlider   = document.getElementById("topWallBetaSlider");
let topWallBetaDiv      = document.getElementById("topWallBetaDiv");
let topWallBetaSpan     = document.getElementById("topWallBetaSpan");
let topColRuleSelect = document.getElementById("topColRuleSelect");
// bottom wall
let botColRule, botWallAlpha, botWallBeta;
let botWallAlphaSlider  = document.getElementById("botWallAlphaSlider");
let botWallAlphaDiv     = document.getElementById("botWallAlphaDiv");
let botWallAlphaSpan    = document.getElementById("botWallAlphaSpan");
let botWallBetaSlider   = document.getElementById("botWallBetaSlider");
let botWallBetaDiv      = document.getElementById("botWallBetaDiv");
let botWallBetaSpan     = document.getElementById("botWallBetaSpan");
let botColRuleSelect = document.getElementById("botColRuleSelect");

function colParam(rule, alphaSlider, betaSlider) {
    let alpha = 0;
    let beta = 0;
    let alphaVis = "none";
    let alphaSpan = "";
    let betaVis = "none";
    let betaSpan = "";
    if (rule == 'a') {
        alpha = Math.PI * Math.max(0.0001, Math.min(0.9999, 0.5 + alphaSlider / 2 / 256));
        alphaVis =  "";
        alphaSpan = "α: " + alpha.toFixed(4) + " = " + (alpha / Math.PI).toFixed(4) + "π";
    } else if (rule == 'b') {
        alpha = (alphaSlider + 256) * 8 / 64;
        alphaVis = "";
        alphaSpan = "b: " + alpha.toFixed(4);
    } else if (rule == 'c') {
        alpha = (alphaSlider + 256) / 512;
        alphaVis = "";
        alphaSpan = "c: " + alpha.toFixed(4);
    } else if (rule == 'd') {
        alpha = Math.PI * Math.max(0.0001, Math.min(0.9999, 0.5 + alphaSlider / 2 / 256));
        alphaVis = "";
        alphaSpan = "α: " + alpha.toFixed(4) + " = " + (alpha / Math.PI).toFixed(4) + "π";
    } else if (rule == 'm') {
        alpha = alphaSlider  / 128;
        beta  = (betaSlider + 256) * 2 / 512.;
        alphaVis = "";
        alphaSpan = "α: " + alpha.toFixed(4);
        betaVis = "";
        betaSpan = "T<sub>w</sub>: " + beta.toFixed(4);
    } else if (rule == 'n') {
        alpha = alphaSlider * 8 / 256;
        alphaVis = "";
        alphaSpan = "bias: " + alpha.toFixed(4);
    }
    return {
        alpha: alpha,
        beta: beta,
        alphaVis: alphaVis,
        betaVis: betaVis,
        alphaSpan: alphaSpan,
        betaSpan: betaSpan
    };
}

function colChange() {
    topColRule = topColRuleSelect.value;
    let topColParam = colParam(topColRule, topWallAlphaSlider.valueAsNumber, topWallBetaSlider.valueAsNumber);
    topWallAlpha = topColParam.alpha;
    topWallBeta  = topColParam.beta;
    botColRule = (botColRuleSelect.value == "-") ? topColRule : botColRuleSelect.value;
    let botColParam = (botColRuleSelect.value == "-") ?
        topColParam : colParam(botColRule, botWallAlphaSlider.valueAsNumber, botWallBetaSlider.valueAsNumber);
    botWallAlpha = botColParam.alpha;
    botWallBeta  = botColParam.beta;

    topWallAlphaDiv.style.display = topColParam.alphaVis;
    topWallAlphaSpan.innerHTML = topColParam.alphaSpan;
    topWallBetaDiv.style.display = topColParam.betaVis;
    topWallBetaSpan.innerHTML = topColParam.betaSpan;

    if (botColRuleSelect.value == "-") {
        botWallAlphaDiv.style.display = "none";
        botWallAlphaSpan.innerHTML = "";
        botWallBetaDiv.style.display = "none";
        botWallBetaSpan.innerHTML = "";
    } else {
        botWallAlphaDiv.style.display = botColParam.alphaVis;
        botWallAlphaSpan.innerHTML = botColParam.alphaSpan;
        botWallBetaDiv.style.display = botColParam.betaVis;
        botWallBetaSpan.innerHTML = botColParam.betaSpan;
    }

    recalibrate = true;
}

topColRuleSelect.onchange = colChange;
topWallAlphaSlider.oninput = colChange;
topWallBetaSlider.oninput = colChange;
botColRuleSelect.onchange = colChange;
botWallAlphaSlider.oninput = colChange;
botWallBetaSlider.oninput = colChange;

colChange();


let observationsDiv = document.getElementById("observationsDiv");
let averagePointsSpan = document.getElementById("averagePointsSpan");
let averagePointsDiv  = document.getElementById("averagePointsDiv");

let observations;
function resetObservations() {
    let averagePoints = Math.pow(2, averagePointsSlider.valueAsNumber);
    averagePointsSpan.innerHTML = "Approximate number of points for computing averages: " + averagePoints.toString();

    observations = {
        averagePoints: averagePoints,
        // mean free path between collisions is 1 / (8 N r), assuming mean energy 1, observe every:
        frequency : 1. / (realN * r),
        last_time: null,
        time: 0.,
        h_N: 32 + 1,     // vertical is divided into ... "bins"
        energy: 0,
        vel: [],
        vel_hist_max: 4 * 4096,
        hist_po: [],
        hist_po_cur : 0,
        hist_po_max: 4 * 4 * 4096,
        num_col_pp: 0,
        num_col_po: 0,
        num_col: 0
    };
    observationsDiv.innerHTML = "";
}
averagePointsSlider.onchange = resetObservations;
resetObservations();

// bin number for a given y: we place observations.h_N bins in [0.25+r, 0.75-r]
function ybin(y) {
    let bin = Math.floor((y - (0.25 + r)) * observations.h_N / ((0.75 - r) - (0.25 + r)));
    if (bin >= observations.h_N) {
        return observations.h_N - 1;
    } else if (bin < 0) {
        return 0;
    } else {
        return bin;
    }
}
// center of the bin with given number
function ybin_center(b) {
    return (0.25 + r) + (Number(b) + 0.5) * ((0.75 - r) - (0.25 + r)) / observations.h_N;
}

let state = null;

function SQR(x) {
    return x*x;
}

function initTable() {
    worker.postMessage({
        purpose: 'init',
        N: realN,
        r: r,
        topColRule: topColRule,
        botColRule: botColRule,
        gravityX: gravityX,
        gravityY: gravityY,
        alphaBottom: botWallAlpha,
        betaBottom:  botWallBeta,
        alphaTop:    topWallAlpha,
        betaTop:     topWallBeta
    });
}

let prev_timestamp = null;
let prev_dt = null;

function drawParticles(timestamp) {

    if (recalibrate) {
        recalibrate = false;
        prev_dt = null;
        initTable();

        resetObservations();
        return;
    }

    // draw the current state
    f_table2d.drawImage(table,0,0);

    // prepare the next state
    let dt = 0;
    if (prev_timestamp != null) {
        let dtstamp = timestamp - prev_timestamp;
        dt = dtstamp / 1000 * speed / 4;
        if (prev_dt != null) {
            // we did prev_dt in state.processing_time milliseconds
            // should we adjust speed so the next frame is computed in like 20ms or less?
            let dt_20ms = Math.max(prev_dt, 0.00000000001) * 20 / Math.max(state.processing_time, 1);
            dt = Math.min(dt_20ms, prev_dt * 1.03125, dt);
        } else {
            dt = 0.00000001 * dt;
        }
        prev_dt = dt;
    }
    prev_timestamp = timestamp;

    let dtt = state.time + dt;
    if (observations.last_time != null) {
        dtt = Math.min(dtt, observations.last_time + observations.frequency);
    }

    worker.postMessage({
        purpose: 'live',
        time: state.time + dt
    });

    table2d.clearRect(0,0,W,H);

    table2d.fillStyle = "#EEEEEE";
    table2d.fillRect(0,0,W,H);


    function drawDisk(x,y,r) {
        table2d.beginPath();
        table2d.arc(x, y, r, 0, Math.PI*2);
        table2d.closePath();
        table2d.fill();
    }

    function drawArrow(fromx, fromy, tox, toy) {
        let headlen = Math.min(5, 0.25 * Math.sqrt((tox - fromx)*(tox-fromx) + (toy - fromy)*(toy - fromy)));   // length of head in pixels
        let angle = Math.atan2(toy-fromy,tox-fromx);
        table2d.moveTo(fromx, fromy);
        table2d.lineTo(tox, toy);
        table2d.lineTo(tox-headlen*Math.cos(angle-Math.PI/6),toy-headlen*Math.sin(angle-Math.PI/6));
        table2d.moveTo(tox, toy);
        table2d.lineTo(tox-headlen*Math.cos(angle+Math.PI/6),toy-headlen*Math.sin(angle+Math.PI/6));
    }

    function tx(x) {
        return W * x;
    }

    function ty(y) {
        return H * 2 * (0.75 - y);
    }

    if (animateCheckbox.checked) {
        for (let j=0; j<state.N; j++) {
            if (j == 0) {
                table2d.fillStyle = 'rgba(200,0,0,0.5)';
            } else if (j == 1) {
                table2d.fillStyle = 'rgba(0,200,0,0.5)';
            }
            let p = state.particles[j];
            let x = tx(p.x);
            let y = ty(p.y);
            let r = W * p.r;
            drawDisk(x,y,Math.max(1,r));
            if (p.x + p.r > 1) {
                drawDisk(x-W, y, r);
            } else if (p.x - p.r < 0) {
                drawDisk(x+W, y, r);
            }
        }
    }

    // observe
    // some observations happen every step
    // wall collisions
    {
        for (let k in state.hist_po) {
            let p = state.hist_po[k];
            if (observations.hist_po_cur == observations.hist_po.length) {
                observations.hist_po.push([p.vx, p.vy]);
            } else {
                observations.hist_po[observations.hist_po_cur] = [p.vx, p.vy];
            }
            observations.hist_po_cur = (observations.hist_po_cur + 1)
                % observations.hist_po_max;
        }
    }
    // time and energy
    {
        observations.energy = 0.;
        for (let j=0; j<state.N; j++) {
            observations.time = state.time;
            let p = state.particles[j];
            observations.energy += 0.5 * p.m * (SQR(p.vx) + SQR(p.vy));
        }
    }
    // collision count
    {
        observations.num_col_pp = state.num_col_pp;
        observations.num_col_po = state.num_col_po;
        observations.num_col    = state.num_col;
    }

    // some observations happen sometimes
    if (observations.last_time == null || state.time + 0.0000001 >= observations.last_time + observations.frequency) {
        let avx = [];
        for (let j=0; j<state.N; j++) {
            let p = state.particles[j];
            let yy = ybin(p.y);

            // average velocity in the strip at the time
            {
                if (typeof avx[yy] === 'undefined') {
                    avx[yy] = {
                        s: 0,
                        n: 0
                    }
                }
                avx[yy].s += p.vx;
                avx[yy].n++;
            }

            // velocities and related
            {
                if (typeof observations.vel[yy] === 'undefined') {
                    observations.vel[yy] = {
                        vx: 0,
                        vy: 0,
                        sigma2_x: 0,
                        sigma2_y: 0,
                        n: 0,
                        hist: [],
                        hist_cur : 0,
                        avx: 0,
                        avx_sigma2: 0,
                        avx_n: 0
                    };
                }
                let d = Math.min(observations.averagePoints, observations.vel[yy].n + 1);
                observations.vel[yy].sigma2_x = (d-1)/d * observations.vel[yy].sigma2_x + 1/d * SQR(p.vx - observations.vel[yy].vx);
                observations.vel[yy].sigma2_y = (d-1)/d * observations.vel[yy].sigma2_y + 1/d * SQR(p.vy - observations.vel[yy].vy);
                observations.vel[yy].vx = (d-1)/d * observations.vel[yy].vx + 1/d * p.vx;
                observations.vel[yy].vy = (d-1)/d * observations.vel[yy].vy + 1/d * p.vy;
                observations.vel[yy].n ++;

                if (observations.vel[yy].hist_cur == observations.vel[yy].hist.length) {
                    observations.vel[yy].hist.push([p.vx, p.vy]);
                } else {
                    observations.vel[yy].hist[observations.vel[yy].hist_cur] = [p.vx, p.vy];
                }
                observations.vel[yy].hist_cur ++;
                if (observations.vel[yy].hist_cur == observations.vel_hist_max) {
                    observations.vel[yy].hist_cur = 0;
                }
            }
        }

        // covariance of average velocity in strips
        {
            for (let yy in avx) {
                let d = Math.min(observations.averagePoints, observations.vel[yy].avx_n + 1);
                let a = avx[yy].s / avx[yy].n;
                observations.vel[yy].avx_sigma2 = (d-1)/d * observations.vel[yy].avx_sigma2 + 1/d * SQR(a - observations.vel[yy].avx);
                observations.vel[yy].avx        = (d-1)/d * observations.vel[yy].avx        + 1/d * a;
                observations.vel[yy].avx_n++;
            }
        }

        observations.last_time = state.time;
    }

    // occasionally we display the observations
    let seconds = new Date().getTime() / 1000;
    if (graph_last_update + 5 < seconds) {
        graph_last_update = seconds;

        // text
        {
            observationsDiv.innerHTML = "Time: " + observations.time.toFixed(4)
                + ", total energy: " + observations.energy.toFixed(6)
                + ",<br />collision count: " + observations.num_col_po.toLocaleString()
                + " with walls and " + observations.num_col_pp.toLocaleString() + " between particles"
                + "." ;
        }

        // temp and drift and density etc
        {
            let aN = 0;
            for (let yy in observations.vel) {
                aN += observations.vel[yy].n;
            }
            for (let yy in observations.vel) {
                let y0 = ybin_center(yy);
            }

            let data_vfield = [];
            let data_temp = [];
            let data_density = [];
            let data_temp_x = [];
            let data_temp_y = [];
            let data_avx_sigma2 = [];
            let data_avx_sigma2_norm;
            for (let yy in observations.vel) {
                let y0 = ybin_center(yy);
                let dx = observations.vel[yy].vx;
                let dy = observations.vel[yy].vy;
                let Tx = observations.vel[yy].sigma2_x;
                let Ty = observations.vel[yy].sigma2_y;
                let T = Tx + Ty;
                let n = observations.vel[yy].n;
                let density = n / aN * observations.h_N;
                let avx_sigma2 = observations.vel[yy].avx_sigma2;
                let avx_sigma2_norm = avx_sigma2 * n / aN * state.N;
                data_vfield.push({y: y0-0.5, x: dx});
                data_temp_x.push({y: y0-0.5, x: Tx});
                data_temp_y.push({y: y0-0.5, x: Ty});
                data_temp.push({y: y0-0.5, x: Tx + Ty});
                data_density.push({y: y0-0.5, x: density});
                data_avx_sigma2.push({y: y0-0.5, x: avx_sigma2_norm});
            }
            graph_chart.data.datasets[0].data  = data_vfield;
            graph_chart.data.datasets[1].data  = data_temp;
            graph_chart.data.datasets[2].data  = data_density;
            graph_chart.data.datasets[3].data  = data_temp_x;
            graph_chart.data.datasets[4].data  = data_temp_y;
            graph_chart.data.datasets[5].data  = data_avx_sigma2;
        }

        graph_chart.update();
    }

    // plot observations
    if (false)
    for (let yy in observations.vel) {
        let x0 = 0.5;
        let y0 = ybin_center(yy);
        let dx = observations.vel[yy].vx / 8;
        let dy = observations.vel[yy].vy / 8;
        table2d.beginPath();
        drawArrow(tx(x0), ty(y0), tx(x0+dx), ty(y0+dy));
        //table2d.moveTo(tx(x0), ty(y0));
        //table2d.lineTo(tx(x0 + dx), ty(y0 + dy));
        table2d.stroke();
    }
}

function downloadObservations() {
    let s = "strip_data <- list(";
    for (let yy in observations.vel) {
        let y0 = ybin_center(yy);
        let vel = observations.vel[yy];
        s += "list("
            + "y = " + y0.toFixed(3)
            + ", vel_hist = fromJSON('" + JSON.stringify(vel.hist) + "')"
            + ", vx = " + vel.vx.toFixed(5)
            + ", vy = " + vel.vy.toFixed(5)
            + ", sigma2_x = " + vel.sigma2_x.toFixed(5)
            + ", sigma2_y = " + vel.sigma2_y.toFixed(5)
            + "), "
    }
    s = s.substring(0, s.length - 2); // remote last ", "
    s += ");"; // end velocities list

    s += "\nhist_po <- fromJSON('" + JSON.stringify(observations.hist_po) + "');";

    // separately top and bottom walls from hist_po
    let hist_po_top = [];
    let hist_po_bottom = [];
    for (let k=0; k<observations.hist_po.length; k++) {
        let vx = observations.hist_po[k][0];
        let vy = observations.hist_po[k][1];
        if (vy > 0) {
            hist_po_top.push([vx, vy]);
        } else if (vy < 0) {
            hist_po_bottom.push([vx, vy]);
        }
    }
    s += "\nhist_po_top <- fromJSON('" + JSON.stringify(hist_po_top) + "');";
    s += "\nhist_po_bottom <- fromJSON('" + JSON.stringify(hist_po_bottom) + "');";

    downloadText("observations.txt", s);
}

let worker = new Worker('worker.js');
worker.onmessage = function (s) {
    state = s.data;
    //console.log("state log: " + state.log);
    requestAnimationFrame(drawParticles);
}

initTable();
