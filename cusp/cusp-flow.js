"use strict";

let configChanged = false;

let speed;
let speedSlider = document.getElementById("speedSlider");
let speedValue = document.getElementById("speedValue");
speedSlider.oninput = function() {
    speed = Math.exp((this.value-256)/256*Math.log(100));
    speedValue.innerHTML = "Speed: " + speed.toFixed(2);
}
speedSlider.oninput();

let alpha = 0.5;
let beta;
let betaSlider = document.getElementById("betaSlider");
let betaValue = document.getElementById("betaValue");
betaSlider.oninput = function() {
    configChanged = true;
    beta = 1. + 3*Math.pow(this.value / 256, 4);
    betaValue.innerHTML = "beta: " + beta.toFixed(2);
}
betaSlider.oninput();

let f_table = document.querySelector("#table");
let f_table2d = f_table.getContext("2d");

let W = f_table.width;
let H = f_table.height;

let table = document.createElement('canvas');
table.width = W;
table.height = H;
let table2d = table.getContext("2d");

let table_bg = document.createElement('canvas');
table_bg.width = W;
table_bg.height = H;
let table_bg2d = table_bg.getContext("2d");

let p = {
    x: 0.5, 
    y: 0.0, 
    theta: 3 * Math.PI / 4
};
    
function tx(x) {
    return x * (W - 15);
}

function ty(y) {
    return 0.5 * H - y * (H - 15);
}

function drawTable() {
    table_bg2d.clearRect(0,0,W,H);

    table_bg2d.fillStyle = '#9aeaae';
    table_bg2d.fillRect(0,0,W,H);

    table_bg2d.fillStyle = '#FFFFFF';

    // obstables
    table_bg2d.beginPath();
    table_bg2d.moveTo(0,0.5*H);
    for (let x=0; x<=1; x+=1./256) {
        table_bg2d.lineTo(tx(x), ty(alpha*Math.pow(x, beta)));
    }

    for (let x=1; x>=0; x-=1./256) {
        table_bg2d.lineTo(tx(x), ty(-alpha*Math.pow(x, beta)));
    }
    table_bg2d.closePath();
    table_bg2d.fill();
}

function drawParticles(timestamp) {

    // draw the current state
    f_table2d.drawImage(table_bg,0,0);
    f_table2d.drawImage(table,0,0);

    // prepare the next
    if (configChanged) {
        drawTable();
        
        let my = alpha * Math.pow(p.x, beta);
        while (Math.abs(p.y) > my) {
            p.y = Math.sign(p.y) * my;
        }
    }
    p = CUSPt(p, alpha, beta, speed / 32);
    
    table2d.clearRect(0,0,W,H);

    table2d.fillStyle = 'rgba(200,0,0,0.8)';
        
    table2d.beginPath();
    table2d.arc(tx(p.x), ty(p.y), 5, 0, Math.PI*2);
    table2d.closePath();
    table2d.fill();

    requestAnimationFrame(drawParticles);
    //setTimeout(drawParticles,5000);
}

drawParticles(0);
