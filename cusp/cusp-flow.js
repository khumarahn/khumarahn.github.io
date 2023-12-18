"use strict";

let configChanged = true;

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
let f_table_pre = f_table.getContext("bitmaprenderer");

let W = f_table.width;
let H = f_table.height;

let table = new OffscreenCanvas(W, H);
let table2d = table.getContext("2d");
let table_bitmap = table.transferToImageBitmap();

let table_bg = new OffscreenCanvas(W, H);
let table_bg2d = table_bg.getContext("2d", { willReadFrequently: true });
let table_bg2d_image = table_bg2d.getImageData(0,0,W,H);

let R = 5;
let p_canvas_image = new Image();
{
    let p_canvas = document.createElement('canvas');
    let p_canvas2d = p_canvas.getContext("2d");
    p_canvas2d.fillStyle = 'rgba(200,0,0,0.8)';
    p_canvas2d.clearRect(0,0, 2*R + 1, 2*R + 1);
    p_canvas2d.beginPath();
    p_canvas2d.arc(R, R, R, 0, Math.PI*2);
    p_canvas2d.closePath();
    p_canvas2d.fill();
    p_canvas_image.src = p_canvas.toDataURL();
}

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

    table_bg2d_image = table_bg2d.getImageData(0,0,W,H);
}

function drawParticles(timestamp) {
    if (configChanged) {
        drawTable();
        
        let my = alpha * Math.pow(p.x, beta);
        while (Math.abs(p.y) > my) {
            p.y = Math.sign(p.y) * my;
        }

        configChanged = false;
    } else {
        // draw the current state
        f_table_pre.transferFromImageBitmap(table_bitmap);

        // prepare the next
        p = CUSPt(p, alpha, beta, speed / 32);

        table2d.putImageData(table_bg2d_image, 0, 0);
        table2d.drawImage(p_canvas_image, tx(p.x) - R, ty(p.y) - R);

        table_bitmap = table.transferToImageBitmap();
    }

    requestAnimationFrame(drawParticles);
}

window.onresize = resize;

function resize() {
    W = Math.min(document.body.clientWidth, document.body.clientHeight) - 10;
    H = W;
    f_table.width = W;
    f_table.height = H;
    table.width = W;
    table.height = H;
    table_bg.width = W;
    table_bg.height = H;

    configChanged = true;
}

resize();

drawParticles(0);
