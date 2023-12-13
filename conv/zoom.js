"use strict";

let f_table = document.querySelector("#table");
let f_table2d_pre = f_table.getContext("bitmaprenderer");

let W = f_table.width;
let H = f_table.height;

let table = new OffscreenCanvas(W, H);
let table2d = table.getContext("2d");

let graph = [{t:0, x:0, y:0}, {t:0.01, x:0.01, y:0.01}]; // array of {t: ?, x: ?, y: }

let ND = 256;
let maxT = 7/32;

function randn_bm(sigma2) {
        let u = 0, v = 0;
        while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
        while(v === 0) v = Math.random();
        return Math.sqrt( -2.0 * Math.log( u ) * sigma2 ) * Math.cos( 2.0 * Math.PI * v );
}

function updateGraph(zoom) {
    
    // zoom
    let sqZoom = zoom * zoom;
    for (let p = 1; p < graph.length; p++) {
        graph[p].t *= sqZoom;
        graph[p].x *= zoom;
        graph[p].y *= zoom;
        if (graph[p].t > maxT) {
            graph.length = p;
            break;
        }
    }

    // fill in
    {
        let newGraph = [graph[0]];
        for (let p=1; p < graph.length; p++) {
            let dt = graph[p].t - graph[p-1].t;
            let dx = graph[p].x - graph[p-1].x;
            let dy = graph[p].y - graph[p-1].y;
            if (dt > 1 / ND || Math.abs(dx) > 1 / ND || Math.abs(dy) > 1 / ND) {
                let varr = dt / 4;
                newGraph.push({
                    t: graph[p-1].t + dt/2,
                    x: graph[p-1].x + dx/2 + randn_bm(varr),
                    y: graph[p-1].y + dy/2 + randn_bm(varr)
                });
            }
            newGraph.push(graph[p]);
        }
        graph = newGraph;
    }
}

function draw(timestamp) {

    // draw the current state
    let table_bitmap = table.transferToImageBitmap();
    f_table2d_pre.transferFromImageBitmap(table_bitmap);

    // prepare the next
    updateGraph(1.005);
    table2d.clearRect(0,0,W,H);
    
    {   // axes
        table2d.strokeStyle = '#AAAAAA';
        table2d.lineWidth = 1;

        table2d.beginPath();
        table2d.moveTo(0.5*W,0);
        table2d.lineTo(0.5*W,H);
        table2d.moveTo(0,0.5*H);
        table2d.lineTo(W,0.5*H);
        table2d.stroke();
        table2d.closePath();
    }

    table2d.strokeStyle = '#F00000';
    table2d.lineWidth = 1;
    
    table2d.beginPath();
    for (let p=0; p<graph.length; p++) {
        let x = W * (1 + graph[p].x) / 2;
        let y = H * (1 - graph[p].y) / 2;
        if (p==0) {
            table2d.moveTo(x, y);
        } else {
            table2d.lineTo(x, y);
        }
    }
    table2d.stroke();
    table2d.closePath();

    requestAnimationFrame(draw);
}

draw(0);
