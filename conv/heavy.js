"use strict";

let reportDiv = document.getElementById("reportDiv");
let last_report_timestamp = 0;

let f_table = document.querySelector("#table");
let f_table2d_pre = f_table.getContext("bitmaprenderer");

let W = f_table.width;
let H = f_table.height;

let table = new OffscreenCanvas(W, H);
table.width = W;
table.height = H;
let table2d = table.getContext("2d");

let sums = [{n:0, s:0}];
let sums_N = 1;

let ND = 256;
let maxT = 7/32;

function randn_bm(sigma2) {
        let u = 0, v = 0;
        while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
        while(v === 0) v = Math.random();
        return Math.sqrt( -2.0 * Math.log( u ) * sigma2 ) * Math.cos( 2.0 * Math.PI * v );
}

function update_sums(N) {
    let new_sums  = [sums[0]];
    let max_deque = [sums[0]];
    let min_deque = [sums[0]];
    let minmax = "min";
    let dN = Math.max(1, Math.floor(N/W / 2));
    let max_left = dN+1;
    let min_left = dN+1;
    
    function push(x) {
        let tt = [];
        let ttn = 0;
        
        // min
        while(min_deque.length > 0 && min_deque[min_deque.length - 1].n + dN > x.n && x.s <= min_deque[min_deque.length - 1].s) {
            min_deque.pop();
        }
        min_deque.push(x);
        ttn = new_sums[new_sums.length - 1].n;
        while (min_deque.length > 1 && min_deque[0].n + dN <= x.n) {
            if (min_left >= dN || ttn + dN <= min_deque[0].n) {
                tt.push(min_deque[0]);
                ttn = min_deque[0].n;
            }
            min_left = min_deque[1].n - min_deque[0].n;
            min_deque.shift();
        }
        
        // max
        while(max_deque.length > 0 && max_deque[max_deque.length - 1].n + dN > x.n && x.s >= max_deque[max_deque.length - 1].s) {
            max_deque.pop();
        }
        max_deque.push(x);
        ttn = new_sums[new_sums.length - 1].n;
        while (max_deque.length > 1 && max_deque[0].n + dN <= x.n) {
            if (max_left >= dN || ttn + dN <= max_deque[0].n) {
                tt.push(max_deque[0]);
                ttn = max_deque[0].n;
            }
            max_left = max_deque[1].n - max_deque[0].n;
            max_deque.shift();
        }
        
        tt.sort(function(a, b) { return a.n < b.n ? -1 : a.n > b.n ? 1 : 0; });
        for (let t = 0; t < tt.length; t++) {
            if (new_sums[new_sums.length - 1].n < tt[t].n) {
                new_sums.push(tt[t]);
            }
        }
    } // push

    function flush() {
        let tt = [min_deque[0], max_deque[0], min_deque[min_deque.length - 1]];
        tt.sort(function(a, b) { return a.n < b.n ? -1 : a.n > b.n ? 1 : 0; });
        for (let t = 0; t < tt.length; t++) {
            if (new_sums[new_sums.length - 1].n < tt[t].n) {
                new_sums.push(tt[t]);
            }
        }
    } // flush

    for (let t = 0; t < sums.length; t++) {
        push(sums[t]);
    }
    
    let s = sums[sums.length - 1].s;
    for (let n = sums_N; n < N; n++) {
        let r = 0;
        while (r == 0) r = Math.random();
        s += Math.pow(r, -3./4.) - 4;
      push({n: n, s: s});
    }
    flush();
    sums = new_sums;
    sums_N = N;
}

function draw(timestamp) {

    // draw the current state
    let table_bitmap = table.transferToImageBitmap();
    f_table2d_pre.transferFromImageBitmap(table_bitmap);
    
    // report n
    if (timestamp > last_report_timestamp + 2000) {
        reportDiv.innerHTML = "n: " + sums_N.toLocaleString("en");
        last_report_timestamp = timestamp;
    }

    // prepare the next
    update_sums(sums_N + 1 + Math.min(50000, 
        Math.floor(Math.pow(sums_N / 10, 1.1) * 0.01)
    ));
    table2d.clearRect(0,0,W,H);
    
    {   // axes
        table2d.strokeStyle = '#AAAAAA';
        table2d.lineWidth = 2;

        table2d.beginPath();
        table2d.moveTo(0,0.5*H);
        table2d.lineTo(W,0.5*H);
        table2d.stroke();
        table2d.closePath();
    }

    let s_abs_max = Math.max.apply(null, sums.map(function (x) { return Math.abs(x.s); }));
    let norm_s = 0.125 * 0.5 * Math.pow(sums_N, -1/2);
    let extra_norm_s = Math.min(1, 1. / (norm_s * 1.02 * s_abs_max));
    norm_s *= extra_norm_s;
    { // n^{?}
        table2d.strokeStyle = '#CCCCCC';
        table2d.lineWidth = 1;
        
        for (let s = -2; s <= 2; s += 0.5) {
            if (s == 0) {
                continue;
            }
            table2d.beginPath();
            table2d.moveTo(0,0.5*H);
            for (let t = 0.; t <= 1.; t += 1./128.) {
                let x = W * t;
                let y = H * (1. - s * 1 * Math.pow(t, 3./4.) * extra_norm_s) / 2;
                table2d.lineTo(x, y);
            }
            table2d.stroke();
            table2d.closePath();
        }
    }

    table2d.strokeStyle = '#F00000';
    table2d.lineWidth = 1;
    
    table2d.beginPath();
    let norm_t = 1. / sums_N;
    for (let p=0; p<sums.length; p++) {
        let t = W * sums[p].n * norm_t;
        let t1 = W * (sums[p].n + 1) * norm_t;
        let s = H * (1 - sums[p].s * norm_s) / 2;
        if (p==0) {
            table2d.moveTo(t, s);
        }
        if (p + 1 < sums.length && sums[p+1].n == sums[p].n + 1) {
            table2d.lineTo(t, s);
            table2d.lineTo(t1, s);
        } else {
            table2d.lineTo(t, s);
        }
    }
    table2d.stroke();
    table2d.closePath();

    if (sums_N < 1000) {
        requestAnimationFrame(() => {
            setTimeout(draw, Math.min(100, 10000/sums_N));
        });
    } else {
        requestAnimationFrame(draw);
    }
}

draw(0);
