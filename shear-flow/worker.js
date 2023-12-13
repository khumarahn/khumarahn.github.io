"use strict";

var Module = {
    'onRuntimeInitialized': runOrders
};   

importScripts('20.js');

let orders = [];

onmessage = function(e) {
    if (orders.length == 0) {
        orders.push(e.data);
    } else {
        for (let j=0; j<orders.length; j++) {
            if (e.data.purpose == orders[j].purpose) {
                orders[j] = e.data;
            }
        }
    }
    //console.log(orders.length,JSON.stringify(orders));
}

let state = {};

function runOrders() {
    if (orders.length > 0) {
        let start_time = Date.now();
        let e = orders.shift();
        switch (e.purpose) {
            case 'init':
                let N_ = Module.table_N();
                let r_ = Module.table_particle(0).r;
                let col_rule_bottom = e.botColRule.charCodeAt(0);
                let col_rule_top = e.topColRule.charCodeAt(0);
                if (e.N == N_ && e.r == r_) {
                    Module.table_adjust(col_rule_bottom, col_rule_top, e.gravityX, e.gravityY, e.alphaBottom, e.betaBottom, e.alphaTop, e.betaTop);
                } else {
                    Module.table_init(e.N, e.r, col_rule_bottom, col_rule_top, e.gravityX, e.gravityY, e.alphaBottom, e.betaBottom, e.alphaTop, e.betaTop);
                }
                state.N = e.N;
                state.r = e.r;
                state.gravityX = e.gravityX;
                state.gravityY = e.gravityY;
                state.alphaBottom = e.alphaBottom;
                state.betaBottom = e.betaBottom;
                state.alphaTop = e.alphaTop;
                state.betaTop = e.betaTop;
                break;
            case 'live':
                Module.table_live(e.time);
                break;
        }
        state.particles = [];
        for (let j=0; j<state.N; j++) {
            let p = Module.table_particle(j);
            state.particles.push({
                x: p.x,
                y: p.y,
                r: p.r,
                vx: p.vx,
                vy: p.vy,
                m: p.m
            });
        }
        state.hist_po = [];
        for (let j=0; j<Module.hist_po_N(); j++) {
            let p = Module.hist_po(j);
            state.hist_po.push({
                x: p.x,
                y: p.y,
                r: p.r,
                vx: p.vx,
                vy: p.vy,
                m: p.m
            });
        }
        state.time = Module.table_time();
        state.num_col_pp = Module.table_num_col_pp();
        state.num_col_po = Module.table_num_col_po();
        state.num_col    = Module.table_num_col();

        state.processing_time = Date.now() - start_time;
        state.log = Module.table_log();
        postMessage(state);
    }
    setTimeout(runOrders, 5);
}

