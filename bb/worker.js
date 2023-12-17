"use strict";

let module_init_done = false;
var Module = {
    'onRuntimeInitialized': function() {
        module_init_done = true;
        run_orders();
    }
};   
importScripts('32.js');


var bo; // billiard orbit
let table = -1; // table variant, negative is nothing

let orders = [];
let orders_running = true;

onmessage = function(e) {
    orders.push(e.data);
    //console.log("onmessage", orders.length,JSON.stringify(orders));
    if (module_init_done) {
        run_orders();
    }
}

function run_orders() {
    for (let o of orders) {
        if ("table" in o) {
            table = o.table;
            bo = new Module.orbit_c(table);
        }
        if (table >= 0 && "dt" in o) {
            let start_time = Date.now();
            bo.step(o.dt);
            let orbit = [];
            // <= size() because the extra point is the current position
            for (let i=0; i<=bo.size(); i++) {
                orbit.push({x: bo.xy(i).x, y: bo.xy(i).y});
            }
            postMessage({
                orbit: orbit,
                T: bo.T(),
                run_time: Date.now() - start_time
            });
        }
    }
    orders = [];
}
