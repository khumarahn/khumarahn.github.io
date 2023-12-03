"use strict";

var Module = {
    'onRuntimeInitialized': init
};   
importScripts('32.js');

var bo; // billiard orbit
let orbit = [{x: 0, y: 0}];
let orbit_T = 0;
let run_time = 0;
let max_orbit = 128*1024;
let init_complete = 0;

function init() {
    bo = new Module.orbit_c();
    init_complete = 1;
}

onmessage = function(e) {
    postMessage({orbit: orbit, T: orbit_T, run_time: run_time});
    if (init_complete) {
        reOrbit(e.data.dt);
    }
}

function reOrbit(dt) {
    let t = Date.now();
    bo.step(dt);
    orbit.length = 0;
    // <= size() because the extra point is the current position
    for (let i=0; i<=bo.size(); i++) {
        orbit.push({x: bo.xy(i).x, y: bo.xy(i).y});
    }
    orbit_T = bo.T();
    run_time = Date.now() - t;
}
