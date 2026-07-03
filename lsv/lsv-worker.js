"use strict";

var lsv_cpp;
var orders = [];

self.Module = {
    print: function(text) {
        self.postMessage({
            type: 'stdout',
            text: text
        });
    },
    onRuntimeInitialized: function() {
        lsv_cpp = new Module.LSV();
        runOrders();
    }
};

onmessage = function(e) {
    orders.push(e.data);
    runOrders();
}

function runOrders() {
    if (orders.length > 0) {
        let e = orders.shift();
        if (e.type === 'compute-bounds') {
            let gamma = e.gamma;

            lsv_cpp.set_gamma(gamma);
            lsv_cpp.compute_L();
            lsv_cpp.compute_h_meta();
            lsv_cpp.compute_h_cheb();
            lsv_cpp.compute_F();
            lsv_cpp.compute_derivative_signs_right();
            lsv_cpp.compute_derivative_bounds();

            postMessage({
                type:               'bounds',
                gamma:              lsv_cpp.oracle("gamma"),
                gamma_minus:        lsv_cpp.oracle("gamma-"),
                gamma_plus:         lsv_cpp.oracle("gamma+"),
                alpha:              lsv_cpp.oracle("alpha"),
                alpha_minus:        lsv_cpp.oracle("alpha-"),
                alpha_plus:         lsv_cpp.oracle("alpha+"),
                min_hp_h_prime:     lsv_cpp.oracle("min_hp_h_prime"),
                max_hpp_h_prime:    lsv_cpp.oracle("max_hpp_h_prime")
            });
        } else {
            console.log('worker: unknown order ', e);
        }
    }

    if (orders.length > 0) {
        setTimeout(runOrders, 1000);
    }
}

importScripts('lsv-cpp.js');
