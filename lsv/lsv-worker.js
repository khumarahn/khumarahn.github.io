var orders = [];

onmessage = function(e) {
    orders.push(e.data);
}

function runOrders() {
    if (orders.length > 0) {

        e = orders.shift();

        gamma = e.gamma;

        lsv_cpp.set_gamma(gamma);
        lsv_cpp.compute_L();
        lsv_cpp.compute_h_meta();
        lsv_cpp.compute_h_cheb();
        lsv_cpp.compute_F();
        lsv_cpp.compute_derivative_signs_right();
        lsv_cpp.compute_derivative_bounds();

        postMessage({
            gamma:              lsv_cpp.double_gamma(),
            alpha:              lsv_cpp.double_alpha(),
            min_hp_h_prime:     lsv_cpp.double_min_hp_h_prime(),
            max_hpp_h_prime:    lsv_cpp.double_max_hpp_h_prime()
        });
    }

    setTimeout(runOrders, 100);
}

var lsv_cpp;

var Module = {
    'onRuntimeInitialized': function() {
        lsv_cpp = new Module.LSV();

        runOrders();
    }

};

importScripts('lsv-cpp.js');
