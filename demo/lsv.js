"use strict";

// we work with LSV with alpha, where alpha is the stability index, the left branch is
// x -> x (1 + (2x)^(1 / alpha)

// left branch
function LSV_left(x, alpha) {
    let gamma = 1. / alpha;
    return x * (1 + Math.pow(2 * x, gamma));
}

// right branch
function LSV_right(x, alpha) {
    return 2 * x - 1;
}

// whole map
function LSV(x, alpha) {
    return (x < 0.5) ? LSV_left(x, alpha) : LSV_right(x, alpha);
}

// derivative of left branch
function LSV_left_p(x, alpha) {
    let gamma = 1. / alpha;
    return 1 + (gamma + 1) * Math.pow(2 * x, gamma);
}

// derivative of right branch
function LSV_right_p(x, alpha) {
    return 2;
}

// inverse of left branch
function LSV_left_i(x, alpha) {
    // convex, so Newton should work very well
    let y = x;
    let f,p,z;
    for (;;) {
        f = LSV_left(y, alpha) - x;
        p = LSV_left_p(y, alpha);
        z = y - f / p;
        if (z >= y)
            break;
        y = z;
    }
    return y;
}

// inverse of right branch
function LSV_right_i(x, alpha) {
    return (x + 1) / 2;
}

// transfer operator wrt Lebesgue (L^n v)(x)
function LSV_Ln(v, x, alpha, n) {
    if (n <= 0) {
        return v(x);
    } else {
        let y1 = LSV_left_i(x, alpha),
            p1 = LSV_left_p(y1, alpha),
            y2 = LSV_right_i(x, alpha),
            p2 = LSV_right_p(y2, alpha);
        return LSV_Ln(v, y1, alpha, n-1) / p1 + LSV_Ln(v, y2, alpha, n-1) / p2;
    }
}

