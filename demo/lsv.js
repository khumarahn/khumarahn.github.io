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

// 1st derivative
function LSV_left_p(x, alpha) {
    let gamma = 1. / alpha;
    return 1 + (gamma + 1) * Math.pow(2 * x, gamma);
}
function LSV_right_p(x, alpha) {
        return 2;
}

// 2nd derivative
function LSV_left_pp(x, alpha) {
    let gamma = 1. / alpha;
    return gamma * (gamma + 1) * Math.pow(2, gamma) * Math.pow(x, gamma - 1);
}
function LSV_right_pp(x, alpha) {
    return 0;
}

// 3rd derivative
function LSV_left_ppp(x, alpha) {
    let gamma = 1. / alpha;
    return (gamma - 1) * gamma * (gamma + 1) * Math.pow(2, gamma) * Math.pow(x, gamma - 2);
}
function LSV_left_ppp(x, alpha) {
    return 0;
}

function LSV_left_w(x, alpha) {
    return 1. / LSV_left_p(x, alpha);
}
function LSV_right_w(x, alpha) {
    return 1. / LSV_right_p(x, alpha);
}

function LSV_left_wp(x, alpha) {
    let p = LSV_left_p(x, alpha),
        pp = LSV_left_pp(x, alpha);
    return - pp / (p * p);
}
function LSV_right_wp(x, alpha) {
    return 0;
}

function LSV_left_wpp(x, alpha) {
    let p = LSV_left_p(x, alpha),
        pp = LSV_left_pp(x, alpha),
        ppp = LSV_left_ppp(x, alpha);
    return ( - ppp * p + 2 * pp * pp ) / (p * p * p);
}
function LSV_right_wpp(x, alpha) {
    return 0;
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

// transfer operator wrt Lebesgue (L^n v)(x) with first two derivatives
function LSV_Ln_pp(vvv, x, alpha, n) {
    if (n <= 0) {
        return vvv(x);
    } else {
        let y1 = LSV_left_i(x, alpha),
            y2 = LSV_right_i(x, alpha),
            v1 = LSV_Ln_pp(vvv, y1, alpha, n-1),
            v2 = LSV_Ln_pp(vvv, y2, alpha, n-1),
            w1 = LSV_left_w(y1, alpha),
            w2 = LSV_right_w(y2, alpha),
            wp1 = LSV_left_wp(y1, alpha),
            wp2 = LSV_right_wp(y2, alpha),
            wpp1 = LSV_left_wpp(y1, alpha),
            wpp2 = LSV_right_wpp(y2, alpha);
        return [
            v1[0] * w1 +
            v2[0] * w2,
            //
            v1[1] * w1 * w1 + v1[0] * wp1 * w1 +
            v2[1] * w2 * w2 + v2[0] * wp2 * w2,
            //
            v1[2] * w1 * w1 * w1 + 3 * v1[1] * wp1 * w1 * w1 + v1[0] * (wpp1 * w1 * w1 + wp1 * wp1 * w1) +
            v2[2] * w2 * w2 * w2 + 3 * v2[1] * wp2 * w2 * w2 + v2[0] * (wpp2 * w2 * w2 + wp2 * wp2 * w2)
        ];
    }
}
