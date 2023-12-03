// Collision map for a billiard with a cusp,
// bounded by the curves y=alpha*x^beta, y=-alpha*x^beta and x=1.
// INPUT: point in the table {x: x, y: y, theta: theta}, alpha and beta.
// x and y are the coordinates on the plane, and theta is the angle with
// the positive x axes. The point does not have to be on the boundary.
// OUTPUT: a point in the format above, right after the first collision.
// Collisions with x=1 can be identified by x=1 precisely.

function CUSP(p, alpha, beta) {
    function iter(p0) {
        // x is in [0,1], y is nonnegative
        let x0 = p0.x,
            y0 = p0.y,
            a = p0.a,
            b = p0.b;
        let t1 = 2 * (1 + alpha),
            x1 = x0, 
            y1 = y0;
        // check collision with right wall
        if (a > 0) {
            let t = (1 - x0) / a;
            if (t < t1) {
                t1 = t;
                x1 = 1;
                y1 = y0 + b * t1;
            }
        }
        // check collision with floor
        if (b < 0) {
            let t = - y0 / b;
            if (t < t1) {
                t1 = t;
                x1 = x0 + a * t1;
                y1 = 0;
            }
        }
        // estimate collision with the curve
        let f0 = y0 - alpha * Math.pow(x0, beta),
            fp0 = b - alpha * beta * a * Math.pow(x0, beta - 1);
        if (fp0 > 0) {
            let t = Math.max(0, -f0 / fp0);
            if (t < t1) {
                t1 = t;
                x1 = x0 + a * t1;
                y1 = y0 + b * t1;
            }
        }
        return {
            x: x1,
            y: y1,
            a: a,
            b: b
        }
    }

    function iter2(p0) {
        if (p0.y > 0 || (p0.y == 0 && p0.b >= 0)) {
            return iter(p0);
        } else {
            let p1 = iter({
                x: p0.x, 
                y: -p0.y,
                a: p0.a,
                b: -p0.b
            });
            return {
                x: p1.x,
                y: -p1.y,
                a: p1.a,
                b: -p1.b
            };
        }
    }
    
    // iterate to collision
    let ppp = {
        x: p.x, 
        y: p.y, 
        a: Math.cos(p.theta), 
        b: Math.sin(p.theta)
    };
    let ppp_iter = 0,
        ppp_max_iter = 256;
    do {
        ppp = iter2(ppp);
        ppp_iter ++;
    } while (
        ! (ppp.x == 1 && ppp.a > 0)
        && ppp_iter < ppp_max_iter 
        && alpha * Math.pow(ppp.x, beta) - Math.abs(ppp.y) > 256 * Number.EPSILON
    );

    if (ppp_iter >= ppp_max_iter) {
        console.log("Max iterations!");
    }

    // reflect!
    let phi;
    if (ppp.x == 1 && ppp.a > 0) {
        ppp.a *= -1;
        phi = Math.atan2(-ppp.b, -ppp.a);
    } else { 
        let n = {
            x: alpha * beta * Math.pow(ppp.x, beta - 1),
            y: - Math.sign(ppp.y)
        }
        let n2 = n.x*n.x + n.y*n.y;
        ppp = {
            x: ppp.x,
            y: ppp.y,
            a: ppp.a - 2 * (ppp.a * n.x + ppp.b * n.y) * n.x / n2,
            b: ppp.b - 2 * (ppp.a * n.x + ppp.b * n.y) * n.y / n2
        }
        phi = Math.atan2(ppp.b, ppp.a) - Math.atan2(n.y, n.x);
        while (phi < - Math.PI / 2) {
            phi += Math.PI * 2;
        }
        while (phi > + Math.PI / 2) {
            phi -= Math.PI * 2;
        }
    }

    return {
        x: ppp.x,
        y: ppp.y,
        theta: Math.atan2(ppp.b, ppp.a),
        phi: phi
    };
}

// continuous version: fly until time t.
function CUSPt(p, alpha, beta, t) {
    let pp = p;
    let tt = 0;

    while (true) {
        let ppp = CUSP(pp, alpha, beta);
        let dt = Math.sqrt((pp.x-ppp.x)*(pp.x-ppp.x) + (pp.y-ppp.y)*(pp.y-ppp.y));
        if (tt + dt < t) {
            tt += dt;
            pp = ppp;
        } else {
            return {
                x: pp.x + Math.cos(pp.theta) * (t - tt),
                y: pp.y + Math.sin(pp.theta) * (t - tt),
                theta: pp.theta
            };
        }
    }
}
