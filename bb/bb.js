"use strict";

let reportDiv = document.getElementById("reportDiv");
let last_report_timestamp = 0;

let orbit = [];
let orbit_T = 0;
let orbit_run_time = 0;

let delta_t = 0.001;


function getJsonFromUrl(url) {
    if(!url) url = location.search;
    var query = url.substr(1);
    var result = {};
    query.split("&").forEach(function(part) {
        var item = part.split("=");
        result[item[0]] = decodeURIComponent(item[1]);
    });
    return result;
}

let table = Math.floor(Math.random() * 2);
{
    let opts = getJsonFromUrl();
    if ("table" in opts) {
        if (opts.table == 1) {
            table = 1;
        } else {
            table = 0;
        }
    }
}

const worker = new Worker("worker.js");
let worker_responded = 0;
worker.onmessage = (e) => {
    worker_responded = 1;
    orbit = e.data.orbit;
    orbit_T = e.data.T;
    orbit_run_time = e.data.run_time;
    //console.log("Message received from worker: " + JSON.stringify(e.data));
};
worker.postMessage({table: table, dt: 0.0});

const app = new PIXI.Application({
    background: '#FFFFFF',
    width: 120,
    height: 120,
    antialias: true
});

let obs;
const obs_H = 1024;

function obs_update() {
    const obs_gr = new PIXI.Graphics();

    const obs_color = 0xCCCCCC;
    
    if (table == 0) {
        const obs_R = 0.4 * obs_H;
        const obs_r = 0.2 * obs_H;
        obs_gr.beginFill(obs_color);
        obs_gr.drawCircle(obs_H / 2, obs_H / 2, obs_r);
        obs_gr.endFill();
        for (let x=0; x <= obs_H; x+= obs_H) {
            for (let y=0; y <= obs_H; y+= obs_H) {
                obs_gr.beginFill(obs_color);
                obs_gr.drawCircle(x, y, obs_R);
                obs_gr.endFill();
            }
        }
    } else {
        // table == 1
        { // long arcs
            const ax =  1.25;
            const ay =  3./16;
            const bx = -1.25;
            const by = 13./16;

            const mx = (ax + bx) * 0.5;
            const my = (ay + by) * 0.5;
            const dab = Math.sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by));
            const nx = (by - ay) / dab;
            const ny = (ax - bx) / dab;

            const R = Math.sqrt(16*16 + dab*dab/4.);

            const q = dab * dab / 4 / Math.sqrt(R * R - dab * dab / 4);

            const xc1 = mx + nx * q;
            const yc1 = my + ny * q;

            const xc2 = mx - nx * q;
            const yc2 = my - ny * q;

            for (let i = -1; i <= 2; i++) {
                obs_gr.beginFill(obs_color);
                obs_gr.moveTo((ax + i) * obs_H, (1 - ay) * obs_H);
                obs_gr.arcTo((xc1 + i) * obs_H, (1 - yc1) * obs_H, (bx + i) * obs_H, (1 - by) * obs_H, R * obs_H);
                obs_gr.arcTo((xc2 + i) * obs_H, (1 - yc2) * obs_H, (ax + i) * obs_H, (1 - ay) * obs_H, R * obs_H);
                obs_gr.endFill();
            }
        }
        { // obstacles: short arcs
            const R = Math.sqrt(0.1875 * 0.1875 + 1 * 1);
            const phi = Math.acos(1. / R);

            const ax = 0.5;
            const ay = - R * Math.sin(phi);

            const bx = ax;
            const by = - ay;

            const xc1 = 0.5 + ay * ay;
            const yc1 = 0;
            const xc2 = 0.5 - ay * ay;
            const yc2 = 0;
            
            for (let j = 0; j <= 1; j++) {
                obs_gr.beginFill(obs_color);
                obs_gr.moveTo(ax * obs_H, (1 - ay - j) * obs_H);
                obs_gr.arcTo(xc1 * obs_H, (1 - yc1 - j) * obs_H, bx * obs_H, (1 - by - j) * obs_H, R * obs_H);
                obs_gr.arcTo(xc2 * obs_H, (1 - yc2 - j) * obs_H, ax * obs_H, (1 - ay - j) * obs_H, R * obs_H);
                obs_gr.endFill();
            }
        }
    }

    const obs_tx = app.renderer.generateTexture(obs_gr, {region: new PIXI.Rectangle(0, 0, obs_H, obs_H)});
    obs = new PIXI.TilingSprite(obs_tx);
    app.stage.addChild(obs);
};
obs_update();

const bm = new PIXI.Graphics();
app.stage.addChild(bm);

function resize() {
    app.renderer.resize(
        Math.max(120, window.innerWidth),
        Math.max(120, window.innerHeight)
    );
    
    if (!obs.destroyed) {
        obs.width = app.screen.width;
        obs.height = app.screen.height;
        obs.tilePosition.x = app.screen.width / 2;
        obs.tilePosition.y = app.screen.height / 2;
    }
}
resize();

window.addEventListener('resize', resize);

document.body.appendChild(app.view);

window.addEventListener('click', function (e) {
    window.location.reload();
});


let scale_factor = 1.0;
let last_tick = Date.now();

app.ticker.add((delta) =>
{
    if (worker_responded) {
    
        let now_tick = Date.now();
        let delta_tick = now_tick - last_tick;
        last_tick = now_tick;

        if (orbit_run_time < 10 && delta_t < 0.2 * orbit_T * delta_tick) {
            delta_t *= 1.001;
        } else {
            delta_t *= 0.999;
        }

        worker_responded = 0;
        worker.postMessage({dt: delta_t});

        let x_max = 0;
        let y_max = 0;
        for (let i=0; i<orbit.length; i++) {
            x_max = Math.max(x_max, Math.abs(orbit[i].x));
            y_max = Math.max(y_max, Math.abs(orbit[i].y));
        }

        let scale = Math.min(
            Math.min(app.screen.width, app.screen.height) * scale_factor * Math.pow(orbit_T, -0.5),
            app.screen.width / 2.05 / x_max,
            app.screen.height / 2.05 / y_max
        );
        if ( x_max * scale > 0.85 * app.screen.width / 2 || y_max * scale > 0.85 * app.screen.height / 2 ) {
            scale_factor /= 1.002;
        } else if ( x_max * scale < 0.7 * app.screen.width / 2 && y_max * scale < 0.7 * app.screen.height / 2 ) {
            scale_factor *= 1.0005;
        }

        bm.clear();
        bm.lineStyle({width: 2, join: PIXI.LINE_JOIN.BEVEL, color: 0xFF0000, alpha: 1});


        if (!obs.destroyed) {
            obs.tileScale.x = scale / obs_H;
            obs.tileScale.y = scale / obs_H;
            obs.alpha = Math.min(1.0, Math.pow((scale - 5) / 22, 1.1));
            if (scale < 5.) {
                obs.destroy();
                app.renderer.background.color = 0xFFFFFF;
            }
        }

        for (let i=0; i<orbit.length; i++) {
            let x = app.screen.width  / 2 + orbit[i].x * scale;
            let y = app.screen.height / 2 - orbit[i].y * scale;
            if (i==0) {
                bm.moveTo(x, y);
            }
            bm.lineTo(x, y);
        }
    } // worker responded
});
