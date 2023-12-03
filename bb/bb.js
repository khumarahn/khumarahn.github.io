"use strict";

let reportDiv = document.getElementById("reportDiv");
let last_report_timestamp = 0;

let orbit = [];
let orbit_T = 0;
let orbit_run_time = 0;

let delta_t = 0.001;

const worker = new Worker("worker.js");
let worker_responded = 0;
worker.onmessage = (e) => {
    worker_responded = 1;
    orbit = e.data.orbit;
    orbit_T = e.data.T;
    orbit_run_time = e.data.run_time;
    //console.log("Message received from worker: " + JSON.stringify(e.data));
};
worker.postMessage({dt: 0.0});

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

    const obs_R = 0.4 * obs_H;
    const obs_r = 0.2 * obs_H;
    const obs_color = 0xCCCCCC;
    obs_gr.beginFill(obs_color);
    obs_gr.moveTo(0, obs_H);
    obs_gr.lineTo(obs_R, obs_H);
    obs_gr.arcTo(obs_R, obs_H - obs_R, 0, obs_H - obs_R, obs_R); 
    obs_gr.lineTo(0, obs_H);
    obs_gr.endFill();
    obs_gr.beginFill(obs_color);
    obs_gr.moveTo(obs_H, 0);
    obs_gr.lineTo(obs_H, obs_R);
    obs_gr.arcTo(obs_H - obs_R, obs_R, obs_H - obs_R, 0, obs_R); 
    obs_gr.lineTo(obs_H, 0);
    obs_gr.endFill();
    obs_gr.beginFill(obs_color);
    obs_gr.moveTo(0, 0);
    obs_gr.lineTo(0, obs_R);
    obs_gr.arcTo(obs_R, obs_R, obs_R, 0, obs_R);
    obs_gr.lineTo(0, 0);
    obs_gr.endFill();
    obs_gr.beginFill(obs_color);
    obs_gr.moveTo(obs_H, obs_H);
    obs_gr.lineTo(obs_H - obs_R, obs_H);
    obs_gr.arcTo(obs_H - obs_R, obs_H - obs_R, obs_H, obs_H - obs_R, obs_R); 
    obs_gr.lineTo(obs_H, obs_H);
    obs_gr.endFill();
    obs_gr.beginFill(obs_color);
    obs_gr.drawCircle(obs_H / 2, obs_H / 2, obs_r);
    obs_gr.endFill();

    const obs_tx = app.renderer.generateTexture(obs_gr);
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
