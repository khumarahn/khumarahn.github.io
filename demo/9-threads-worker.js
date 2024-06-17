let all = 0;
let inside = 0;

function calc(iterations) {
    for (let i = 0; i < iterations; i++) {
        let x = Math.random();
        let y = Math.random();
        if (x * x + y * y <= 1) { 
            inside++;
        }
        all++;
    }
}

onmessage = (e) => {
    calc(e.data.iterations);
    
    postMessage({
        pi: 4 * inside / all,
        inside: inside,
        all: all
    });
};
