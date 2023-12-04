"use strict";

function p_var_backbone(path_size, p, path_dist) {
	if (path_size == 0) {
		return -Infinity;
	}
	else if (path_size == 1) {
		return 0.0;
	}

	let run_p_var = [];
	run_p_var.length = path_size;
	run_p_var.fill(0.0);

	let s = path_size - 1;

	let N = 1;
	while (s >> N) {
		N++;
	}

	let ind = [];
	ind.length = s;
	ind.fill(0.0);
	let ind_n = (j, n) => (s >> n) + (j >> n);
	let ind_k = (j, n) => Math.min(((j >> n) << n) + (1 << (n-1)), s);

	let max_p_var = 0.0;

	let point_links = [];
	point_links.length = path_size;
	point_links.fill(0);

	for (let j = 0; j < path_size; j++) {
		for (let n = 1; n <= N; n++) {
			if (!(j >> n == s >> n && (s >> (n-1)) % 2 == 0)) {
				ind[ind_n(j, n)] = Math.max(ind[ind_n(j, n)], path_dist(ind_k(j, n), j));
			}
		}
		if (j == 0) {
			continue;
		}

		let m = j - 1;
		let delta = 0;
		let delta_m = j;
		for (let n=0;;) {
			while (n > 0 && m >> n == s >> n && (s >> (n-1)) % 2 == 0) {
				n--;
			}

			let skip = false;
			if (n > 0) {
				let id = ind[ind_n(m, n)] + path_dist(ind_k(m, n), j);
				if (delta >= id) {
					skip = true;
				}
				else if (m < delta_m) {
					delta = Math.pow(max_p_var - run_p_var[m], 1. / p);
					delta_m = m;
					if (delta >= id) {
						skip = true;
					}
				}
			}

			if (skip) {
				let k = (m >> n) << n;
				if (k > 0) {
					m = k - 1;
					while (n < N && (k>>n) % 2 == 0) {
						n++;
					}
				}
				else {
					break;
				}
			}
			else {
				if (n > 1) {
					n--;
				}
				else {
					let d = path_dist(m, j);
					if (d >= delta) {
						let new_p_var = run_p_var[m] + Math.pow(d, p);
						if (new_p_var >= max_p_var) {
							max_p_var = new_p_var;
							point_links[j] = m;
						}
					}

					if (m > 0) {
						while (n < N  &&  (m>>n) % 2 == 0) {
							n++;
						}
						m--;
					}
					else {
						break;
					}
				}
			}
		}

		run_p_var[j] = max_p_var;
	}

	let points = [];
	for (let q = s; ; q = point_links[q]) {
		points.push(q);
		if (q == 0) {
			break;
		}
	}
	points.reverse();

	return {
		"p_var": run_p_var[run_p_var.length - 1],
		"points": points,
		"run_p_var": run_p_var
	};
}
