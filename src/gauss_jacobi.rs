//! Definitions of Gauss-Jacobi quadrature rules.
//!
//! Adapted from the C++ code by Chris Richardson
//! <https://github.com/FEniCS/basix/blob/main/cpp/basix/quadrature.cpp>

use crate::types::NumericalQuadratureDefinition;
use itertools::izip;

/// Evaluate the nth Jacobi polynomial and derivatives with weight
/// parameters (a, 0) at points x
fn compute_deriv(a: f64, n: usize, nderiv: usize, x: &[f64]) -> Vec<f64> {
    let mut j_all = vec![0.0; (nderiv + 1) * (n + 1) * x.len()];

    for i in 0..=nderiv {
        if i == 0 {
            for j in 0..x.len() {
                j_all[i * (n + 1) * x.len() + j] = 1.0;
            }
        }

        if n > 0 {
            if i == 0 {
                for (j, x_j) in x.iter().enumerate() {
                    j_all[i * (n + 1) * x.len() + x.len() + j] = (*x_j * (a + 2.0) + a) / 2.0;
                }
            } else if i == 1 {
                for j in 0..x.len() {
                    j_all[i * (n + 1) * x.len() + x.len() + j] = a / 2.0 + 1.0;
                }
            }
        }
        for j in 2..=n {
            let j_t = j as f64;
            let a1 = 2.0 * j_t * (j_t + a) * (2.0 * j_t + a - 2.0);
            let a2 = (2.0 * j_t + a - 1.0) * (a * a) / a1;
            let a3 = (2.0 * j_t + a - 1.0) * (2.0 * j_t + a) / (2.0 * j_t * (j_t + a));
            let a4 = 2.0 * (j_t + a - 1.0) * (j_t - 1.0) * (2.0 * j_t + a) / a1;
            for (k, x_k) in x.iter().enumerate() {
                j_all[i * (n + 1) * x.len() + j * x.len() + k] =
                    j_all[i * (n + 1) * x.len() + (j - 1) * x.len() + k] * (*x_k * a3 + a2)
                        - j_all[i * (n + 1) * x.len() + (j - 2) * x.len() + k] * a4;
            }
            if i > 0 {
                for (k, _) in x.iter().enumerate() {
                    let add =
                        i as f64 * a3 * j_all[(i - 1) * (n + 1) * x.len() + (j - 1) * x.len() + k];
                    j_all[i * (n + 1) * x.len() + j * x.len() + k] += add;
                }
            }
        }
    }

    let mut j = vec![0.0; (nderiv + 1) * x.len()];
    for i0 in 0..=nderiv {
        for i1 in 0..x.len() {
            j[i0 * x.len() + i1] = j_all[i0 * (n + 1) * x.len() + n * x.len() + i1];
        }
    }

    j
}

/// Computes the m roots of \f$P_{m}^{a,0}\f$ on [-1,1] by Newton's
/// method. The initial guesses are the Chebyshev points.  Algorithm
/// implemented from the pseudocode given by Karniadakis and Sherwin.
fn compute_points(a: f64, m: usize) -> Vec<f64> {
    let eps = 1.0e-8;
    let max_iter = 100;

    let mut x = vec![0.0; m];

    for k in 0..m {
        // Initial guess
        x[k] = -f64::cos((2 * k + 1) as f64 * std::f64::consts::PI / (2.0 * m as f64));
        if k > 0 {
            x[k] = (x[k] + x[k - 1]) / 2.0;
        }

        for _ in 0..max_iter {
            let s = x.iter().take(k).map(|i| 1.0 / (x[k] - *i)).sum::<f64>();
            let f = compute_deriv(a, m, 1, &x[k..k + 1]);
            let delta: f64 = f[0] / (f[1] - f[0] * s);
            x[k] -= delta;
            if delta.abs() < eps {
                break;
            }
        }
    }
    x
}

/// Note: computes on [-1, 1]
fn compute_rule(a: f64, m: usize) -> (Vec<f64>, Vec<f64>) {
    let pts = compute_points(a, m);
    let j_d = compute_deriv(a, m, 1, &pts);
    let a1 = f64::powf(2.0, a + 1.0);
    let wts = pts
        .iter()
        .enumerate()
        .map(|(i, x)| a1 / (1.0 - x.powi(2)) / j_d[pts.len() + i].powi(2))
        .collect::<Vec<_>>();

    (pts, wts)
}

/// Get the points and weights for a Gauss-Jacobi quadrature rule on an interval
pub fn gauss_jacobi_interval(m: usize) -> NumericalQuadratureDefinition {
    let (mut pts, mut wts) = compute_rule(0.0, m);

    for p in pts.iter_mut() {
        *p = 0.5 * (*p + 1.0);
    }
    for w in wts.iter_mut() {
        *w *= 0.5;
    }
    NumericalQuadratureDefinition {
        dim: 1,
        order: m,
        npoints: wts.len(),
        weights: wts,
        points: pts,
    }
}

/// Get the points and weights for a collapsed Gauss-Jacobi quadrature rule on a triangle
pub fn gauss_jacobi_triangle(m: usize) -> NumericalQuadratureDefinition {
    let (ptx, wtx) = compute_rule(0.0, m);
    let (pty, wty) = compute_rule(1.0, m);

    let mut pts = vec![0.0; m.pow(2) * 2];
    let mut wts = vec![0.0; m.pow(2)];

    for (i, (px, wx)) in izip!(&ptx, &wtx).enumerate() {
        for (j, (py, wy)) in izip!(&pty, &wty).enumerate() {
            let index = i * wty.len() + j;
            pts[2 * index] = 0.25 * (1.0 + *px) * (1.0 - *py);
            pts[2 * index + 1] = 0.5 * (1.0 + *py);
            wts[index] = *wx * *wy * 0.125;
        }
    }
    NumericalQuadratureDefinition {
        dim: 2,
        order: m,
        npoints: wts.len(),
        weights: wts,
        points: pts,
    }
}

/// Get the points and weights for a collapsed Gauss-Jacobi quadrature rule on a tetrahedron
pub fn gauss_jacobi_tetrahedron(m: usize) -> NumericalQuadratureDefinition {
    let (ptx, wtx) = compute_rule(0.0, m);
    let (pty, wty) = compute_rule(1.0, m);
    let (ptz, wtz) = compute_rule(2.0, m);

    let mut pts = vec![0.0; m.pow(3) * 3];
    let mut wts = vec![0.0; m.pow(3)];

    for (i, (px, wx)) in izip!(&ptx, &wtx).enumerate() {
        for (j, (py, wy)) in izip!(&pty, &wty).enumerate() {
            for (k, (pz, wz)) in izip!(&ptz, &wtz).enumerate() {
                let index = i * wty.len() * wtz.len() + j * wtz.len() + k;
                pts[3 * index] = 0.125 * (1.0 + *px) * (1.0 - *py) * (1.0 - *pz);
                pts[3 * index + 1] = 0.25 * (1.0 + *py) * (1.0 - *pz);
                pts[3 * index + 2] = 0.5 * (1.0 + *pz);
                wts[index] = *wx * *wy * *wz * 0.015625;
            }
        }
    }
    NumericalQuadratureDefinition {
        dim: 3,
        order: m,
        npoints: wts.len(),
        weights: wts,
        points: pts,
    }
}

/// Get the points and weights for a collapsed Gauss-Jacobi quadrature rule on a quadrilateral
pub fn gauss_jacobi_quadrilateral(m: usize) -> NumericalQuadratureDefinition {
    let rule = gauss_jacobi_interval(m);
    let mut pts = vec![0.0; m.pow(2) * 2];
    let mut wts = vec![0.0; m.pow(2)];

    for (i, (pi, wi)) in izip!(&rule.points, &rule.weights).enumerate() {
        for (j, (pj, wj)) in izip!(&rule.points, &rule.weights).enumerate() {
            let index = i * m + j;
            pts[2 * index] = *pi;
            pts[2 * index + 1] = *pj;
            wts[index] = *wi * *wj;
        }
    }
    NumericalQuadratureDefinition {
        dim: 2,
        order: m,
        npoints: wts.len(),
        weights: wts,
        points: pts,
    }
}

/// Get the points and weights for a collapsed Gauss-Jacobi quadrature rule on a hexahedron
pub fn gauss_jacobi_hexahedron(m: usize) -> NumericalQuadratureDefinition {
    let rule = gauss_jacobi_interval(m);
    let mut pts = vec![0.0; m.pow(3) * 3];
    let mut wts = vec![0.0; m.pow(3)];

    for (i, (pi, wi)) in izip!(&rule.points, &rule.weights).enumerate() {
        for (j, (pj, wj)) in izip!(&rule.points, &rule.weights).enumerate() {
            for (k, (pk, wk)) in izip!(&rule.points, &rule.weights).enumerate() {
                let index = i * m.pow(2) + j * m + k;
                pts[3 * index] = *pi;
                pts[3 * index + 1] = *pj;
                pts[3 * index + 2] = *pk;
                wts[index] = *wi * *wj * *wk;
            }
        }
    }

    NumericalQuadratureDefinition {
        dim: 3,
        order: m,
        npoints: wts.len(),
        weights: wts,
        points: pts,
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::*;

    #[test]
    fn test_interval_3() {
        let rule = gauss_jacobi_interval(3);

        for (p, q) in izip!(rule.points, [0.1127016653792583, 0.5, 0.8872983346207417]) {
            assert_relative_eq!(p, q);
        }
        for (w, v) in izip!(
            rule.weights,
            [0.2777777777777777, 0.4444444444444444, 0.2777777777777777]
        ) {
            assert_relative_eq!(w, v);
        }
    }
}
