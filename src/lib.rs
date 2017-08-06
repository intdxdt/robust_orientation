extern crate robust_sum;
extern crate robust_scale;
extern crate robust_subtract;
extern crate two_product;

use robust_sum::robust_sum as rsum;
use robust_subtract::robust_subtract as rdiff;
use robust_scale::robust_scale as rscale;
use two_product::two_product as tprod;

const EPSILON: f64 = 1.1102230246251565e-16;
const ERRBOUND3: f64 = (3.0 + 16.0 * EPSILON) * EPSILON;
const ERRBOUND4: f64 = (7.0 + 56.0 * EPSILON) * EPSILON;


///orientation in 2d space
/// < 0 if ccw - c is on left of segment(a, b)
/// > 0 if cw - c is on right of segment(a, b)
/// = 0 if a, b, and c are coplanar
pub fn orientation_2d(a: &[f64], b: &[f64], c: &[f64]) -> f64 {
    let l = (a[1] - c[1]) * (b[0] - c[0]);
    let r = (a[0] - c[0]) * (b[1] - c[1]);
    let det = l - r;
    let s;
    if l > 0.0 {
        if r <= 0.0 {
            return det;
        } else {
            s = l + r;
        }
    } else if l < 0.0 {
        if r >= 0.0 {
            return det;
        } else {
            s = -(l + r);
        }
    } else {
        return det;
    }
    let tol = ERRBOUND3 * s;
    if det >= tol || det <= -tol {
        return det;
    }
    orientation3_exact(a, b, c)
}

///orientation in 3d space
/// < 0 if ccw - d is on left of (a, b, c)
/// > 0 if cw  - d is on right of (a, b, c)
/// = 0 if a, b, c and d are coplanar
pub fn orientation_3d(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> f64 {
    let adx = a[0] - d[0];
    let bdx = b[0] - d[0];
    let cdx = c[0] - d[0];
    let ady = a[1] - d[1];
    let bdy = b[1] - d[1];
    let cdy = c[1] - d[1];
    let adz = a[2] - d[2];
    let bdz = b[2] - d[2];
    let cdz = c[2] - d[2];
    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;

    let det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * adz.abs() + (cdxady.abs() + adxcdy.abs()) * bdz.abs() +
        (adxbdy.abs() + bdxady.abs()) * cdz.abs();

    let tol = ERRBOUND4 * permanent;
    if (det > tol) || (-det > tol) {
        return det;
    }
    orientation4_exact(a, b, c, d)
}

///orientation 2d exact
fn orientation3_exact(m0: &[f64], m1: &[f64], m2: &[f64]) -> f64 {
    let p = rsum(
        &rsum(&tprod(m1[1], m2[0]), &tprod(-m2[1], m1[0])),
        &rsum(&tprod(m0[1], m1[0]), &tprod(-m1[1], m0[0])),
    );
    let n = rsum(&tprod(m0[1], m2[0]), &tprod(-m2[1], m0[0]));
    let d = rdiff(&p, &n);
    d[d.len() - 1]
}

///orientation 3d exact
fn orientation4_exact(m0: &[f64], m1: &[f64], m2: &[f64], m3: &[f64]) -> f64 {
    let p = rsum(
        &rsum(
            &rscale(&rsum(&tprod(m2[1], m3[0]), &tprod(-m3[1], m2[0])), m1[2]),
            &rsum(
                &rscale(&rsum(&tprod(m1[1], m3[0]), &tprod(-m3[1], m1[0])), -m2[2]),
                &rscale(&rsum(&tprod(m1[1], m2[0]), &tprod(-m2[1], m1[0])), m3[2]),
            ),
        ),
        &rsum(
            &rscale(&rsum(&tprod(m1[1], m3[0]), &tprod(-m3[1], m1[0])), m0[2]),
            &rsum(
                &rscale(&rsum(&tprod(m0[1], m3[0]), &tprod(-m3[1], m0[0])), -m1[2]),
                &rscale(&rsum(&tprod(m0[1], m1[0]), &tprod(-m1[1], m0[0])), m3[2]),
            ),
        ),
    );

    let n = rsum(
        &rsum(
            &rscale(&rsum(&tprod(m2[1], m3[0]), &tprod(-m3[1], m2[0])), m0[2]),
            &rsum(
                &rscale(&rsum(&tprod(m0[1], m3[0]), &tprod(-m3[1], m0[0])), -m2[2]),
                &rscale(&rsum(&tprod(m0[1], m2[0]), &tprod(-m2[1], m0[0])), m3[2]),
            ),
        ),
        &rsum(
            &rscale(&rsum(&tprod(m1[1], m2[0]), &tprod(-m2[1], m1[0])), m0[2]),
            &rsum(
                &rscale(&rsum(&tprod(m0[1], m2[0]), &tprod(-m2[1], m0[0])), -m1[2]),
                &rscale(&rsum(&tprod(m0[1], m1[0]), &tprod(-m1[1], m0[0])), m2[2]),
            ),
        ),
    );
    let d = rdiff(&p, &n);

    d[d.len() - 1]
}

#[cfg(test)]
mod robust_orientation {
    #[test]
    fn test_ro() {
        assert!(true);
    }
}

