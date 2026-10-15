#![allow(dead_code)]

use dashu::float::DBig;
use dashu::integer::IBig;
use std::str::FromStr;

pub const PRECISION: usize = 200;

#[derive(Clone, Debug)]
pub struct TargetConstant {
    pub name: &'static str,
    pub latex: &'static str,
    pub value: DBig,
}

#[allow(dead_code)]
pub struct WeightConstantPool {
    pub w3: Vec<TargetConstant>,
    pub w4: Vec<TargetConstant>,
    pub w5: Vec<TargetConstant>,
}

pub struct NovelConstantPool {
    pub w7: Vec<TargetConstant>,
    pub beta: Vec<TargetConstant>,
    pub w4_cubic: Vec<TargetConstant>,
    pub w5_cubic: Vec<TargetConstant>,
}

/// Cohen-Villegas-Zagier (CVZ) alternating series acceleration algorithm.
/// Computes S = sum_{k=0}^\infty (-1)^k a_k to >205 decimal digits with n = 270 steps.
pub fn cvz_accelerate<F>(n: usize, prec: usize, mut a_k: F) -> DBig
where
    F: FnMut(usize) -> DBig,
{
    let mut t = Vec::with_capacity(n + 1);
    t.push(IBig::ONE);

    for j in 1..=n {
        let num = IBig::from(2 * (n as i64 + j as i64 - 1) * (n as i64 - j as i64 + 1));
        let den = IBig::from(j as i64 * (2 * j as i64 - 1));
        let next_t = &t[j - 1] * num / den;
        t.push(next_t);
    }

    let mut d = IBig::ZERO;
    for val in &t {
        d += val;
    }

    let d_dbig = DBig::from_parts(d.clone(), 0).with_precision(prec).value();
    let mut c = d;
    let mut s = DBig::ZERO.with_precision(prec).value();

    for k in 0..n {
        c -= &t[k];
        let c_dbig = DBig::from_parts(c.clone(), 0).with_precision(prec).value();
        let a = a_k(k);
        let term = &c_dbig * &a;
        if k % 2 == 0 {
            s = &s + &term;
        } else {
            s = &s - &term;
        }
    }

    (&s / &d_dbig).with_precision(prec).value()
}

pub fn compute_pi(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let sixteen = DBig::from_parts(IBig::from(16), 0).with_precision(prec).value();
    let four = DBig::from_parts(IBig::from(4), 0).with_precision(prec).value();

    let arctan = |inv_x: i64, steps: usize| -> DBig {
        let mut sum = DBig::ZERO.with_precision(prec).value();
        let x_d = DBig::from_parts(IBig::from(inv_x), 0).with_precision(prec).value();
        let inv_x_sq = &one / &(&x_d * &x_d);
        let mut cur_power = &one / &x_d;

        for m in 0..steps {
            let denom_val = 2 * m as i64 + 1;
            let denom_d = DBig::from_parts(IBig::from(denom_val), 0).with_precision(prec).value();
            let term = &cur_power / &denom_d;
            if m % 2 == 0 {
                sum = &sum + &term;
            } else {
                sum = &sum - &term;
            }
            cur_power = &cur_power * &inv_x_sq;
        }
        sum
    };

    let s1 = (prec * 170) / 200 + 15;
    let s2 = (prec * 80) / 200 + 15;
    let a1 = arctan(5, s1);
    let a2 = arctan(239, s2);
    (&(&sixteen * &a1) - &(&four * &a2)).with_precision(prec).value()
}

pub fn compute_ln2(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let cvz_n = if prec > 200 { 310 } else { 270 };
    cvz_accelerate(cvz_n, prec, |k| {
        let den = DBig::from_parts(IBig::from(k as i64 + 1), 0).with_precision(prec).value();
        &one / &den
    })
}

pub fn compute_sqrt3(prec: usize) -> DBig {
    let two = DBig::from_parts(IBig::from(2), 0).with_precision(prec).value();
    let three = DBig::from_parts(IBig::from(3), 0).with_precision(prec).value();
    let mut x = DBig::from_str("1.7320508075688772935274463415058723669428").unwrap().with_precision(prec).value();
    for _ in 0..14 {
        let inv = &three / &x;
        x = &(&x + &inv) / &two;
    }
    x.with_precision(prec).value()
}

pub fn compute_sqrt2(prec: usize) -> DBig {
    let two = DBig::from_parts(IBig::from(2), 0).with_precision(prec).value();
    let mut x = DBig::from_str("1.41421356237309504880168872420969807856967187537694807317667973799").unwrap().with_precision(prec).value();
    for _ in 0..14 {
        let inv = &two / &x;
        x = &(&x + &inv) / &two;
    }
    x.with_precision(prec).value()
}

pub fn compute_catalan(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let cvz_n = if prec > 200 { 310 } else { 270 };
    cvz_accelerate(cvz_n, prec, |k| {
        let odd = 2 * k as i64 + 1;
        let den = DBig::from_parts(IBig::from(odd * odd), 0).with_precision(prec).value();
        &one / &den
    })
}

pub fn compute_beta4(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let cvz_n = if prec > 200 { 310 } else { 270 };
    cvz_accelerate(cvz_n, prec, |k| {
        let odd = 2 * k as i64 + 1;
        let odd2 = odd * odd;
        let odd4 = IBig::from(odd2) * IBig::from(odd2);
        let den = DBig::from_parts(odd4, 0).with_precision(prec).value();
        &one / &den
    })
}

pub fn compute_beta6(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let cvz_n = if prec > 200 { 310 } else { 270 };
    cvz_accelerate(cvz_n, prec, |k| {
        let odd = 2 * k as i64 + 1;
        let odd2 = odd * odd;
        let odd4 = IBig::from(odd2) * IBig::from(odd2);
        let odd6 = &odd4 * IBig::from(odd2);
        let den = DBig::from_parts(odd6, 0).with_precision(prec).value();
        &one / &den
    })
}

pub fn compute_zeta3(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let four_thirds = DBig::from_parts(IBig::from(4), 0).with_precision(prec).value()
        / DBig::from_parts(IBig::from(3), 0).with_precision(prec).value();
    let cvz_n = if prec > 200 { 310 } else { 270 };
    let eta3 = cvz_accelerate(cvz_n, prec, |k| {
        let m = k as i64 + 1;
        let den_i = IBig::from(m * m) * IBig::from(m);
        let den = DBig::from_parts(den_i, 0).with_precision(prec).value();
        &one / &den
    });
    (&four_thirds * &eta3).with_precision(prec).value()
}

pub fn compute_zeta5(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let sixteen_fifteenths = DBig::from_parts(IBig::from(16), 0).with_precision(prec).value()
        / DBig::from_parts(IBig::from(15), 0).with_precision(prec).value();
    let eta5 = cvz_accelerate(270, prec, |k| {
        let m = k as i64 + 1;
        let m2 = m * m;
        let den_i = IBig::from(m2 * m2) * IBig::from(m);
        let den = DBig::from_parts(den_i, 0).with_precision(prec).value();
        &one / &den
    });
    (&sixteen_fifteenths * &eta5).with_precision(prec).value()
}

pub fn compute_zeta7(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let sixtyfour_sixtythirds = DBig::from_parts(IBig::from(64), 0).with_precision(prec).value()
        / DBig::from_parts(IBig::from(63), 0).with_precision(prec).value();
    let eta7 = cvz_accelerate(270, prec, |k| {
        let m = k as i64 + 1;
        let m2 = m * m;
        let m4 = IBig::from(m2) * IBig::from(m2);
        let m7 = &m4 * IBig::from(m2 * m);
        let den = DBig::from_parts(m7, 0).with_precision(prec).value();
        &one / &den
    });
    (&sixtyfour_sixtythirds * &eta7).with_precision(prec).value()
}

pub fn compute_li4_half(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let mut sum = DBig::ZERO.with_precision(prec).value();
    let mut two_pow = IBig::ONE;
    for n in 1..=700 {
        two_pow = &two_pow * &IBig::from(2);
        let n_i = n as i64;
        let n4 = IBig::from(n_i * n_i) * IBig::from(n_i * n_i);
        let den_i = &two_pow * &n4;
        let den = DBig::from_parts(den_i, 0).with_precision(prec).value();
        sum = &sum + &(&one / &den);
    }
    sum.with_precision(prec).value()
}

pub fn compute_li5_half(prec: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let mut sum = DBig::ZERO.with_precision(prec).value();
    let mut two_pow = IBig::ONE;
    for n in 1..=700 {
        two_pow = &two_pow * &IBig::from(2);
        let n_i = n as i64;
        let n2 = n_i * n_i;
        let n5 = IBig::from(n2 * n2) * IBig::from(n_i);
        let den_i = &two_pow * &n5;
        let den = DBig::from_parts(den_i, 0).with_precision(prec).value();
        sum = &sum + &(&one / &den);
    }
    sum.with_precision(prec).value()
}

/// Independent Koecher-Leshchiner accelerated series for zeta(5) verification:
/// sum_{n=1}^\infty (-1)^n / binom(2n, n) * [ 9/n^5 - 5*H_n^(2)/n^3 ] = -2 * zeta(5)
pub fn compute_zeta5_koecher(prec: usize, steps: usize) -> DBig {
    let one = DBig::ONE.with_precision(prec).value();
    let two = DBig::from_parts(IBig::from(2), 0).with_precision(prec).value();
    let mut sum = DBig::ZERO.with_precision(prec).value();
    let mut binom = IBig::ONE;
    let mut h2 = DBig::ZERO.with_precision(prec).value();

    for n in 1..=steps {
        binom = &binom * IBig::from(2 * (2 * n as i64 - 1)) / IBig::from(n as i64);

        let n_i = IBig::from(n as i64);
        let n_d = DBig::from_parts(n_i, 0).with_precision(prec).value();
        let inv_n = &one / &n_d;
        let inv_n2 = &inv_n * &inv_n;
        let inv_n3 = &inv_n2 * &inv_n;
        let inv_n5 = &inv_n3 * &inv_n2;
        h2 = &h2 + &inv_n2;

        let binom_d = DBig::from_parts(binom.clone(), 0).with_precision(prec).value();
        let u_n = &one / &binom_d;

        let nine = DBig::from_parts(IBig::from(9), 0).with_precision(prec).value();
        let five = DBig::from_parts(IBig::from(5), 0).with_precision(prec).value();

        let t1 = &nine * &inv_n5;
        let t2 = &(&five * &h2) * &inv_n3;
        let bracket = &t1 - &t2;
        let term = &u_n * &bracket;

        if n % 2 == 1 {
            sum = &sum - &term;
        } else {
            sum = &sum + &term;
        }
    }

    let neg_sum = -sum;
    (&neg_sum / &two).with_precision(prec).value()
}

impl WeightConstantPool {
    pub fn new() -> Self {
        let prec = PRECISION;

        let pi = compute_pi(prec);
        let ln2 = compute_ln2(prec);
        let g = compute_catalan(prec);
        let sqrt3 = compute_sqrt3(prec);
        let zeta3 = compute_zeta3(prec);
        let beta4 = compute_beta4(prec);
        let li4_half = compute_li4_half(prec);
        let zeta5 = compute_zeta5(prec);
        let li5_half = compute_li5_half(prec);

        let pi2 = &pi * &pi;
        let pi3 = &pi2 * &pi;
        let pi4 = &pi2 * &pi2;

        let ln2_sq = &ln2 * &ln2;
        let ln2_cube = &ln2_sq * &ln2;
        let ln2_4 = &ln2_sq * &ln2_sq;
        let ln2_5 = &ln2_4 * &ln2;

        let make_tc = |name: &'static str, latex: &'static str, val: DBig| -> TargetConstant {
            TargetConstant {
                name,
                latex,
                value: val,
            }
        };

        // --- Grup W3 (Ağırlık-3 Sabitleri - 5 Adet) ---
        let pi2_ln2 = &pi2 * &ln2;
        let pi3_over_sqrt3 = (&pi3 / &sqrt3).with_precision(prec).value();
        let pi_g = &pi * &g;

        let w3 = vec![
            make_tc("ζ(3)", "\\zeta(3)", zeta3.clone()),
            make_tc("π²ln(2)", "\\pi^2 \\ln 2", pi2_ln2),
            make_tc("ln³(2)", "\\ln^3 2", ln2_cube.clone()),
            make_tc("π³/√3", "\\frac{\\pi^3}{\\sqrt{3}}", pi3_over_sqrt3),
            make_tc("π·G", "\\pi G", pi_g),
        ];

        // --- Grup W4 (Ağırlık-4 Sabitleri - 6 Adet) ---
        let zeta3_ln2 = &zeta3 * &ln2;
        let pi2_ln2_sq = &pi2 * &ln2_sq;

        let w4 = vec![
            make_tc("π⁴", "\\pi^4", pi4.clone()),
            make_tc("ζ(3)ln(2)", "\\zeta(3) \\ln 2", zeta3_ln2),
            make_tc("π²ln²(2)", "\\pi^2 \\ln^2 2", pi2_ln2_sq),
            make_tc("ln⁴(2)", "\\ln^4 2", ln2_4),
            make_tc("Li₄(1/2)", "\\text{Li}_4(1/2)", li4_half),
            make_tc("β(4)", "\\beta(4)", beta4),
        ];

        // --- Grup W5 (Ağırlık-5 Sabitleri - 7 Adet) ---
        let pi2_zeta3 = &pi2 * &zeta3;
        let pi4_ln2 = &pi4 * &ln2;
        let zeta3_ln2_sq = &zeta3 * &ln2_sq;
        let pi2_ln2_cube = &pi2 * &ln2_cube;

        let w5 = vec![
            make_tc("ζ(5)", "\\zeta(5)", zeta5),
            make_tc("π²ζ(3)", "\\pi^2 \\zeta(3)", pi2_zeta3),
            make_tc("π⁴ln(2)", "\\pi^4 \\ln 2", pi4_ln2),
            make_tc("ζ(3)ln²(2)", "\\zeta(3) \\ln^2 2", zeta3_ln2_sq),
            make_tc("π²ln³(2)", "\\pi^2 \\ln^3 2", pi2_ln2_cube),
            make_tc("ln⁵(2)", "\\ln^5 2", ln2_5),
            make_tc("Li₅(1/2)", "\\text{Li}_5(1/2)", li5_half),
        ];

        Self { w3, w4, w5 }
    }
}

impl NovelConstantPool {
    pub fn new() -> Self {
        let prec = PRECISION;

        let pi = compute_pi(prec);
        let ln2 = compute_ln2(prec);
        let g = compute_catalan(prec);
        let zeta3 = compute_zeta3(prec);
        let zeta5 = compute_zeta5(prec);
        let zeta7 = compute_zeta7(prec);
        let beta4 = compute_beta4(prec);
        let li4_half = compute_li4_half(prec);

        let pi2 = &pi * &pi;
        let pi3 = &pi2 * &pi;
        let pi4 = &pi2 * &pi2;
        let pi6 = &pi4 * &pi2;

        let ln2_sq = &ln2 * &ln2;
        let _ln2_cube = &ln2_sq * &ln2;
        let ln2_4 = &ln2_sq * &ln2_sq;
        let ln2_5 = &ln2_4 * &ln2;

        let g_sq = &g * &g;
        let pi_g = &pi * &g;
        let pi2_g = &pi2 * &g;
        let pi3_g = &pi3 * &g;

        let make_tc = |name: &'static str, latex: &'static str, val: DBig| -> TargetConstant {
            TargetConstant {
                name,
                latex,
                value: val,
            }
        };

        // --- 1. Grup W7 (ζ(7) Avcısı - 6 Adet) ---
        let pi2_zeta5 = &pi2 * &zeta5;
        let pi4_zeta3 = &pi4 * &zeta3;
        let pi6_ln2 = &pi6 * &ln2;
        let zeta3_sq_ln2 = &(&zeta3 * &zeta3) * &ln2;
        let pi2_g_zeta3 = &pi2_g * &zeta3;

        let w7 = vec![
            make_tc("ζ(7)", "\\zeta(7)", zeta7),
            make_tc("π²ζ(5)", "\\pi^2 \\zeta(5)", pi2_zeta5),
            make_tc("π⁴ζ(3)", "\\pi^4 \\zeta(3)", pi4_zeta3),
            make_tc("π⁶ln 2", "\\pi^6 \\ln 2", pi6_ln2),
            make_tc("ζ(3)²ln 2", "\\zeta(3)^2 \\ln 2", zeta3_sq_ln2),
            make_tc("π²Gζ(3)", "\\pi^2 G \\zeta(3)", pi2_g_zeta3),
        ];

        // --- 2. Grup Beta (Tekil Payda Grubu - 5 Adet) ---
        let pi3_ln2 = &pi3 * &ln2;

        let beta = vec![
            make_tc("β(4)", "\\beta(4)", beta4),
            make_tc("π⁴", "\\pi^4", pi4.clone()),
            make_tc("π²G", "\\pi^2 G", pi2_g.clone()),
            make_tc("G²", "G^2", g_sq),
            make_tc("π³ln 2", "\\pi^3 \\ln 2", pi3_ln2),
        ];

        // --- 3. Grup W4_Cubic (Kübik Binom Ağırlık-4 - 6 Adet) ---
        let zeta3_ln2 = &zeta3 * &ln2;
        let pi2_ln2_sq = &pi2 * &ln2_sq;
        let pi_g_ln2 = &pi_g * &ln2;

        let w4_cubic = vec![
            make_tc("π⁴", "\\pi^4", pi4),
            make_tc("ζ(3)ln 2", "\\zeta(3) \\ln 2", zeta3_ln2),
            make_tc("π²ln² 2", "\\pi^2 \\ln^2 2", pi2_ln2_sq),
            make_tc("ln⁴ 2", "\\ln^4 2", ln2_4),
            make_tc("Li₄(1/2)", "\\text{Li}_4(1/2)", li4_half),
            make_tc("πGln 2", "\\pi G \\ln 2", pi_g_ln2),
        ];

        // --- 4. Grup W5_Cubic (Kübik Binom Ağırlık-5 - 5 Adet) ---
        let pi2_zeta3 = &pi2 * &zeta3;
        let pi_g_ln2_sq = &pi_g * &ln2_sq;

        let w5_cubic = vec![
            make_tc("ζ(5)", "\\zeta(5)", zeta5),
            make_tc("π²ζ(3)", "\\pi^2 \\zeta(3)", pi2_zeta3),
            make_tc("πGln² 2", "\\pi G \\ln^2 2", pi_g_ln2_sq),
            make_tc("π³G", "\\pi^3 G", pi3_g),
            make_tc("ln⁵ 2", "\\ln^5 2", ln2_5),
        ];

        Self {
            w7,
            beta,
            w4_cubic,
            w5_cubic,
        }
    }
}

pub struct RamanujanGuilleraPool {
    pub targets: Vec<TargetConstant>,
    pub pi_powers: Vec<TargetConstant>,
    pub radicals: Vec<TargetConstant>,
    pub catalan_zeta: Vec<TargetConstant>,
}

impl RamanujanGuilleraPool {
    pub fn new() -> Self {
        let prec = PRECISION;
        let one = DBig::ONE.with_precision(prec).value();
        let pi = compute_pi(prec);
        let sqrt2 = compute_sqrt2(prec);
        let sqrt3 = compute_sqrt3(prec);
        let g = compute_catalan(prec);
        let zeta3 = compute_zeta3(prec);

        let pi2 = (&pi * &pi).with_precision(prec).value();
        let pi3 = (&pi2 * &pi).with_precision(prec).value();
        let pi4 = (&pi2 * &pi2).with_precision(prec).value();

        let inv_pi = (&one / &pi).with_precision(prec).value();
        let inv_pi2 = (&one / &pi2).with_precision(prec).value();
        let inv_pi3 = (&one / &pi3).with_precision(prec).value();
        let inv_pi4 = (&one / &pi4).with_precision(prec).value();

        let sqrt2_over_pi = (&sqrt2 / &pi).with_precision(prec).value();
        let sqrt3_over_pi = (&sqrt3 / &pi).with_precision(prec).value();

        let pi_sqrt2 = (&pi * &sqrt2).with_precision(prec).value();
        let pi_sqrt3 = (&pi * &sqrt3).with_precision(prec).value();
        let inv_pi_sqrt2 = (&one / &pi_sqrt2).with_precision(prec).value();
        let inv_pi_sqrt3 = (&one / &pi_sqrt3).with_precision(prec).value();

        let sqrt2_over_pi2 = (&sqrt2 / &pi2).with_precision(prec).value();
        let sqrt3_over_pi2 = (&sqrt3 / &pi2).with_precision(prec).value();

        let g_over_pi = (&g / &pi).with_precision(prec).value();
        let zeta3_over_pi3 = (&zeta3 / &pi3).with_precision(prec).value();

        let make_tc = |name: &'static str, latex: &'static str, val: DBig| -> TargetConstant {
            TargetConstant {
                name,
                latex,
                value: val,
            }
        };

        let targets = vec![
            make_tc("1/π", "\\frac{1}{\\pi}", inv_pi.clone()),
            make_tc("1/π²", "\\frac{1}{\\pi^2}", inv_pi2.clone()),
            make_tc("1/π³", "\\frac{1}{\\pi^3}", inv_pi3.clone()),
            make_tc("1/π⁴", "\\frac{1}{\\pi^4}", inv_pi4.clone()),
            make_tc("√2/π", "\\frac{\\sqrt{2}}{\\pi}", sqrt2_over_pi.clone()),
            make_tc("√3/π", "\\frac{\\sqrt{3}}{\\pi}", sqrt3_over_pi.clone()),
            make_tc("1/(π√2)", "\\frac{1}{\\pi\\sqrt{2}}", inv_pi_sqrt2),
            make_tc("1/(π√3)", "\\frac{1}{\\pi\\sqrt{3}}", inv_pi_sqrt3),
            make_tc("√2/π²", "\\frac{\\sqrt{2}}{\\pi^2}", sqrt2_over_pi2.clone()),
            make_tc("√3/π²", "\\frac{\\sqrt{3}}{\\pi^2}", sqrt3_over_pi2.clone()),
            make_tc("G/π", "\\frac{G}{\\pi}", g_over_pi.clone()),
            make_tc("ζ(3)/π³", "\\frac{\\zeta(3)}{\\pi^3}", zeta3_over_pi3.clone()),
        ];

        let pi_powers = vec![
            make_tc("1/π", "\\frac{1}{\\pi}", inv_pi),
            make_tc("1/π²", "\\frac{1}{\\pi^2}", inv_pi2),
            make_tc("1/π³", "\\frac{1}{\\pi^3}", inv_pi3),
            make_tc("1/π⁴", "\\frac{1}{\\pi^4}", inv_pi4),
        ];

        let radicals = vec![
            make_tc("√2/π", "\\frac{\\sqrt{2}}{\\pi}", sqrt2_over_pi),
            make_tc("√3/π", "\\frac{\\sqrt{3}}{\\pi}", sqrt3_over_pi),
            make_tc("√2/π²", "\\frac{\\sqrt{2}}{\\pi^2}", sqrt2_over_pi2),
            make_tc("√3/π²", "\\frac{\\sqrt{3}}{\\pi^2}", sqrt3_over_pi2),
        ];

        let catalan_zeta = vec![
            make_tc("G/π", "\\frac{G}{\\pi}", g_over_pi),
            make_tc("ζ(3)/π³", "\\frac{\\zeta(3)}{\\pi^3}", zeta3_over_pi3),
        ];

        Self {
            targets,
            pi_powers,
            radicals,
            catalan_zeta,
        }
    }
}

pub const FRONTIER_PRECISION: usize = 220;

pub struct FrontierConstantPool {
    pub reference_targets: Vec<TargetConstant>,
    pub w2_targets: Vec<TargetConstant>,
    pub w3_targets: Vec<TargetConstant>,
    pub w4_targets: Vec<TargetConstant>,
    pub catalan_targets: Vec<TargetConstant>,
}

impl FrontierConstantPool {
    pub fn new() -> Self {
        let prec = FRONTIER_PRECISION;
        let one = DBig::ONE.with_precision(prec).value();
        let pi = compute_pi(prec);
        let sqrt2 = compute_sqrt2(prec);
        let sqrt3 = compute_sqrt3(prec);
        let sqrt6 = (&sqrt2 * &sqrt3).with_precision(prec).value();
        let g = compute_catalan(prec);
        let zeta3 = compute_zeta3(prec);

        let pi2 = (&pi * &pi).with_precision(prec).value();
        let pi3 = (&pi2 * &pi).with_precision(prec).value();
        let pi4 = (&pi2 * &pi2).with_precision(prec).value();

        // 1. Modüler Referans & Kontrol (Ağırlık-1): 1/π, √2/π, √3/π, √6/π
        let inv_pi = (&one / &pi).with_precision(prec).value();
        let sqrt2_over_pi = (&sqrt2 / &pi).with_precision(prec).value();
        let sqrt3_over_pi = (&sqrt3 / &pi).with_precision(prec).value();
        let sqrt6_over_pi = (&sqrt6 / &pi).with_precision(prec).value();

        // 2. Ağırlık-2: 1/π², √2/π², √3/π², √6/π²
        let inv_pi2 = (&one / &pi2).with_precision(prec).value();
        let sqrt2_over_pi2 = (&sqrt2 / &pi2).with_precision(prec).value();
        let sqrt3_over_pi2 = (&sqrt3 / &pi2).with_precision(prec).value();
        let sqrt6_over_pi2 = (&sqrt6 / &pi2).with_precision(prec).value();

        // 3. Ağırlık-3: 1/π³, √2/π³, √3/π³, √6/π³, ζ(3)/π³
        let inv_pi3 = (&one / &pi3).with_precision(prec).value();
        let sqrt2_over_pi3 = (&sqrt2 / &pi3).with_precision(prec).value();
        let sqrt3_over_pi3 = (&sqrt3 / &pi3).with_precision(prec).value();
        let sqrt6_over_pi3 = (&sqrt6 / &pi3).with_precision(prec).value();
        let zeta3_over_pi3 = (&zeta3 / &pi3).with_precision(prec).value();

        // 4. Ağırlık-4: 1/π⁴, √2/π⁴, √3/π⁴, √6/π⁴
        let inv_pi4 = (&one / &pi4).with_precision(prec).value();
        let sqrt2_over_pi4 = (&sqrt2 / &pi4).with_precision(prec).value();
        let sqrt3_over_pi4 = (&sqrt3 / &pi4).with_precision(prec).value();
        let sqrt6_over_pi4 = (&sqrt6 / &pi4).with_precision(prec).value();

        // 5. Catalan Hedefleri: G/π, √2G/π, √3G/π, √6G/π, G/π²
        let g_over_pi = (&g / &pi).with_precision(prec).value();
        let sqrt2_g_over_pi = (&(&sqrt2 * &g) / &pi).with_precision(prec).value();
        let sqrt3_g_over_pi = (&(&sqrt3 * &g) / &pi).with_precision(prec).value();
        let sqrt6_g_over_pi = (&(&sqrt6 * &g) / &pi).with_precision(prec).value();
        let g_over_pi2 = (&g / &pi2).with_precision(prec).value();

        let make_tc = |name: &'static str, latex: &'static str, val: DBig| -> TargetConstant {
            TargetConstant {
                name,
                latex,
                value: val,
            }
        };

        let reference_targets = vec![
            make_tc("1/π", "\\frac{1}{\\pi}", inv_pi),
            make_tc("√2/π", "\\frac{\\sqrt{2}}{\\pi}", sqrt2_over_pi),
            make_tc("√3/π", "\\frac{\\sqrt{3}}{\\pi}", sqrt3_over_pi),
            make_tc("√6/π", "\\frac{\\sqrt{6}}{\\pi}", sqrt6_over_pi),
        ];

        let w2_targets = vec![
            make_tc("1/π²", "\\frac{1}{\\pi^2}", inv_pi2),
            make_tc("√2/π²", "\\frac{\\sqrt{2}}{\\pi^2}", sqrt2_over_pi2),
            make_tc("√3/π²", "\\frac{\\sqrt{3}}{\\pi^2}", sqrt3_over_pi2),
            make_tc("√6/π²", "\\frac{\\sqrt{6}}{\\pi^2}", sqrt6_over_pi2),
        ];

        let w3_targets = vec![
            make_tc("1/π³", "\\frac{1}{\\pi^3}", inv_pi3),
            make_tc("√2/π³", "\\frac{\\sqrt{2}}{\\pi^3}", sqrt2_over_pi3),
            make_tc("√3/π³", "\\frac{\\sqrt{3}}{\\pi^3}", sqrt3_over_pi3),
            make_tc("√6/π³", "\\frac{\\sqrt{6}}{\\pi^3}", sqrt6_over_pi3),
            make_tc("ζ(3)/π³", "\\frac{\\zeta(3)}{\\pi^3}", zeta3_over_pi3),
        ];

        let w4_targets = vec![
            make_tc("1/π⁴", "\\frac{1}{\\pi^4}", inv_pi4),
            make_tc("√2/π⁴", "\\frac{\\sqrt{2}}{\\pi^4}", sqrt2_over_pi4),
            make_tc("√3/π⁴", "\\frac{\\sqrt{3}}{\\pi^4}", sqrt3_over_pi4),
            make_tc("√6/π⁴", "\\frac{\\sqrt{6}}{\\pi^4}", sqrt6_over_pi4),
        ];

        let catalan_targets = vec![
            make_tc("G/π", "\\frac{G}{\\pi}", g_over_pi),
            make_tc("√2G/π", "\\frac{\\sqrt{2}G}{\\pi}", sqrt2_g_over_pi),
            make_tc("√3G/π", "\\frac{\\sqrt{3}G}{\\pi}", sqrt3_g_over_pi),
            make_tc("√6G/π", "\\frac{\\sqrt{6}G}{\\pi}", sqrt6_g_over_pi),
            make_tc("G/π²", "\\frac{G}{\\pi^2}", g_over_pi2),
        ];

        Self {
            reference_targets,
            w2_targets,
            w3_targets,
            w4_targets,
            catalan_targets,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zeta5_precision() {
        let prec = PRECISION;
        let zeta5_cvz = compute_zeta5(prec);
        let zeta5_koecher = compute_zeta5_koecher(prec, 250);

        let diff = if zeta5_cvz >= zeta5_koecher {
            &zeta5_cvz - &zeta5_koecher
        } else {
            &zeta5_koecher - &zeta5_cvz
        };

        let thresh = DBig::from_parts(IBig::ONE, -140).with_precision(prec).value();
        assert!(
            diff < thresh,
            "CVZ zeta(5) and Koecher zeta(5) must agree to >140 decimal places! diff: {:?}",
            diff
        );
    }

    #[test]
    fn test_zeta7_precision() {
        let prec = PRECISION;
        let zeta7_cvz = compute_zeta7(prec);
        // Direct sum of 100 terms
        let one = DBig::ONE.with_precision(prec).value();
        let mut direct = DBig::ZERO.with_precision(prec).value();
        for k in 1..=100 {
            let k_i = k as i64;
            let k2 = k_i * k_i;
            let k4 = IBig::from(k2) * IBig::from(k2);
            let k7 = &k4 * IBig::from(k2 * k_i);
            let den = DBig::from_parts(k7, 0).with_precision(prec).value();
            direct = &direct + &(&one / &den);
        }
        let diff = if zeta7_cvz >= direct {
            &zeta7_cvz - &direct
        } else {
            &direct - &zeta7_cvz
        };
        // Tail of sum 1/k^7 is ~ 1/(6 * 100^6) = 1/(6 * 10^12) ~ 1.6e-13
        let thresh = DBig::from_parts(IBig::ONE, -11).with_precision(prec).value();
        assert!(diff < thresh, "Zeta(7) CVZ must agree with direct sum to tail accuracy");
    }

    #[test]
    fn test_beta4_precision() {
        let prec = PRECISION;
        let beta4_cvz = compute_beta4(prec);
        // Direct sum of 100 terms
        let one = DBig::ONE.with_precision(prec).value();
        let mut direct = DBig::ZERO.with_precision(prec).value();
        for k in 0..100 {
            let odd = 2 * k as i64 + 1;
            let odd2 = odd * odd;
            let odd4 = IBig::from(odd2) * IBig::from(odd2);
            let den = DBig::from_parts(odd4, 0).with_precision(prec).value();
            let term = &one / &den;
            if k % 2 == 0 {
                direct = &direct + &term;
            } else {
                direct = &direct - &term;
            }
        }
        let diff = if beta4_cvz >= direct {
            &beta4_cvz - &direct
        } else {
            &direct - &beta4_cvz
        };
        // Tail is ~ 1/(2 * 200^4) ~ 3e-10
        let thresh = DBig::from_parts(IBig::ONE, -8).with_precision(prec).value();
        assert!(diff < thresh, "Beta(4) CVZ must agree with direct sum to tail accuracy");
    }

    #[test]
    fn test_sqrt2_precision() {
        let prec = PRECISION;
        let sqrt2 = compute_sqrt2(prec);
        let two = DBig::from_parts(IBig::from(2), 0).with_precision(prec).value();
        let sq = (&sqrt2 * &sqrt2).with_precision(prec).value();
        let diff = if sq >= two { &sq - &two } else { &two - &sq };
        let thresh = DBig::from_parts(IBig::ONE, -190).with_precision(prec).value();
        assert!(diff < thresh, "sqrt(2)^2 must equal 2 to >190 decimal places");
    }
}
