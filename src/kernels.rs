#![allow(dead_code)]

use dashu::float::DBig;
use dashu::integer::IBig;
use dashu::base::BitTest;
use crate::constants::PRECISION;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum KernelType {
    CentralBinomial,        // (2n choose n)
    SquaredCentralBinomial, // (2n choose n)^2
    CubicBinomial,          // (3n choose n)
    QuarticBinomial,        // (4n choose 2n)
}

#[derive(Clone, Debug)]
pub struct SeriesKernel {
    pub a: i64,
    pub kernel_type: KernelType,
}

impl SeriesKernel {
    pub fn new(a: i64, kernel_type: KernelType) -> Self {
        Self { a, kernel_type }
    }
}

#[derive(Clone, Debug)]
pub struct KernelWeightBases {
    pub kernel: SeriesKernel,
    pub w3: [DBig; 5],
    pub w4: [DBig; 6],
    pub w5: [DBig; 8],
}

pub struct HarmonicTables {
    pub h: Vec<DBig>,   // H_m = sum_{j=1}^m 1/j
    pub h2: Vec<DBig>,  // H_m^(2) = sum_{j=1}^m 1/j^2
    pub h3: Vec<DBig>,  // H_m^(3) = sum_{j=1}^m 1/j^3
    pub h4: Vec<DBig>,  // H_m^(4) = sum_{j=1}^m 1/j^4
    pub h5: Vec<DBig>,  // H_m^(5) = sum_{j=1}^m 1/j^5
    pub o: Vec<DBig>,   // O_m = sum_{j=1}^m 1/(2j-1)
    pub o2: Vec<DBig>,  // O_m^(2) = sum_{j=1}^m 1/(2j-1)^2
    pub n_dbig: Vec<DBig>,
}

impl HarmonicTables {
    pub fn build(max_m: usize, prec: usize) -> Self {
        let mut h = Vec::with_capacity(max_m + 1);
        let mut h2 = Vec::with_capacity(max_m + 1);
        let mut h3 = Vec::with_capacity(max_m + 1);
        let mut h4 = Vec::with_capacity(max_m + 1);
        let mut h5 = Vec::with_capacity(max_m + 1);
        let mut o = Vec::with_capacity(max_m + 1);
        let mut o2 = Vec::with_capacity(max_m + 1);
        let mut n_dbig = Vec::with_capacity(max_m + 1);

        let one_d = DBig::ONE.with_precision(prec).value();
        let mut cur_h = DBig::ZERO.with_precision(prec).value();
        let mut cur_h2 = DBig::ZERO.with_precision(prec).value();
        let mut cur_h3 = DBig::ZERO.with_precision(prec).value();
        let mut cur_h4 = DBig::ZERO.with_precision(prec).value();
        let mut cur_h5 = DBig::ZERO.with_precision(prec).value();
        let mut cur_o = DBig::ZERO.with_precision(prec).value();
        let mut cur_o2 = DBig::ZERO.with_precision(prec).value();

        h.push(cur_h.clone());
        h2.push(cur_h2.clone());
        h3.push(cur_h3.clone());
        h4.push(cur_h4.clone());
        h5.push(cur_h5.clone());
        o.push(cur_o.clone());
        o2.push(cur_o2.clone());
        n_dbig.push(DBig::ZERO.with_precision(prec).value());

        for m in 1..=max_m {
            let m_i = IBig::from(m as i64);
            let m_d = DBig::from_parts(m_i, 0).with_precision(prec).value();
            let inv_m = &one_d / &m_d;
            let inv_m2 = &inv_m * &inv_m;
            let inv_m3 = &inv_m2 * &inv_m;
            let inv_m4 = &inv_m2 * &inv_m2;
            let inv_m5 = &inv_m4 * &inv_m;

            cur_h = &cur_h + &inv_m;
            cur_h2 = &cur_h2 + &inv_m2;
            cur_h3 = &cur_h3 + &inv_m3;
            cur_h4 = &cur_h4 + &inv_m4;
            cur_h5 = &cur_h5 + &inv_m5;

            // Odd harmonic terms (2m - 1)
            let odd_val = 2 * m as i64 - 1;
            let odd_i = IBig::from(odd_val);
            let odd_d = DBig::from_parts(odd_i, 0).with_precision(prec).value();
            let inv_odd = &one_d / &odd_d;
            let inv_odd2 = &inv_odd * &inv_odd;

            cur_o = &cur_o + &inv_odd;
            cur_o2 = &cur_o2 + &inv_odd2;

            h.push(cur_h.clone());
            h2.push(cur_h2.clone());
            h3.push(cur_h3.clone());
            h4.push(cur_h4.clone());
            h5.push(cur_h5.clone());
            o.push(cur_o.clone());
            o2.push(cur_o2.clone());
            n_dbig.push(m_d);
        }

        Self {
            h,
            h2,
            h3,
            h4,
            h5,
            o,
            o2,
            n_dbig,
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════════
//  ALAN A: AĞIRLIK-7 (ζ(7)) HOMOJEN BAZLARI (10 Baz Seri)
//  C(n) = a^n \binom{2n}{n}, a \in {-2, -1, 1, 2}
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Debug)]
pub struct W7KernelBases {
    pub a: i64,
    pub bases: [DBig; 10],
}

pub fn compute_w7_bases(a: i64, steps: usize, tables: &HarmonicTables) -> W7KernelBases {
    let prec = PRECISION;
    let mut bases = [
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
    ];

    let one_d = DBig::ONE.with_precision(prec).value();
    let mut a_pow = IBig::ONE;
    let a_ibig = IBig::from(a);
    let mut binom_2n = IBig::ONE;

    for n in 1..=steps {
        a_pow = &a_pow * &a_ibig;
        binom_2n = &binom_2n * IBig::from(2 * (2 * n as i64 - 1)) / IBig::from(n as i64);

        let den_ibig = &a_pow * &binom_2n;
        if den_ibig.bit_len() > 720 {
            break;
        }

        let den_d = DBig::from_parts(den_ibig, 0).with_precision(prec).value();
        let u_n = &one_d / &den_d; // 1 / C(n)

        let h_n = &tables.h[n];
        let h_2n = &tables.h[2 * n];
        let h2_n = &tables.h2[n];
        let h2_2n = &tables.h2[2 * n];
        let h3_n = &tables.h3[n];
        let h4_n = &tables.h4[n];

        let h_h2 = h_n * h2_n;
        let h2_sq = h2_n * h2_n;
        let h_h3 = h_n * h3_n;

        let n_d = &tables.n_dbig[n];
        let inv_n = &one_d / n_d;
        let inv_n2 = &inv_n * &inv_n;
        let inv_n3 = &inv_n2 * &inv_n;
        let inv_n4 = &inv_n2 * &inv_n2;
        let inv_n5 = &inv_n4 * &inv_n;
        let inv_n6 = &inv_n5 * &inv_n;
        let inv_n7 = &inv_n6 * &inv_n;

        // 0: 1 / n^7
        bases[0] = &bases[0] + &(&u_n * &inv_n7);
        // 1: H_n / n^6
        bases[1] = &bases[1] + &(&(&u_n * h_n) * &inv_n6);
        // 2: H_{2n} / n^6
        bases[2] = &bases[2] + &(&(&u_n * h_2n) * &inv_n6);
        // 3: H_n^(2) / n^5
        bases[3] = &bases[3] + &(&(&u_n * h2_n) * &inv_n5);
        // 4: H_{2n}^(2) / n^5
        bases[4] = &bases[4] + &(&(&u_n * h2_2n) * &inv_n5);
        // 5: H_n^(3) / n^4
        bases[5] = &bases[5] + &(&(&u_n * h3_n) * &inv_n4);
        // 6: (H_n * H_n^(2)) / n^4
        bases[6] = &bases[6] + &(&(&u_n * &h_h2) * &inv_n4);
        // 7: H_n^(4) / n^3
        bases[7] = &bases[7] + &(&(&u_n * h4_n) * &inv_n3);
        // 8: (H_n^(2))^2 / n^3
        bases[8] = &bases[8] + &(&(&u_n * &h2_sq) * &inv_n3);
        // 9: (H_n * H_n^(3)) / n^3
        bases[9] = &bases[9] + &(&(&u_n * &h_h3) * &inv_n3);
    }

    W7KernelBases { a, bases }
}

pub fn format_w7_lhs_latex(a: i64, c: &[i64]) -> String {
    let denom = match a {
        1 => "\\binom{2n}{n}".to_string(),
        -1 => "(-1)^n \\binom{2n}{n}".to_string(),
        val if val < 0 => format!("({val})^n \\binom{{2n}}{{n}}"),
        val => format!("{val}^n \\binom{{2n}}{{n}}"),
    };

    let has_harm = c[1..10].iter().any(|&v| v != 0);
    if !has_harm {
        return format!("\\sum_{{n=1}}^\\infty \\frac{{{}}}{{n^7 {}}}", c[0], denom);
    }

    let mut inner_parts = Vec::new();
    for (i, &coeff) in c.iter().enumerate() {
        if coeff != 0 {
            let sign = if coeff < 0 {
                if inner_parts.is_empty() { "-" } else { "- " }
            } else if inner_parts.is_empty() {
                ""
            } else {
                "+ "
            };
            let abs_c = coeff.abs();
            let term_str = match i {
                0 => if abs_c == 1 { "\\frac{1}{n^7}".to_string() } else { format!("\\frac{{{abs_c}}}{{n^7}}") },
                1 => if abs_c == 1 { "\\frac{H_n}{n^6}".to_string() } else { format!("\\frac{{{abs_c}H_n}}{{n^6}}") },
                2 => if abs_c == 1 { "\\frac{H_{2n}}{n^6}".to_string() } else { format!("\\frac{{{abs_c}H_{{2n}}}}{{n^6}}") },
                3 => if abs_c == 1 { "\\frac{H_n^{(2)}}{n^5}".to_string() } else { format!("\\frac{{{abs_c}H_n^{{(2)}}}}{{n^5}}") },
                4 => if abs_c == 1 { "\\frac{H_{2n}^{(2)}}{n^5}".to_string() } else { format!("\\frac{{{abs_c}H_{{2n}}^{{(2)}}}}{{n^5}}") },
                5 => if abs_c == 1 { "\\frac{H_n^{(3)}}{n^4}".to_string() } else { format!("\\frac{{{abs_c}H_n^{{(3)}}}}{{n^4}}") },
                6 => if abs_c == 1 { "\\frac{H_n H_n^{(2)}}{n^4}".to_string() } else { format!("\\frac{{{abs_c}H_n H_n^{{(2)}}}}{{n^4}}") },
                7 => if abs_c == 1 { "\\frac{H_n^{(4)}}{n^3}".to_string() } else { format!("\\frac{{{abs_c}H_n^{{(4)}}}}{{n^3}}") },
                8 => if abs_c == 1 { "\\frac{(H_n^{(2)})^2}{n^3}".to_string() } else { format!("\\frac{{{abs_c}(H_n^{{(2)}})^2}}{{n^3}}") },
                9 => if abs_c == 1 { "\\frac{H_n H_n^{(3)}}{n^3}".to_string() } else { format!("\\frac{{{abs_c}H_n H_n^{{(3)}}}}{{n^3}}") },
                _ => String::new(),
            };
            inner_parts.push(format!("{}{}", sign, term_str));
        }
    }

    format!("\\sum_{{n=1}}^\\infty \\frac{{1}}{{{}}} \\left[ {} \\right]", denom, inner_parts.join(" "))
}

// ════════════════════════════════════════════════════════════════════════════════
//  ALAN B: KÜBİK BİNOM \binom{3n}{n} MERDİVENİ (Ağırlık-4 ve Ağırlık-5)
//  C(n) = a^n \binom{3n}{n}, a \in {2, -2}
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Debug)]
pub struct CubicBinomialBases {
    pub a: i64,
    pub w4: [DBig; 6],
    pub w5: [DBig; 7],
}

pub fn compute_cubic_bases(a: i64, steps: usize, tables: &HarmonicTables) -> CubicBinomialBases {
    let prec = PRECISION;

    let mut w4 = [
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
    ];

    let mut w5 = [
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
    ];

    let one_d = DBig::ONE.with_precision(prec).value();
    let mut a_pow = IBig::ONE;
    let a_ibig = IBig::from(a);
    let mut binom_3n = IBig::ONE;

    for n in 1..=steps {
        a_pow = &a_pow * &a_ibig;
        let num = IBig::from(3 * (3 * n as i64 - 1) * (3 * n as i64 - 2));
        let den = IBig::from(2 * n as i64 * (2 * n as i64 - 1));
        binom_3n = &binom_3n * num / den;

        let den_ibig = &a_pow * &binom_3n;
        if den_ibig.bit_len() > 720 {
            break;
        }

        let den_d = DBig::from_parts(den_ibig, 0).with_precision(prec).value();
        let u_n = &one_d / &den_d;

        let h_n = &tables.h[n];
        let h_2n = &tables.h[2 * n];
        let h_3n = &tables.h[3 * n];
        let h2_n = &tables.h2[n];
        let h3_n = &tables.h3[n];
        let h_sq = h_n * h_n;
        let h_h2 = h_n * h2_n;

        let n_d = &tables.n_dbig[n];
        let inv_n = &one_d / n_d;
        let inv_n2 = &inv_n * &inv_n;
        let inv_n3 = &inv_n2 * &inv_n;
        let inv_n4 = &inv_n2 * &inv_n2;
        let inv_n5 = &inv_n4 * &inv_n;

        // --- Ağırlık-4 Bazları (6 baz) ---
        // 0: 1 / n^4
        w4[0] = &w4[0] + &(&u_n * &inv_n4);
        // 1: H_n / n^3
        w4[1] = &w4[1] + &(&(&u_n * h_n) * &inv_n3);
        // 2: H_{2n} / n^3
        w4[2] = &w4[2] + &(&(&u_n * h_2n) * &inv_n3);
        // 3: H_{3n} / n^3
        w4[3] = &w4[3] + &(&(&u_n * h_3n) * &inv_n3);
        // 4: H_n^(2) / n^2
        w4[4] = &w4[4] + &(&(&u_n * h2_n) * &inv_n2);
        // 5: H_n^2 / n^2
        w4[5] = &w4[5] + &(&(&u_n * &h_sq) * &inv_n2);

        // --- Ağırlık-5 Bazları (7 baz) ---
        // 0: 1 / n^5
        w5[0] = &w5[0] + &(&u_n * &inv_n5);
        // 1: H_n / n^4
        w5[1] = &w5[1] + &(&(&u_n * h_n) * &inv_n4);
        // 2: H_{2n} / n^4
        w5[2] = &w5[2] + &(&(&u_n * h_2n) * &inv_n4);
        // 3: H_{3n} / n^4
        w5[3] = &w5[3] + &(&(&u_n * h_3n) * &inv_n4);
        // 4: H_n^(2) / n^3
        w5[4] = &w5[4] + &(&(&u_n * h2_n) * &inv_n3);
        // 5: H_n^(3) / n^2
        w5[5] = &w5[5] + &(&(&u_n * h3_n) * &inv_n2);
        // 6: (H_n H_n^(2)) / n^2
        w5[6] = &w5[6] + &(&(&u_n * &h_h2) * &inv_n2);
    }

    CubicBinomialBases { a, w4, w5 }
}

pub fn format_cubic_w4_lhs_latex(a: i64, c: &[i64]) -> String {
    let denom = if a < 0 {
        format!("({a})^n \\binom{{3n}}{{n}}")
    } else {
        format!("{a}^n \\binom{{3n}}{{n}}")
    };

    let has_harm = c[1..6].iter().any(|&v| v != 0);
    if !has_harm {
        return format!("\\sum_{{n=1}}^\\infty \\frac{{{}}}{{n^4 {}}}", c[0], denom);
    }

    let mut inner_parts = Vec::new();
    for (i, &coeff) in c.iter().enumerate() {
        if coeff != 0 {
            let sign = if coeff < 0 {
                if inner_parts.is_empty() { "-" } else { "- " }
            } else if inner_parts.is_empty() {
                ""
            } else {
                "+ "
            };
            let abs_c = coeff.abs();
            let term_str = match i {
                0 => if abs_c == 1 { "\\frac{1}{n^4}".to_string() } else { format!("\\frac{{{abs_c}}}{{n^4}}") },
                1 => if abs_c == 1 { "\\frac{H_n}{n^3}".to_string() } else { format!("\\frac{{{abs_c}H_n}}{{n^3}}") },
                2 => if abs_c == 1 { "\\frac{H_{2n}}{n^3}".to_string() } else { format!("\\frac{{{abs_c}H_{{2n}}}}{{n^3}}") },
                3 => if abs_c == 1 { "\\frac{H_{3n}}{n^3}".to_string() } else { format!("\\frac{{{abs_c}H_{{3n}}}}{{n^3}}") },
                4 => if abs_c == 1 { "\\frac{H_n^{(2)}}{n^2}".to_string() } else { format!("\\frac{{{abs_c}H_n^{{(2)}}}}{{n^2}}") },
                5 => if abs_c == 1 { "\\frac{H_n^2}{n^2}".to_string() } else { format!("\\frac{{{abs_c}H_n^2}}{{n^2}}") },
                _ => String::new(),
            };
            inner_parts.push(format!("{}{}", sign, term_str));
        }
    }

    format!("\\sum_{{n=1}}^\\infty \\frac{{1}}{{{}}} \\left[ {} \\right]", denom, inner_parts.join(" "))
}

pub fn format_cubic_w5_lhs_latex(a: i64, c: &[i64]) -> String {
    let denom = if a < 0 {
        format!("({a})^n \\binom{{3n}}{{n}}")
    } else {
        format!("{a}^n \\binom{{3n}}{{n}}")
    };

    let has_harm = c[1..7].iter().any(|&v| v != 0);
    if !has_harm {
        return format!("\\sum_{{n=1}}^\\infty \\frac{{{}}}{{n^5 {}}}", c[0], denom);
    }

    let mut inner_parts = Vec::new();
    for (i, &coeff) in c.iter().enumerate() {
        if coeff != 0 {
            let sign = if coeff < 0 {
                if inner_parts.is_empty() { "-" } else { "- " }
            } else if inner_parts.is_empty() {
                ""
            } else {
                "+ "
            };
            let abs_c = coeff.abs();
            let term_str = match i {
                0 => if abs_c == 1 { "\\frac{1}{n^5}".to_string() } else { format!("\\frac{{{abs_c}}}{{n^5}}") },
                1 => if abs_c == 1 { "\\frac{H_n}{n^4}".to_string() } else { format!("\\frac{{{abs_c}H_n}}{{n^4}}") },
                2 => if abs_c == 1 { "\\frac{H_{2n}}{n^4}".to_string() } else { format!("\\frac{{{abs_c}H_{{2n}}}}{{n^4}}") },
                3 => if abs_c == 1 { "\\frac{H_{3n}}{n^4}".to_string() } else { format!("\\frac{{{abs_c}H_{{3n}}}}{{n^4}}") },
                4 => if abs_c == 1 { "\\frac{H_n^{(2)}}{n^3}".to_string() } else { format!("\\frac{{{abs_c}H_n^{{(2)}}}}{{n^3}}") },
                5 => if abs_c == 1 { "\\frac{H_n^{(3)}}{n^2}".to_string() } else { format!("\\frac{{{abs_c}H_n^{{(3)}}}}{{n^2}}") },
                6 => if abs_c == 1 { "\\frac{H_n H_n^{(2)}}{n^2}".to_string() } else { format!("\\frac{{{abs_c}H_n H_n^{{(2)}}}}{{n^2}}") },
                _ => String::new(),
            };
            inner_parts.push(format!("{}{}", sign, term_str));
        }
    }

    format!("\\sum_{{n=1}}^\\infty \\frac{{1}}{{{}}} \\left[ {} \\right]", denom, inner_parts.join(" "))
}

// ════════════════════════════════════════════════════════════════════════════════
//  ALAN C: TEKİL PAYDALI ÇEKİRDEKLER ((2n+1)^k ve \beta(4))
//  D(n) = a^n (2n+1)^k \binom{2n}{n}, a \in {-4, -1, 1, 4}, k \in {1, 2, 3, 4}
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Debug)]
pub struct OddDenomBases {
    pub a: i64,
    pub k: usize,
    pub bases: [DBig; 4],
}

pub fn compute_odd_denom_bases(a: i64, k: usize, steps: usize, tables: &HarmonicTables) -> OddDenomBases {
    let prec = PRECISION;
    let mut bases = [
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
        DBig::ZERO.with_precision(prec).value(),
    ];

    let one_d = DBig::ONE.with_precision(prec).value();

    // n = 0: a^0 * 1^k * binom(0,0) = 1.
    // O_0 = 0, H_0 = 0, O_0^(2) = 0.
    bases[0] = &bases[0] + &one_d;

    let mut a_pow = IBig::ONE;
    let a_ibig = IBig::from(a);
    let mut binom_2n = IBig::ONE;

    for n in 1..=steps {
        a_pow = &a_pow * &a_ibig;
        binom_2n = &binom_2n * IBig::from(2 * (2 * n as i64 - 1)) / IBig::from(n as i64);

        let odd_val = 2 * n as i64 + 1;
        let odd_i = IBig::from(odd_val);
        let mut odd_k = IBig::ONE;
        for _ in 0..k {
            odd_k = &odd_k * &odd_i;
        }

        let den_ibig = &(&a_pow * &odd_k) * &binom_2n;
        if den_ibig.bit_len() > 720 {
            break;
        }

        let den_d = DBig::from_parts(den_ibig, 0).with_precision(prec).value();
        let u_n = &one_d / &den_d;

        let o_n = &tables.o[n];
        let h_n = &tables.h[n];
        let o2_n = &tables.o2[n];

        // 0: 1 / D(n)
        bases[0] = &bases[0] + &u_n;
        // 1: O_n / D(n)
        bases[1] = &bases[1] + &(&u_n * o_n);
        // 2: H_n / D(n)
        bases[2] = &bases[2] + &(&u_n * h_n);
        // 3: O_n^(2) / D(n)
        bases[3] = &bases[3] + &(&u_n * o2_n);
    }

    OddDenomBases { a, k, bases }
}

pub fn format_odd_denom_lhs_latex(a: i64, k: usize, c: &[i64]) -> String {
    let a_str = match a {
        1 => String::new(),
        -1 => "(-1)^n ".to_string(),
        val if val < 0 => format!("({val})^n "),
        val => format!("{val}^n "),
    };

    let k_str = if k == 1 {
        "(2n+1)".to_string()
    } else {
        format!("(2n+1)^{k}")
    };

    let denom = format!("{a_str}{k_str} \\binom{{2n}}{{n}}");

    let has_harm = c[1..4].iter().any(|&v| v != 0);
    if !has_harm {
        return format!("\\sum_{{n=0}}^\\infty \\frac{{{}}}{{{}}}", c[0], denom);
    }

    let mut inner_parts = Vec::new();
    for (i, &coeff) in c.iter().enumerate() {
        if coeff != 0 {
            let sign = if coeff < 0 {
                if inner_parts.is_empty() { "-" } else { "- " }
            } else if inner_parts.is_empty() {
                ""
            } else {
                "+ "
            };
            let abs_c = coeff.abs();
            let term_str = match i {
                0 => format!("{abs_c}"),
                1 => if abs_c == 1 { "O_n".to_string() } else { format!("{abs_c} O_n") },
                2 => if abs_c == 1 { "H_n".to_string() } else { format!("{abs_c} H_n") },
                3 => if abs_c == 1 { "O_n^{(2)}".to_string() } else { format!("{abs_c} O_n^{{(2)}}") },
                _ => String::new(),
            };
            inner_parts.push(format!("{}{}", sign, term_str));
        }
    }

    format!("\\sum_{{n=0}}^\\infty \\frac{{1}}{{{}}} \\left[ {} \\right]", denom, inner_parts.join(" "))
}

// ════════════════════════════════════════════════════════════════════════════════
//  GENEL ÇEKİRDEK LİSTESİ VE GERİYE UYUMLULUK
// ════════════════════════════════════════════════════════════════════════════════

pub fn generate_all_kernels() -> Vec<SeriesKernel> {
    let a_values = [-16, -8, -4, -2, -1, 1, 2, 3, 4, 8, 16];
    let ktypes = [
        KernelType::CentralBinomial,
        KernelType::SquaredCentralBinomial,
        KernelType::CubicBinomial,
        KernelType::QuarticBinomial,
    ];

    let mut kernels = Vec::with_capacity(a_values.len() * ktypes.len());
    for &ktype in &ktypes {
        for &a in &a_values {
            kernels.push(SeriesKernel::new(a, ktype));
        }
    }
    kernels
}

// ════════════════════════════════════════════════════════════════════════════════
//  RAMANUJAN / GUILLERA SÜPER-BİNOM ÇEKİRDEKLERİ VE POLİNOM BAZLARI
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SuperBinomialType {
    B1, // binom(2n, n)^3
    B2, // binom(2n, n)^4
    B3, // binom(2n, n)^5
    B4, // binom(2n, n)^2 * binom(3n, n)
    B5, // binom(2n, n) * binom(4n, 2n)
    B6, // binom(2n, n)^2 * binom(4n, 2n)
    B7, // binom(3n, n)^2
}

impl SuperBinomialType {
    pub const ALL: [SuperBinomialType; 7] = [
        SuperBinomialType::B1,
        SuperBinomialType::B2,
        SuperBinomialType::B3,
        SuperBinomialType::B4,
        SuperBinomialType::B5,
        SuperBinomialType::B6,
        SuperBinomialType::B7,
    ];

    pub fn code(&self) -> &'static str {
        match self {
            Self::B1 => "B_1",
            Self::B2 => "B_2",
            Self::B3 => "B_3",
            Self::B4 => "B_4",
            Self::B5 => "B_5",
            Self::B6 => "B_6",
            Self::B7 => "B_7",
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Self::B1 => "B_1(n) = binom(2n,n)^3",
            Self::B2 => "B_2(n) = binom(2n,n)^4",
            Self::B3 => "B_3(n) = binom(2n,n)^5",
            Self::B4 => "B_4(n) = binom(2n,n)^2 * binom(3n,n)",
            Self::B5 => "B_5(n) = binom(2n,n) * binom(4n,2n)",
            Self::B6 => "B_6(n) = binom(2n,n)^2 * binom(4n,2n)",
            Self::B7 => "B_7(n) = binom(3n,n)^2",
        }
    }

    pub fn latex(&self) -> &'static str {
        match self {
            Self::B1 => "\\binom{2n}{n}^3",
            Self::B2 => "\\binom{2n}{n}^4",
            Self::B3 => "\\binom{2n}{n}^5",
            Self::B4 => "\\binom{2n}{n}^2 \\binom{3n}{n}",
            Self::B5 => "\\binom{2n}{n} \\binom{4n}{2n}",
            Self::B6 => "\\binom{2n}{n}^2 \\binom{4n}{2n}",
            Self::B7 => "\\binom{3n}{n}^2",
        }
    }

    pub fn growth_ratio(&self) -> f64 {
        match self {
            Self::B1 => 64.0,
            Self::B2 => 256.0,
            Self::B3 => 1024.0,
            Self::B4 => 108.0,
            Self::B5 => 64.0,
            Self::B6 => 256.0,
            Self::B7 => 729.0 / 16.0, // 45.5625
        }
    }
}

pub const STANDARD_BASES: [i64; 26] = [
    16, -16, 27, -27, 32, -32, 64, -64, 81, -81, 108, -108, 128, -128,
    216, -216, 256, -256, 512, -512, 1024, -1024, 4096, -4096, 8192, -8192,
];

pub const GUILLERA_5_BASES: [i64; 2] = [1048576, -1048576]; // 2^20, -2^20

pub fn get_bases_for_kernel(kernel: SuperBinomialType) -> Vec<i64> {
    let mut bases = STANDARD_BASES.to_vec();
    if kernel == SuperBinomialType::B3 {
        bases.extend_from_slice(&GUILLERA_5_BASES);
    }
    bases
}

pub fn compute_super_binomial_bases(
    kernel: SuperBinomialType,
    a: i64,
    prec: usize,
    max_steps: usize,
) -> Option<[DBig; 4]> {
    let growth = kernel.growth_ratio();
    if (a.abs() as f64) <= growth + 0.1 {
        return None;
    }

    let prec_calc = prec + 25;
    let mut s0 = DBig::ZERO.with_precision(prec_calc).value();
    let mut s1 = DBig::ZERO.with_precision(prec_calc).value();
    let mut s2 = DBig::ZERO.with_precision(prec_calc).value();
    let mut s3 = DBig::ZERO.with_precision(prec_calc).value();

    // n = 0 term: B(0) = 1, a^0 = 1 => 1
    s0 = &s0 + DBig::ONE.with_precision(prec_calc).value();

    let mut c2n = IBig::ONE;
    let mut c3n = IBig::ONE;
    let mut c4n = IBig::ONE;
    let mut a_pow = IBig::ONE;
    let a_ibig = IBig::from(a);

    let eps = DBig::from_parts(IBig::ONE, -215).with_precision(prec_calc).value();

    for n in 1..=max_steps {
        let n_i = n as i64;
        c2n = &c2n * IBig::from(2 * (2 * n_i - 1)) / IBig::from(n_i);
        c3n = &c3n * IBig::from(3 * (3 * n_i - 1) * (3 * n_i - 2))
            / IBig::from(2 * n_i * (2 * n_i - 1));
        c4n = &c4n * IBig::from(2 * (4 * n_i - 1) * (4 * n_i - 3))
            / IBig::from(n_i * (2 * n_i - 1));

        a_pow = &a_pow * &a_ibig;

        let b_n: IBig = match kernel {
            SuperBinomialType::B1 => &c2n * &c2n * &c2n,
            SuperBinomialType::B2 => {
                let sq = &c2n * &c2n;
                &sq * &sq
            }
            SuperBinomialType::B3 => {
                let sq = &c2n * &c2n;
                let q = &sq * &sq;
                &q * &c2n
            }
            SuperBinomialType::B4 => {
                let sq = &c2n * &c2n;
                &sq * &c3n
            }
            SuperBinomialType::B5 => &c2n * &c4n,
            SuperBinomialType::B6 => {
                let sq = &c2n * &c2n;
                &sq * &c4n
            }
            SuperBinomialType::B7 => &c3n * &c3n,
        };

        let b_d = DBig::from_parts(b_n, 0).with_precision(prec_calc).value();
        let a_d = DBig::from_parts(a_pow.clone(), 0).with_precision(prec_calc).value();
        let term = (&b_d / &a_d).with_precision(prec_calc).value();

        if n == 15 {
            let abs_t = if term < DBig::ZERO { -term.clone() } else { term.clone() };
            if abs_t >= DBig::ONE.with_precision(prec_calc).value() {
                return None;
            }
        }

        let n_d = DBig::from_parts(IBig::from(n_i), 0).with_precision(prec_calc).value();
        let n2_d = (&n_d * &n_d).with_precision(prec_calc).value();
        let n3_d = (&n2_d * &n_d).with_precision(prec_calc).value();

        s0 = &s0 + &term;
        s1 = &s1 + &(&n_d * &term);
        s2 = &s2 + &(&n2_d * &term);
        s3 = &s3 + &(&n3_d * &term);

        let abs_t = if term < DBig::ZERO { -term } else { term };
        if abs_t < eps {
            break;
        }
    }

    Some([
        s0.with_precision(prec).value(),
        s1.with_precision(prec).value(),
        s2.with_precision(prec).value(),
        s3.with_precision(prec).value(),
    ])
}

pub fn format_polynomial_text(coeffs: &[i64; 4]) -> String {
    let [c0, c1, c2, c3] = *coeffs;
    let mut terms = Vec::new();

    if c3 != 0 {
        terms.push((c3, "n^3"));
    }
    if c2 != 0 {
        terms.push((c2, "n^2"));
    }
    if c1 != 0 {
        terms.push((c1, "n"));
    }
    if c0 != 0 {
        terms.push((c0, ""));
    }

    if terms.is_empty() {
        return "0".to_string();
    }

    let mut res = String::new();
    for (i, (c, var)) in terms.iter().enumerate() {
        let abs_c = c.abs();
        let sign = if *c < 0 {
            if i == 0 { "-" } else { " - " }
        } else if i == 0 {
            ""
        } else {
            " + "
        };

        let term_str = if var.is_empty() {
            format!("{abs_c}")
        } else if abs_c == 1 {
            var.to_string()
        } else {
            format!("{abs_c}*{var}")
        };
        res.push_str(&format!("{sign}{term_str}"));
    }
    res
}

pub fn format_polynomial_latex(coeffs: &[i64; 4]) -> String {
    let [c0, c1, c2, c3] = *coeffs;
    let mut terms = Vec::new();

    if c3 != 0 {
        terms.push((c3, "n^3"));
    }
    if c2 != 0 {
        terms.push((c2, "n^2"));
    }
    if c1 != 0 {
        terms.push((c1, "n"));
    }
    if c0 != 0 {
        terms.push((c0, ""));
    }

    if terms.is_empty() {
        return "0".to_string();
    }

    let mut res = String::new();
    for (i, (c, var)) in terms.iter().enumerate() {
        let abs_c = c.abs();
        let sign = if *c < 0 {
            if i == 0 { "-" } else { " - " }
        } else if i == 0 {
            ""
        } else {
            " + "
        };

        let term_str = if var.is_empty() {
            format!("{abs_c}")
        } else if abs_c == 1 {
            var.to_string()
        } else {
            format!("{abs_c} {var}")
        };
        res.push_str(&format!("{sign}{term_str}"));
    }
    res
}

// ════════════════════════════════════════════════════════════════════════════════
//  SINIR ALAN ÇEKİRDEKLERİ (DOMB, ALMKVIST-ZUDILIN, B8, B9)
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FrontierKernelType {
    Domb,            // D_n
    AlmkvistZudilin, // Z_n
    B8,              // binom(2n, n)^6
    B9,              // binom(2n, n)^3 * binom(3n, n)^2
}

impl FrontierKernelType {
    pub const ALL: [FrontierKernelType; 4] = [
        FrontierKernelType::Domb,
        FrontierKernelType::AlmkvistZudilin,
        FrontierKernelType::B8,
        FrontierKernelType::B9,
    ];

    pub fn code(&self) -> &'static str {
        match self {
            Self::Domb => "D",
            Self::AlmkvistZudilin => "Z",
            Self::B8 => "B_8",
            Self::B9 => "B_9",
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Self::Domb => "Domb Dizisi D_n = (-1)^n sum_{k=0}^n binom(n,k)^2 binom(2k,k) binom(2n-2k, n-k)",
            Self::AlmkvistZudilin => "Almkvist-Zudilin Z_n = sum_{k=0}^{n/3} (-1)^k 3^{n-3k} (3k)!/(k!)^3 binom(n,3k) binom(n+k,k)",
            Self::B8 => "B_8(n) = binom(2n,n)^6 (6. Mertebe Süper-Binom)",
            Self::B9 => "B_9(n) = binom(2n,n)^3 * binom(3n,n)^2 (Kübik-Merkezi Kuantik Melez)",
        }
    }

    pub fn latex(&self) -> &'static str {
        match self {
            Self::Domb => "D_n",
            Self::AlmkvistZudilin => "Z_n",
            Self::B8 => "\\binom{2n}{n}^6",
            Self::B9 => "\\binom{2n}{n}^3 \\binom{3n}{n}^2",
        }
    }

    pub fn growth_ratio(&self) -> f64 {
        match self {
            Self::Domb => 16.0,
            Self::AlmkvistZudilin => 27.0,
            Self::B8 => 4096.0,
            Self::B9 => 2916.0,
        }
    }
}

pub const DOMB_BASES: [i64; 28] = [
    -32, 32, -64, 64, -128, 128, -256, 256, -288, 288, -300, 300, -576, 576, -768, 768,
    -1080, 1080, -1584, 1584, -2880, 2880, -3969, 3969, -5832, 5832, -11664, 11664,
];

pub const AZ_BASES: [i64; 28] = [
    -72, 72, -81, 81, -108, 108, -216, 216, -324, 324, -432, 432, -576, 576, -864, 864,
    -1296, 1296, -1728, 1728, -2916, 2916, -3888, 3888, -5184, 5184, -46656, 46656,
];

pub const B8_B9_BASES: [i64; 18] = [
    -8192, 8192, -32768, 32768, -65536, 65536, -262144, 262144, -1048576, 1048576, -2097152,
    2097152, -4194304, 4194304, -16777216, 16777216, -67108864, 67108864,
];

pub fn get_frontier_bases(kernel: FrontierKernelType) -> &'static [i64] {
    match kernel {
        FrontierKernelType::Domb => &DOMB_BASES,
        FrontierKernelType::AlmkvistZudilin => &AZ_BASES,
        FrontierKernelType::B8 | FrontierKernelType::B9 => &B8_B9_BASES,
    }
}

pub fn compute_domb_table(max_n: usize) -> Vec<IBig> {
    let mut central = Vec::with_capacity(max_n + 1);
    let mut c = IBig::ONE;
    central.push(c.clone());
    for m in 1..=max_n {
        c = &c * IBig::from(2 * (2 * m as i64 - 1)) / IBig::from(m as i64);
        central.push(c.clone());
    }

    let mut pascal: Vec<Vec<IBig>> = Vec::with_capacity(max_n + 1);
    for n in 0..=max_n {
        let mut row = Vec::with_capacity(n + 1);
        row.push(IBig::ONE);
        for k in 1..n {
            let val = &pascal[n - 1][k - 1] + &pascal[n - 1][k];
            row.push(val);
        }
        if n > 0 {
            row.push(IBig::ONE);
        }
        pascal.push(row);
    }

    let mut domb = Vec::with_capacity(max_n + 1);
    for n in 0..=max_n {
        let mut sum = IBig::ZERO;
        for k in 0..=n {
            let b_nk = &pascal[n][k];
            let b_nk_sq = b_nk * b_nk;
            let c_k = &central[k];
            let c_nmk = &central[n - k];
            let term = (&b_nk_sq * c_k) * c_nmk;
            sum += term;
        }
        if n % 2 == 1 {
            sum = -sum;
        }
        domb.push(sum);
    }
    domb
}

pub fn compute_az_table(max_n: usize) -> Vec<IBig> {
    let max_pascal = max_n + (max_n / 3) + 5;
    let mut pascal: Vec<Vec<IBig>> = Vec::with_capacity(max_pascal + 1);
    for n in 0..=max_pascal {
        let mut row = Vec::with_capacity(n + 1);
        row.push(IBig::ONE);
        for k in 1..n {
            let val = &pascal[n - 1][k - 1] + &pascal[n - 1][k];
            row.push(val);
        }
        if n > 0 {
            row.push(IBig::ONE);
        }
        pascal.push(row);
    }

    let max_k = max_n / 3 + 2;
    let mut w = Vec::with_capacity(max_k + 1);
    for k in 0..=max_k {
        let val = &pascal[3 * k][k] * &pascal[2 * k][k];
        w.push(val);
    }

    let mut powers_of_3 = Vec::with_capacity(max_n + 1);
    let mut p = IBig::ONE;
    powers_of_3.push(p.clone());
    for _ in 1..=max_n {
        p = &p * IBig::from(3);
        powers_of_3.push(p.clone());
    }

    let mut az = Vec::with_capacity(max_n + 1);
    for n in 0..=max_n {
        let mut sum = IBig::ZERO;
        let k_limit = n / 3;
        for k in 0..=k_limit {
            let p3 = &powers_of_3[n - 3 * k];
            let w_k = &w[k];
            let b1 = &pascal[n][3 * k];
            let b2 = &pascal[n + k][k];
            let term = ((p3 * w_k) * b1) * b2;
            if k % 2 == 0 {
                sum += term;
            } else {
                sum -= term;
            }
        }
        az.push(sum);
    }
    az
}

pub fn compute_frontier_bases(
    kernel: FrontierKernelType,
    a: i64,
    prec: usize,
    max_steps: usize,
    domb_table: &[IBig],
    az_table: &[IBig],
) -> Option<[DBig; 5]> {
    let growth = kernel.growth_ratio();
    if (a.abs() as f64) <= growth + 0.1 {
        return None;
    }

    let prec_calc = prec + 25;
    let mut s0 = DBig::ZERO.with_precision(prec_calc).value();
    let mut s1 = DBig::ZERO.with_precision(prec_calc).value();
    let mut s2 = DBig::ZERO.with_precision(prec_calc).value();
    let mut s3 = DBig::ZERO.with_precision(prec_calc).value();
    let mut s4 = DBig::ZERO.with_precision(prec_calc).value();

    // n = 0 term: K(0) = 1, a^0 = 1 => 1
    s0 = &s0 + DBig::ONE.with_precision(prec_calc).value();

    let mut c2n = IBig::ONE;
    let mut c3n = IBig::ONE;
    let mut a_pow = IBig::ONE;
    let a_ibig = IBig::from(a);

    let eps = DBig::from_parts(IBig::ONE, -235).with_precision(prec_calc).value();

    for n in 1..=max_steps {
        let n_i = n as i64;
        a_pow = &a_pow * &a_ibig;

        let k_n: IBig = match kernel {
            FrontierKernelType::Domb => {
                if n < domb_table.len() {
                    domb_table[n].clone()
                } else {
                    return None;
                }
            }
            FrontierKernelType::AlmkvistZudilin => {
                if n < az_table.len() {
                    az_table[n].clone()
                } else {
                    return None;
                }
            }
            FrontierKernelType::B8 => {
                c2n = &c2n * IBig::from(2 * (2 * n_i - 1)) / IBig::from(n_i);
                let sq = &c2n * &c2n;
                (&sq * &sq) * &sq
            }
            FrontierKernelType::B9 => {
                c2n = &c2n * IBig::from(2 * (2 * n_i - 1)) / IBig::from(n_i);
                c3n = &c3n * IBig::from(3 * (3 * n_i - 1) * (3 * n_i - 2))
                    / IBig::from(2 * n_i * (2 * n_i - 1));
                let c2_3 = (&c2n * &c2n) * &c2n;
                let c3_2 = &c3n * &c3n;
                &c2_3 * &c3_2
            }
        };

        let k_d = DBig::from_parts(k_n, 0).with_precision(prec_calc).value();
        let a_d = DBig::from_parts(a_pow.clone(), 0).with_precision(prec_calc).value();
        let term = (&k_d / &a_d).with_precision(prec_calc).value();

        let abs_t = if term < DBig::ZERO { -term.clone() } else { term.clone() };
        if n > 50 && abs_t > DBig::from_parts(IBig::from(1000000), 0).with_precision(prec_calc).value() {
            return None;
        }

        let n_d = DBig::from_parts(IBig::from(n_i), 0).with_precision(prec_calc).value();
        let n2_d = (&n_d * &n_d).with_precision(prec_calc).value();
        let n3_d = (&n2_d * &n_d).with_precision(prec_calc).value();
        let n4_d = (&n2_d * &n2_d).with_precision(prec_calc).value();

        s0 = &s0 + &term;
        s1 = &s1 + &(&n_d * &term);
        s2 = &s2 + &(&n2_d * &term);
        s3 = &s3 + &(&n3_d * &term);
        s4 = &s4 + &(&n4_d * &term);

        let abs_t = if term < DBig::ZERO { -term } else { term };
        if abs_t < eps {
            break;
        }
    }

    Some([
        s0.with_precision(prec).value(),
        s1.with_precision(prec).value(),
        s2.with_precision(prec).value(),
        s3.with_precision(prec).value(),
        s4.with_precision(prec).value(),
    ])
}

pub fn format_polynomial5_text(coeffs: &[i64; 5]) -> String {
    let [c0, c1, c2, c3, c4] = *coeffs;
    let mut terms = Vec::new();

    if c4 != 0 {
        terms.push((c4, "n^4"));
    }
    if c3 != 0 {
        terms.push((c3, "n^3"));
    }
    if c2 != 0 {
        terms.push((c2, "n^2"));
    }
    if c1 != 0 {
        terms.push((c1, "n"));
    }
    if c0 != 0 {
        terms.push((c0, ""));
    }

    if terms.is_empty() {
        return "0".to_string();
    }

    let mut res = String::new();
    for (i, (c, var)) in terms.iter().enumerate() {
        let abs_c = c.abs();
        let sign = if *c < 0 {
            if i == 0 { "-" } else { " - " }
        } else if i == 0 {
            ""
        } else {
            " + "
        };

        let term_str = if var.is_empty() {
            format!("{abs_c}")
        } else if abs_c == 1 {
            var.to_string()
        } else {
            format!("{abs_c}*{var}")
        };
        res.push_str(&format!("{sign}{term_str}"));
    }
    res
}

pub fn format_polynomial5_latex(coeffs: &[i64; 5]) -> String {
    let [c0, c1, c2, c3, c4] = *coeffs;
    let mut terms = Vec::new();

    if c4 != 0 {
        terms.push((c4, "n^4"));
    }
    if c3 != 0 {
        terms.push((c3, "n^3"));
    }
    if c2 != 0 {
        terms.push((c2, "n^2"));
    }
    if c1 != 0 {
        terms.push((c1, "n"));
    }
    if c0 != 0 {
        terms.push((c0, ""));
    }

    if terms.is_empty() {
        return "0".to_string();
    }

    let mut res = String::new();
    for (i, (c, var)) in terms.iter().enumerate() {
        let abs_c = c.abs();
        let sign = if *c < 0 {
            if i == 0 { "-" } else { " - " }
        } else if i == 0 {
            ""
        } else {
            " + "
        };

        let term_str = if var.is_empty() {
            format!("{abs_c}")
        } else if abs_c == 1 {
            var.to_string()
        } else {
            format!("{abs_c} {var}")
        };
        res.push_str(&format!("{sign}{term_str}"));
    }
    res
}

// ════════════════════════════════════════════════════════════════════════════════
//  4. TERS CALABI-YAU / SÜPER-BİNOM KÖŞEGEN MOTORU (B1..B14)
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum InverseKernelType {
    B1,  // binom(2n, n) * binom(3n, n)
    B2,  // binom(2n, n)^2 * binom(3n, n)
    B3,  // binom(2n, n) * binom(3n, n)^2
    B4,  // binom(2n, n)^2 * binom(3n, n)^2
    B5,  // binom(2n, n) * binom(4n, 2n)
    B6,  // binom(2n, n)^2 * binom(4n, 2n)
    B7,  // binom(2n, n) * binom(3n, n) * binom(6n, 3n)
    B8,  // binom(2n, n)^2
    B9,  // binom(2n, n)^3
    B10, // binom(2n, n)^4
    B11, // binom(2n, n)^5
    B12, // binom(2n, n)^7
    B13, // binom(3n, n)^2
    B14, // binom(3n, n)^3
}

impl InverseKernelType {
    pub const ALL: [InverseKernelType; 14] = [
        Self::B1, Self::B2, Self::B3, Self::B4, Self::B5, Self::B6, Self::B7,
        Self::B8, Self::B9, Self::B10, Self::B11, Self::B12, Self::B13, Self::B14,
    ];

    pub fn code(&self) -> &'static str {
        match self {
            Self::B1 => "B_1",
            Self::B2 => "B_2",
            Self::B3 => "B_3",
            Self::B4 => "B_4",
            Self::B5 => "B_5",
            Self::B6 => "B_6",
            Self::B7 => "B_7",
            Self::B8 => "B_8",
            Self::B9 => "B_9",
            Self::B10 => "B_{10}",
            Self::B11 => "B_{11}",
            Self::B12 => "B_{12}",
            Self::B13 => "B_{13}",
            Self::B14 => "B_{14}",
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Self::B1 => "Deuring Kübik B_1(n) = binom(2n,n)*binom(3n,n)",
            Self::B2 => "Kübik-Merkezi Melez B_2(n) = binom(2n,n)^2*binom(3n,n)",
            Self::B3 => "Çift Kübik-Merkezi B_3(n) = binom(2n,n)*binom(3n,n)^2",
            Self::B4 => "Seviye-3 Calabi-Yau B_4(n) = binom(2n,n)^2*binom(3n,n)^2",
            Self::B5 => "Kuantik-Merkezi İkiz B_5(n) = binom(2n,n)*binom(4n,2n)",
            Self::B6 => "Quintic Calabi-Yau B_6(n) = binom(2n,n)^2*binom(4n,2n)",
            Self::B7 => "Sextic Calabi-Yau B_7(n) = binom(2n,n)*binom(3n,n)*binom(6n,3n)",
            Self::B8 => "Merkezi Kare B_8(n) = binom(2n,n)^2",
            Self::B9 => "Merkezi Küp B_9(n) = binom(2n,n)^3",
            Self::B10 => "Merkezi 4. Güç B_{10}(n) = binom(2n,n)^4",
            Self::B11 => "Merkezi 5. Güç B_{11}(n) = binom(2n,n)^5",
            Self::B12 => "Merkezi 7. Güç B_{12}(n) = binom(2n,n)^7",
            Self::B13 => "Saf Kübik Kare B_{13}(n) = binom(3n,n)^2",
            Self::B14 => "Saf Kübik Küp B_{14}(n) = binom(3n,n)^3",
        }
    }

    pub fn latex(&self) -> &'static str {
        match self {
            Self::B1 => "\\binom{2n}{n} \\binom{3n}{n}",
            Self::B2 => "\\binom{2n}{n}^2 \\binom{3n}{n}",
            Self::B3 => "\\binom{2n}{n} \\binom{3n}{n}^2",
            Self::B4 => "\\binom{2n}{n}^2 \\binom{3n}{n}^2",
            Self::B5 => "\\binom{2n}{n} \\binom{4n}{2n}",
            Self::B6 => "\\binom{2n}{n}^2 \\binom{4n}{2n}",
            Self::B7 => "\\binom{2n}{n} \\binom{3n}{n} \\binom{6n}{3n}",
            Self::B8 => "\\binom{2n}{n}^2",
            Self::B9 => "\\binom{2n}{n}^3",
            Self::B10 => "\\binom{2n}{n}^4",
            Self::B11 => "\\binom{2n}{n}^5",
            Self::B12 => "\\binom{2n}{n}^7",
            Self::B13 => "\\binom{3n}{n}^2",
            Self::B14 => "\\binom{3n}{n}^3",
        }
    }

    pub fn growth_ratio(&self) -> f64 {
        match self {
            Self::B1 => 27.0,
            Self::B2 => 108.0,
            Self::B3 => 182.25,
            Self::B4 => 729.0,
            Self::B5 => 64.0,
            Self::B6 => 256.0,
            Self::B7 => 1728.0,
            Self::B8 => 16.0,
            Self::B9 => 64.0,
            Self::B10 => 256.0,
            Self::B11 => 1024.0,
            Self::B12 => 16384.0,
            Self::B13 => 45.5625,
            Self::B14 => 307.546875,
        }
    }

    pub fn has_cubic_pole(&self) -> bool {
        matches!(
            self,
            Self::B1 | Self::B2 | Self::B3 | Self::B4 | Self::B7 | Self::B13 | Self::B14
        )
    }

    pub fn eval_at(&self, n: usize, tables: &BinomialTables) -> IBig {
        let c2 = &tables.c2n[n];
        let c3 = &tables.c3n[n];
        let c4 = &tables.c4n[n];
        let c6 = &tables.c6n[n];

        match self {
            Self::B1 => c2 * c3,
            Self::B2 => (c2 * c2) * c3,
            Self::B3 => c2 * (c3 * c3),
            Self::B4 => (c2 * c2) * (c3 * c3),
            Self::B5 => c2 * c4,
            Self::B6 => (c2 * c2) * c4,
            Self::B7 => (c2 * c3) * c6,
            Self::B8 => c2 * c2,
            Self::B9 => (c2 * c2) * c2,
            Self::B10 => {
                let sq = c2 * c2;
                &sq * &sq
            }
            Self::B11 => {
                let sq = c2 * c2;
                (&sq * &sq) * c2
            }
            Self::B12 => {
                let sq = c2 * c2;
                let q4 = &sq * &sq;
                (&q4 * &sq) * c2
            }
            Self::B13 => c3 * c3,
            Self::B14 => (c3 * c3) * c3,
        }
    }
}

pub const INVERSE_CANDIDATE_BASES: [i64; 26] = [
    -1, 1, -2, 2, -3, 3, -4, 4, -8, 8, -9, 9, -12, 12, -16, 16, -27, 27, -32, 32, -64, 64,
    -81, 81, -256, 256,
];

pub struct BinomialTables {
    pub c2n: Vec<IBig>,
    pub c3n: Vec<IBig>,
    pub c4n: Vec<IBig>,
    pub c6n: Vec<IBig>,
}

pub fn compute_binomial_tables(max_n: usize) -> BinomialTables {
    let mut c2n = Vec::with_capacity(max_n + 1);
    let mut c3n = Vec::with_capacity(max_n + 1);
    let mut c4n = Vec::with_capacity(max_n + 1);
    let mut c6n = Vec::with_capacity(max_n + 1);

    let mut b2 = IBig::ONE;
    let mut b3 = IBig::ONE;
    let mut b4 = IBig::ONE;
    let mut b6 = IBig::ONE;

    c2n.push(b2.clone());
    c3n.push(b3.clone());
    c4n.push(b4.clone());
    c6n.push(b6.clone());

    for n in 1..=max_n {
        let n_i = n as i64;
        b2 = &b2 * IBig::from(2 * (2 * n_i - 1)) / IBig::from(n_i);
        b3 = &b3 * IBig::from(3 * (3 * n_i - 1) * (3 * n_i - 2))
            / IBig::from(2 * n_i * (2 * n_i - 1));
        b4 = &b4 * IBig::from(2 * (4 * n_i - 1) * (4 * n_i - 3))
            / IBig::from(n_i * (2 * n_i - 1));
        b6 = &b6 * IBig::from(2 * (6 * n_i - 1) * (6 * n_i - 3) * (6 * n_i - 5))
            / IBig::from(n_i * (3 * n_i - 1) * (3 * n_i - 2));

        c2n.push(b2.clone());
        c3n.push(b3.clone());
        c4n.push(b4.clone());
        c6n.push(b6.clone());
    }

    BinomialTables { c2n, c3n, c4n, c6n }
}

/// Computes cyclotomic pole basis vectors U_k^(m) for k = 1..6:
/// U_k^(m) = sum_{n=1}^\infty a^n / (n^k * (2n-1)^m * B_i(n))
pub fn compute_inverse_cyclotomic_bases(
    kernel: InverseKernelType,
    a: i64,
    m: usize,
    prec: usize,
    max_steps: usize,
    tables: &BinomialTables,
) -> Option<[DBig; 6]> {
    let growth = kernel.growth_ratio();
    if (a.abs() as f64) > (growth / 2.5) + 1e-6 {
        return None;
    }

    let prec_calc = prec + 30;
    let mut u = [
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
    ];

    let mut a_pow = IBig::ONE;
    let a_ibig = IBig::from(a);
    let eps = DBig::from_parts(IBig::ONE, -240).with_precision(prec_calc).value();

    for n in 1..=max_steps {
        let n_i = n as i64;
        a_pow = &a_pow * &a_ibig;

        let b_val = kernel.eval_at(n, tables);

        let mut pole_factor = IBig::ONE;
        if m > 0 {
            let odd = IBig::from(2 * n_i - 1);
            for _ in 0..m {
                pole_factor = &pole_factor * &odd;
            }
        }
        let den_i = &b_val * &pole_factor;

        let a_d = DBig::from_parts(a_pow.clone(), 0).with_precision(prec_calc).value();
        let den_d = DBig::from_parts(den_i, 0).with_precision(prec_calc).value();
        let term_base = (&a_d / &den_d).with_precision(prec_calc).value();

        let mut n_pow = DBig::from_parts(IBig::from(n_i), 0).with_precision(prec_calc).value();
        for k in 0..6 {
            let term_k = (&term_base / &n_pow).with_precision(prec_calc).value();
            u[k] = &u[k] + &term_k;
            if k < 5 {
                n_pow = (&n_pow * DBig::from_parts(IBig::from(n_i), 0).with_precision(prec_calc).value())
                    .with_precision(prec_calc)
                    .value();
            }
        }

        let abs_t = if term_base < DBig::ZERO { -term_base } else { term_base };
        if abs_t < eps {
            break;
        }
    }

    let [u0, u1, u2, u3, u4, u5] = u;
    Some([
        u0.with_precision(prec).value(),
        u1.with_precision(prec).value(),
        u2.with_precision(prec).value(),
        u3.with_precision(prec).value(),
        u4.with_precision(prec).value(),
        u5.with_precision(prec).value(),
    ])
}

/// Computes cubic pole basis vectors V_k for k = 1..5:
/// V_k = sum_{n=1}^\infty a^n / (n^k * (3n-1)(3n-2) * B_i(n))
pub fn compute_inverse_cubic_pole_bases(
    kernel: InverseKernelType,
    a: i64,
    prec: usize,
    max_steps: usize,
    tables: &BinomialTables,
) -> Option<[DBig; 5]> {
    let growth = kernel.growth_ratio();
    if (a.abs() as f64) > (growth / 2.5) + 1e-6 {
        return None;
    }

    let prec_calc = prec + 30;
    let mut v = [
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
    ];

    let mut a_pow = IBig::ONE;
    let a_ibig = IBig::from(a);
    let eps = DBig::from_parts(IBig::ONE, -240).with_precision(prec_calc).value();

    for n in 1..=max_steps {
        let n_i = n as i64;
        a_pow = &a_pow * &a_ibig;

        let b_val = kernel.eval_at(n, tables);
        let pole_factor = IBig::from((3 * n_i - 1) * (3 * n_i - 2));
        let den_i = &b_val * &pole_factor;

        let a_d = DBig::from_parts(a_pow.clone(), 0).with_precision(prec_calc).value();
        let den_d = DBig::from_parts(den_i, 0).with_precision(prec_calc).value();
        let term_base = (&a_d / &den_d).with_precision(prec_calc).value();

        let mut n_pow = DBig::from_parts(IBig::from(n_i), 0).with_precision(prec_calc).value();
        for k in 0..5 {
            let term_k = (&term_base / &n_pow).with_precision(prec_calc).value();
            v[k] = &v[k] + &term_k;
            if k < 4 {
                n_pow = (&n_pow * DBig::from_parts(IBig::from(n_i), 0).with_precision(prec_calc).value())
                    .with_precision(prec_calc)
                    .value();
            }
        }

        let abs_t = if term_base < DBig::ZERO { -term_base } else { term_base };
        if abs_t < eps {
            break;
        }
    }

    let [v0, v1, v2, v3, v4] = v;
    Some([
        v0.with_precision(prec).value(),
        v1.with_precision(prec).value(),
        v2.with_precision(prec).value(),
        v3.with_precision(prec).value(),
        v4.with_precision(prec).value(),
    ])
}

// ════════════════════════════════════════════════════════════════════════════════
//  5. GERÇEK 7-BİNOMLU (1/pi^3) VE 5-BİNOMLU (1/pi^2) PAY ÇEKİRDEKLERİ (K1..K5)
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NumeratorKernelType {
    K1, // binom(2n, n)^7 (7-binom -> 1/pi^3)
    K2, // binom(2n, n)^5 * binom(3n, n)^2 (7-binom -> 1/pi^3)
    K3, // binom(2n, n)^5 * binom(4n, 2n)^2 (7-binom -> 1/pi^3)
    K4, // binom(2n, n)^4 * binom(4n, 2n) (5-binom -> 1/pi^2)
    K5, // binom(2n, n)^3 * binom(4n, 2n)^2 (5-binom -> 1/pi^2)
}

impl NumeratorKernelType {
    pub const ALL: [NumeratorKernelType; 5] = [
        Self::K1, Self::K2, Self::K3, Self::K4, Self::K5,
    ];

    pub fn code(&self) -> &'static str {
        match self {
            Self::K1 => "K_1",
            Self::K2 => "K_2",
            Self::K3 => "K_3",
            Self::K4 => "K_4",
            Self::K5 => "K_5",
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Self::K1 => "7-Binomlu K_1(n) = binom(2n,n)^7 (1/pi^3)",
            Self::K2 => "7-Binomlu Melez K_2(n) = binom(2n,n)^5 * binom(3n,n)^2 (1/pi^3)",
            Self::K3 => "7-Binomlu Melez K_3(n) = binom(2n,n)^5 * binom(4n,2n)^2 (1/pi^3)",
            Self::K4 => "5-Binomlu Melez K_4(n) = binom(2n,n)^4 * binom(4n,2n) (1/pi^2)",
            Self::K5 => "5-Binomlu Melez K_5(n) = binom(2n,n)^3 * binom(4n,2n)^2 (1/pi^2)",
        }
    }

    pub fn latex(&self) -> &'static str {
        match self {
            Self::K1 => "\\binom{2n}{n}^7",
            Self::K2 => "\\binom{2n}{n}^5 \\binom{3n}{n}^2",
            Self::K3 => "\\binom{2n}{n}^5 \\binom{4n}{2n}^2",
            Self::K4 => "\\binom{2n}{n}^4 \\binom{4n}{2n}",
            Self::K5 => "\\binom{2n}{n}^3 \\binom{4n}{2n}^2",
        }
    }

    pub fn growth_ratio(&self) -> f64 {
        match self {
            Self::K1 => 16384.0,
            Self::K2 => 46656.0,
            Self::K3 => 65536.0,
            Self::K4 => 16384.0,
            Self::K5 => 262144.0,
        }
    }

    pub fn poly_degree(&self) -> usize {
        match self {
            Self::K1 | Self::K2 | Self::K3 => 3, // S0, S1, S2, S3
            Self::K4 | Self::K5 => 2,           // S0, S1, S2
        }
    }

    pub fn eval_at(&self, n: usize, tables: &BinomialTables) -> IBig {
        let c2 = &tables.c2n[n];
        let c3 = &tables.c3n[n];
        let c4 = &tables.c4n[n];

        match self {
            Self::K1 => {
                let sq = c2 * c2;
                let q4 = &sq * &sq;
                (&q4 * &sq) * c2
            }
            Self::K2 => {
                let sq = c2 * c2;
                let q5 = (&sq * &sq) * c2;
                let c3_sq = c3 * c3;
                q5 * c3_sq
            }
            Self::K3 => {
                let sq = c2 * c2;
                let q5 = (&sq * &sq) * c2;
                let c4_sq = c4 * c4;
                q5 * c4_sq
            }
            Self::K4 => {
                let sq = c2 * c2;
                let q4 = &sq * &sq;
                q4 * c4
            }
            Self::K5 => {
                let sq = c2 * c2;
                let q3 = &sq * c2;
                let c4_sq = c4 * c4;
                q3 * c4_sq
            }
        }
    }
}

pub const K1_BASES: [i64; 6] = [
    65536, -65536, 1048576, -1048576, 16777216, -16777216,
];

pub const K2_K3_BASES: [i64; 16] = [
    -65536, 65536, -131072, 131072, -262144, 262144, -524288, 524288,
    -1048576, 1048576, -2097152, 2097152, -4194304, 4194304, -16777216, 16777216,
];

pub const K4_K5_BASES: [i64; 12] = [
    -32768, 32768, -65536, 65536, -131072, 131072, -262144, 262144,
    -524288, 524288, -1048576, 1048576,
];

pub fn get_numerator_bases(kernel: NumeratorKernelType) -> &'static [i64] {
    match kernel {
        NumeratorKernelType::K1 => &K1_BASES,
        NumeratorKernelType::K2 | NumeratorKernelType::K3 => &K2_K3_BASES,
        NumeratorKernelType::K4 | NumeratorKernelType::K5 => &K4_K5_BASES,
    }
}

/// Computes basis sums S_k = sum_{n=0}^\infty n^k * K(n) / a^n up to k=3
pub fn compute_numerator_bases(
    kernel: NumeratorKernelType,
    a: i64,
    prec: usize,
    max_steps: usize,
    tables: &BinomialTables,
) -> Option<Vec<DBig>> {
    let growth = kernel.growth_ratio();
    if (a.abs() as f64) <= growth + 0.1 {
        return None;
    }

    let deg = kernel.poly_degree();
    let num_bases = deg + 1; // S0..S_deg

    let prec_calc = prec + 25;
    let mut s = vec![DBig::ZERO.with_precision(prec_calc).value(); num_bases];

    // n=0 term: K(0) = 1, a^0 = 1, n^0 = 1 => S0 gets 1, others 0
    s[0] = &s[0] + DBig::ONE.with_precision(prec_calc).value();

    let mut a_pow = IBig::ONE;
    let a_ibig = IBig::from(a);
    let eps = DBig::from_parts(IBig::ONE, -235).with_precision(prec_calc).value();

    for n in 1..=max_steps {
        let n_i = n as i64;
        a_pow = &a_pow * &a_ibig;

        let k_val = kernel.eval_at(n, tables);

        let k_d = DBig::from_parts(k_val, 0).with_precision(prec_calc).value();
        let a_d = DBig::from_parts(a_pow.clone(), 0).with_precision(prec_calc).value();
        let term = (&k_d / &a_d).with_precision(prec_calc).value();

        let mut n_mult = DBig::ONE.with_precision(prec_calc).value();
        let n_d = DBig::from_parts(IBig::from(n_i), 0).with_precision(prec_calc).value();

        for k in 0..num_bases {
            s[k] = &s[k] + &(&term * &n_mult);
            n_mult = (&n_mult * &n_d).with_precision(prec_calc).value();
        }

        let abs_t = if term < DBig::ZERO { -term } else { term };
        if abs_t < eps {
            break;
        }
    }

    Some(s.into_iter().map(|v| v.with_precision(prec).value()).collect())
}

pub fn format_descending_polynomial_text(coeffs: &[i64]) -> String {
    // coeffs[0] is n^{d-1}, coeffs[d-1] is constant term
    let d = coeffs.len();
    let mut terms = Vec::new();

    for (i, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let power = d - 1 - i;
        let var = match power {
            0 => "".to_string(),
            1 => "n".to_string(),
            p => format!("n^{p}"),
        };
        terms.push((c, var));
    }

    if terms.is_empty() {
        return "0".to_string();
    }

    let mut res = String::new();
    for (i, (c, var)) in terms.iter().enumerate() {
        let abs_c = c.abs();
        let sign = if *c < 0 {
            if i == 0 { "-" } else { " - " }
        } else if i == 0 {
            ""
        } else {
            " + "
        };

        let term_str = if var.is_empty() {
            format!("{abs_c}")
        } else if abs_c == 1 {
            var.to_string()
        } else {
            format!("{abs_c}*{var}")
        };
        res.push_str(&format!("{sign}{term_str}"));
    }
    res
}

pub fn format_descending_polynomial_latex(coeffs: &[i64]) -> String {
    let d = coeffs.len();
    let mut terms = Vec::new();

    for (i, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let power = d - 1 - i;
        let var = match power {
            0 => "".to_string(),
            1 => "n".to_string(),
            p => format!("n^{p}"),
        };
        terms.push((c, var));
    }

    if terms.is_empty() {
        return "0".to_string();
    }

    let mut res = String::new();
    for (i, (c, var)) in terms.iter().enumerate() {
        let abs_c = c.abs();
        let sign = if *c < 0 {
            if i == 0 { "-" } else { " - " }
        } else if i == 0 {
            ""
        } else {
            " + "
        };

        let term_str = if var.is_empty() {
            format!("{abs_c}")
        } else if abs_c == 1 {
            var.to_string()
        } else {
            format!("{abs_c} {var}")
        };
        res.push_str(&format!("{sign}{term_str}"));
    }
    res
}

// ════════════════════════════════════════════════════════════════════════════════
//  6. AĞIRLIK-5 (ζ(5), π⁶) VE AĞIRLIK-7 (ζ(7)) DERİN KUTUP MOTORU
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DeepPoleKernelType {
    B11, // binom(2n, n)^5 (growth: 1024.0)
    B12, // binom(2n, n)^7 (growth: 16384.0)
    DP3, // binom(2n, n)^3 * binom(3n, n)^2 (growth: 2916.0)
    DP4, // binom(2n, n)^3 * binom(4n, 2n)^2 (growth: 16384.0)
    B7,  // binom(2n, n) * binom(3n, n) * binom(6n, 3n) (growth: 1728.0)
    B15, // binom(2n, n)^2 * binom(3n, n)^3 (growth: 4920.75)
    B16, // binom(2n, n)^4 * binom(3n, n) (growth: 1728.0)
}

impl DeepPoleKernelType {
    pub const ALL: [DeepPoleKernelType; 7] = [
        Self::B11, Self::B12, Self::DP3, Self::DP4, Self::B7, Self::B15, Self::B16,
    ];

    pub fn code(&self) -> &'static str {
        match self {
            Self::B11 => "B_{11}",
            Self::B12 => "B_{12}",
            Self::DP3 => "B_9^{\\text{CY}}",
            Self::DP4 => "B_6^{\\text{Q}}",
            Self::B7  => "B_7",
            Self::B15 => "B_{15}",
            Self::B16 => "B_{16}",
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Self::B11 => "Merkezi 5. Güç B_{11}(n) = binom(2n,n)^5",
            Self::B12 => "Merkezi 7. Güç B_{12}(n) = binom(2n,n)^7",
            Self::DP3 => "5-Binom Calabi-Yau B_9(n) = binom(2n,n)^3 * binom(3n,n)^2",
            Self::DP4 => "5-Binom Kuantik B_6(n) = binom(2n,n)^3 * binom(4n,2n)^2",
            Self::B7  => "Sextic Calabi-Yau B_7(n) = binom(2n,n)*binom(3n,n)*binom(6n,3n)",
            Self::B15 => "5-Binom Kübik B_{15}(n) = binom(2n,n)^2 * binom(3n,n)^3",
            Self::B16 => "5-Binom B_{16}(n) = binom(2n,n)^4 * binom(3n,n)",
        }
    }

    pub fn latex(&self) -> &'static str {
        match self {
            Self::B11 => "\\binom{2n}{n}^5",
            Self::B12 => "\\binom{2n}{n}^7",
            Self::DP3 => "\\binom{2n}{n}^3 \\binom{3n}{n}^2",
            Self::DP4 => "\\binom{2n}{n}^3 \\binom{4n}{2n}^2",
            Self::B7  => "\\binom{2n}{n} \\binom{3n}{n} \\binom{6n}{3n}",
            Self::B15 => "\\binom{2n}{n}^2 \\binom{3n}{n}^3",
            Self::B16 => "\\binom{2n}{n}^4 \\binom{3n}{n}",
        }
    }

    pub fn growth_ratio(&self) -> f64 {
        match self {
            Self::B11 => 1024.0,
            Self::B12 => 16384.0,
            Self::DP3 => 2916.0,
            Self::DP4 => 16384.0,
            Self::B7  => 1728.0,
            Self::B15 => 4920.75,
            Self::B16 => 1728.0,
        }
    }

    pub fn eval_at(&self, n: usize, tables: &BinomialTables) -> IBig {
        let c2 = &tables.c2n[n];
        let c3 = &tables.c3n[n];
        let c4 = &tables.c4n[n];
        let c6 = &tables.c6n[n];

        match self {
            Self::B11 => {
                let sq = c2 * c2;
                (&sq * &sq) * c2
            }
            Self::B12 => {
                let sq = c2 * c2;
                let q4 = &sq * &sq;
                (&q4 * &sq) * c2
            }
            Self::DP3 => {
                let c2_3 = (c2 * c2) * c2;
                let c3_2 = c3 * c3;
                c2_3 * c3_2
            }
            Self::DP4 => {
                let c2_3 = (c2 * c2) * c2;
                let c4_2 = c4 * c4;
                c2_3 * c4_2
            }
            Self::B7 => (c2 * c3) * c6,
            Self::B15 => {
                let c2_2 = c2 * c2;
                let c3_3 = (c3 * c3) * c3;
                c2_2 * c3_3
            }
            Self::B16 => {
                let sq = c2 * c2;
                let c2_4 = &sq * &sq;
                c2_4 * c3
            }
        }
    }
}

pub const DEEP_POLE_BASES: [i64; 16] = [
    -1, 1, -4, 4, -8, 8, -16, 16, -27, 27, -64, 64, -256, 256, -1024, 1024,
];

/// Computes extended deep cyclotomic pole basis vectors U_k^(m) for k = 1..8:
/// U_k^(m) = sum_{n=1}^\infty a^n / (n^k * (2n-1)^m * B(n))
pub fn compute_deep_pole_bases(
    kernel: DeepPoleKernelType,
    a: i64,
    m: usize,
    prec: usize,
    max_steps: usize,
    tables: &BinomialTables,
) -> Option<[DBig; 8]> {
    let growth = kernel.growth_ratio();
    if (a.abs() as f64) > (growth / 2.5) + 1e-6 {
        return None;
    }

    let prec_calc = prec + 30;
    let mut u = [
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
        DBig::ZERO.with_precision(prec_calc).value(),
    ];

    let mut a_pow = IBig::ONE;
    let a_ibig = IBig::from(a);
    let eps = DBig::from_parts(IBig::ONE, -250).with_precision(prec_calc).value();

    for n in 1..=max_steps {
        let n_i = n as i64;
        a_pow = &a_pow * &a_ibig;

        let b_val = kernel.eval_at(n, tables);

        let mut pole_factor = IBig::ONE;
        if m > 0 {
            let odd = IBig::from(2 * n_i - 1);
            for _ in 0..m {
                pole_factor = &pole_factor * &odd;
            }
        }
        let den_i = &b_val * &pole_factor;

        let a_d = DBig::from_parts(a_pow.clone(), 0).with_precision(prec_calc).value();
        let den_d = DBig::from_parts(den_i, 0).with_precision(prec_calc).value();
        let term_base = (&a_d / &den_d).with_precision(prec_calc).value();

        let mut n_pow = DBig::from_parts(IBig::from(n_i), 0).with_precision(prec_calc).value();
        for k in 0..8 {
            let term_k = (&term_base / &n_pow).with_precision(prec_calc).value();
            u[k] = &u[k] + &term_k;
            if k < 7 {
                n_pow = (&n_pow * DBig::from_parts(IBig::from(n_i), 0).with_precision(prec_calc).value())
                    .with_precision(prec_calc)
                    .value();
            }
        }

        let abs_t = if term_base < DBig::ZERO { -term_base } else { term_base };
        if abs_t < eps {
            break;
        }
    }

    let [u0, u1, u2, u3, u4, u5, u6, u7] = u;
    Some([
        u0.with_precision(prec).value(),
        u1.with_precision(prec).value(),
        u2.with_precision(prec).value(),
        u3.with_precision(prec).value(),
        u4.with_precision(prec).value(),
        u5.with_precision(prec).value(),
        u6.with_precision(prec).value(),
        u7.with_precision(prec).value(),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::compute_zeta3;

    #[test]
    fn test_control_point_1_b11_amdeberhan_zeilberger() {
        let prec = 220;
        let tables = compute_binomial_tables(250);
        let u = compute_inverse_cyclotomic_bases(InverseKernelType::B11, -1, 0, prec, 250, &tables)
            .expect("B11 bases must compute");

        // 205*U3 - 160*U4 + 32*U5 + 2*zeta(3) == 0
        let zeta3 = compute_zeta3(prec);
        let c205 = DBig::from_parts(IBig::from(205), 0).with_precision(prec).value();
        let c160 = DBig::from_parts(IBig::from(160), 0).with_precision(prec).value();
        let c32 = DBig::from_parts(IBig::from(32), 0).with_precision(prec).value();
        let c2 = DBig::from_parts(IBig::from(2), 0).with_precision(prec).value();

        let term1 = (&c205 * &u[2]).with_precision(prec).value();
        let term2 = (&c160 * &u[3]).with_precision(prec).value();
        let term3 = (&c32 * &u[4]).with_precision(prec).value();
        let rhs = (&c2 * &zeta3).with_precision(prec).value();

        let sum = (&(&term1 - &term2) + &term3).with_precision(prec).value();
        let diff = if sum >= -rhs.clone() {
            &sum + &rhs
        } else {
            -(&sum + &rhs)
        };

        let thresh = DBig::from_parts(IBig::ONE, -120).with_precision(prec).value();
        assert!(
            diff < thresh,
            "Control Point 1 (B11) failed: diff = {:?}",
            diff
        );
    }

    #[test]
    fn test_control_point_2_b1_deuring_zeta3() {
        let prec = 220;
        let tables = compute_binomial_tables(250);
        let u = compute_inverse_cyclotomic_bases(InverseKernelType::B1, -1, 2, prec, 250, &tables)
            .expect("B1 bases must compute");

        // 56*U1 - 32*U2 + 5*U3 + 4*zeta(3) == 0
        let zeta3 = compute_zeta3(prec);
        let c56 = DBig::from_parts(IBig::from(56), 0).with_precision(prec).value();
        let c32 = DBig::from_parts(IBig::from(32), 0).with_precision(prec).value();
        let c5 = DBig::from_parts(IBig::from(5), 0).with_precision(prec).value();
        let c4 = DBig::from_parts(IBig::from(4), 0).with_precision(prec).value();

        let term1 = (&c56 * &u[0]).with_precision(prec).value();
        let term2 = (&c32 * &u[1]).with_precision(prec).value();
        let term3 = (&c5 * &u[2]).with_precision(prec).value();
        let rhs = (&c4 * &zeta3).with_precision(prec).value();

        let sum = (&(&term1 - &term2) + &term3).with_precision(prec).value();
        let diff = if sum >= -rhs.clone() {
            &sum + &rhs
        } else {
            -(&sum + &rhs)
        };

        let thresh = DBig::from_parts(IBig::ONE, -120).with_precision(prec).value();
        assert!(
            diff < thresh,
            "Control Point 2 (B1) failed: diff = {:?}",
            diff
        );
    }
}

