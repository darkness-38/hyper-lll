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
