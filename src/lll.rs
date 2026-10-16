use dashu::float::DBig;
use dashu::integer::IBig;
use dashu::rational::RBig;

pub const SCALE_EXPONENT: isize = 120;
pub const VERIFY_EXPONENT: isize = -100;
pub const MAX_NORM: f64 = 30000.0;

pub const GUILLERA_SCALE_EXPONENT: isize = 130;
pub const GUILLERA_VERIFY_EXPONENT: isize = -110;
pub const GUILLERA_MAX_NORM: f64 = 80000.0;

pub const FRONTIER_SCALE_EXPONENT: isize = 140;
pub const FRONTIER_VERIFY_EXPONENT: isize = -115;
pub const FRONTIER_MAX_NORM: f64 = 120000.0;

pub const CALABI_YAU_SCALE_EXPONENT: isize = 145;
pub const CALABI_YAU_VERIFY_EXPONENT: isize = -120;
pub const CALABI_YAU_MAX_NORM: f64 = 500000.0;

#[derive(Debug, Clone)]
pub struct IntegerRelation {
    pub m: Vec<i64>,
    #[allow(dead_code)]
    pub norm: f64,
    #[allow(dead_code)]
    pub residual_dbig: DBig,
    #[allow(dead_code)]
    pub residual_f64: f64,
}

pub fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a.abs()
}

pub fn vec_gcd(v: &[i64]) -> i64 {
    let mut g = 0;
    for &x in v {
        g = gcd(g, x);
        if g == 1 {
            break;
        }
    }
    g
}

fn dot_product_rbig(v1: &[RBig], v2: &[RBig]) -> RBig {
    let mut sum = RBig::ZERO;
    for (a, b) in v1.iter().zip(v2.iter()) {
        sum += a * b;
    }
    sum
}

fn round_rational(r: &RBig) -> IBig {
    let n = r.numerator();
    let d = IBig::from(r.denominator().clone());
    let two_n_plus_d = n * 2 + &d;
    let two_d = &d * 2;
    if two_n_plus_d >= IBig::ZERO {
        two_n_plus_d / two_d
    } else {
        (two_n_plus_d - &two_d + 1) / two_d
    }
}

/// Arbitrary dimension LLL lattice reduction.
/// basis is an n x m matrix of IBig (n vectors of dimension m = n + 1).
pub fn lll_reduce_high_prec(basis: &mut Vec<Vec<IBig>>, delta_num: i64, delta_den: i64) {
    let n = basis.len();
    if n <= 1 {
        return;
    }
    let m = basis[0].len();

    let mut b_star: Vec<Vec<RBig>> = vec![vec![RBig::ZERO; m]; n];
    let mut mu: Vec<Vec<RBig>> = vec![vec![RBig::ZERO; n]; n];
    let mut b_star_sq: Vec<RBig> = vec![RBig::ZERO; n];

    // Initial Gram-Schmidt orthogonalization
    for i in 0..n {
        let mut v: Vec<RBig> = basis[i].iter().map(|x| RBig::from(x.clone())).collect();
        for j in 0..i {
            let dot = dot_product_rbig(&v, &b_star[j]);
            let u = dot / &b_star_sq[j];
            mu[i][j] = u.clone();
            for col in 0..m {
                let term = &u * &b_star[j][col];
                v[col] -= term;
            }
        }
        b_star_sq[i] = dot_product_rbig(&v, &v);
        b_star[i] = v;
    }

    let delta_ratio = RBig::from_parts(
        IBig::from(delta_num),
        dashu::integer::UBig::from(delta_den as u64),
    );

    let mut k = 1;
    while k < n {
        let q = round_rational(&mu[k][k - 1]);
        if !q.is_zero() {
            let q_r = RBig::from(q.clone());
            for col in 0..m {
                let sub = &q * &basis[k - 1][col];
                basis[k][col] -= sub;
            }
            for l in 0..(k - 1) {
                let term = &q_r * &mu[k - 1][l];
                mu[k][l] -= term;
            }
            mu[k][k - 1] -= &q_r;
        }

        // Lovász condition: ||b_k*||^2 >= (delta - mu_{k, k-1}^2) * ||b_{k-1}*||^2
        let lhs = &b_star_sq[k];
        let rhs = (&delta_ratio - &mu[k][k - 1].clone().pow(2)) * &b_star_sq[k - 1];

        if lhs >= &rhs {
            for j in (0..k.saturating_sub(1)).rev() {
                let q = round_rational(&mu[k][j]);
                if !q.is_zero() {
                    let q_r = RBig::from(q.clone());
                    for col in 0..m {
                        let sub = &q * &basis[j][col];
                        basis[k][col] -= sub;
                    }
                    for l in 0..j {
                        let term = &q_r * &mu[j][l];
                        mu[k][l] -= term;
                    }
                    mu[k][j] -= &q_r;
                }
            }
            k += 1;
        } else {
            // Swap basis[k] and basis[k-1] with analytical Gram-Schmidt update
            basis.swap(k, k - 1);
            let u = mu[k][k - 1].clone();
            let b_km1_sq = b_star_sq[k - 1].clone();
            let b_k_sq = b_star_sq[k].clone();
            let b_prime = &b_k_sq + &u.clone().pow(2) * &b_km1_sq;

            let new_u_k_km1 = (&u * &b_km1_sq) / &b_prime;
            let new_b_k_sq = (&b_km1_sq * &b_k_sq) / &b_prime;
            let new_b_km1_sq = b_prime.clone();

            let mut new_b_star_km1 = vec![RBig::ZERO; m];
            let mut new_b_star_k = vec![RBig::ZERO; m];
            for col in 0..m {
                new_b_star_km1[col] = &b_star[k][col] + &u * &b_star[k - 1][col];
                new_b_star_k[col] = &b_star[k - 1][col] - &new_u_k_km1 * &new_b_star_km1[col];
            }
            b_star[k - 1] = new_b_star_km1;
            b_star[k] = new_b_star_k;
            b_star_sq[k - 1] = new_b_km1_sq;
            b_star_sq[k] = new_b_k_sq;

            for j in 0..(k - 1) {
                let tmp = mu[k - 1][j].clone();
                mu[k - 1][j] = mu[k][j].clone();
                mu[k][j] = tmp;
            }
            mu[k][k - 1] = new_u_k_km1;

            for i in (k + 1)..n {
                let mu_i_km1 = mu[i][k - 1].clone();
                let mu_i_k = mu[i][k].clone();
                mu[i][k - 1] = (&mu_i_k * &b_k_sq + &mu_i_km1 * &u * &b_km1_sq) / &b_prime;
                mu[i][k] = mu_i_km1 - &u * &mu_i_k;
            }

            k = k.saturating_sub(1).max(1);
        }
    }
}

/// Solves all integer relations for vector `x` with custom LLL parameters.
pub fn find_all_integer_relations_custom(
    x: &[DBig],
    scale_exponent: isize,
    verify_exponent: isize,
    max_norm: f64,
) -> Vec<IntegerRelation> {
    let d = x.len();
    let prec = 200;
    let scale = DBig::from_parts(IBig::ONE, scale_exponent).with_precision(prec).value();
    let half = DBig::from_parts(IBig::from(5), -1).with_precision(prec).value();

    let mut basis: Vec<Vec<IBig>> = Vec::with_capacity(d);
    for i in 0..d {
        let mut row = vec![IBig::ZERO; d + 1];
        row[i] = IBig::ONE;
        let scaled = &x[i] * &scale;
        let rounded = if scaled >= DBig::ZERO {
            &scaled + &half
        } else {
            &scaled - &half
        };
        row[d] = rounded.to_int().value();
        basis.push(row);
    }

    lll_reduce_high_prec(&mut basis, 99, 100);

    let threshold = DBig::from_parts(IBig::ONE, verify_exponent).with_precision(prec).value();
    let mut relations = Vec::new();
    let mut seen_signatures = std::collections::HashSet::new();

    for row in &basis {
        let mut m_vec: Vec<i64> = Vec::with_capacity(d);
        let mut ok = true;
        for i in 0..d {
            match i64::try_from(&row[i]) {
                Ok(v) => m_vec.push(v),
                Err(_) => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            continue;
        }

        // Must not be all zeros
        if m_vec.iter().all(|&v| v == 0) {
            continue;
        }

        // GCD reduction
        let g = vec_gcd(&m_vec);
        if g > 1 {
            for val in &mut m_vec {
                *val /= g;
            }
        }

        // Normalize first non-zero to positive
        if let Some(&first_non_zero) = m_vec.iter().find(|&&v| v != 0) {
            if first_non_zero < 0 {
                for val in &mut m_vec {
                    *val = -*val;
                }
            }
        }

        let sig = format!("{:?}", m_vec);
        if seen_signatures.contains(&sig) {
            continue;
        }

        // Euclidean norm check: norm <= max_norm
        let norm_sq: f64 = m_vec.iter().map(|&v| (v as f64) * (v as f64)).sum();
        let norm = norm_sq.sqrt();
        if norm > max_norm {
            continue;
        }

        // Independent high-precision verification: |sum m_i x_i| < threshold
        let mut residual = DBig::ZERO.with_precision(prec).value();
        for i in 0..d {
            let m_d = DBig::from_parts(IBig::from(m_vec[i]), 0).with_precision(prec).value();
            residual = &residual + &(&m_d * &x[i]);
        }
        let abs_res = if residual < DBig::ZERO { -residual.clone() } else { residual.clone() };

        if abs_res < threshold {
            let res_f64 = abs_res.to_f64().value();
            seen_signatures.insert(sig);
            relations.push(IntegerRelation {
                m: m_vec,
                norm,
                residual_dbig: abs_res,
                residual_f64: res_f64,
            });
        }
    }

    relations
}

/// Solves all integer relations for vector `x` with LLL at 10^120 scale and independent verification at < 10^-100.
pub fn find_all_integer_relations_generic(x: &[DBig]) -> Vec<IntegerRelation> {
    find_all_integer_relations_custom(x, SCALE_EXPONENT, VERIFY_EXPONENT, MAX_NORM)
}

/// Solves all integer relations for Ramanujan-Guillera vector `x` with LLL at 10^130 scale and verification at < 10^-110 (max norm 80000).
pub fn find_all_integer_relations_guillera(x: &[DBig]) -> Vec<IntegerRelation> {
    find_all_integer_relations_custom(
        x,
        GUILLERA_SCALE_EXPONENT,
        GUILLERA_VERIFY_EXPONENT,
        GUILLERA_MAX_NORM,
    )
}

/// Solves all integer relations for Frontier vector `x` with LLL at 10^140 scale and verification at < 10^-115 (max norm 120000).
pub fn find_all_integer_relations_frontier(x: &[DBig]) -> Vec<IntegerRelation> {
    find_all_integer_relations_custom(
        x,
        FRONTIER_SCALE_EXPONENT,
        FRONTIER_VERIFY_EXPONENT,
        FRONTIER_MAX_NORM,
    )
}

/// Solves all integer relations for Calabi-Yau / Super-binomial vector `x` with LLL at 10^145 scale and verification at < 10^-120 (max norm 500000).
pub fn find_all_integer_relations_calabi_yau(x: &[DBig]) -> Vec<IntegerRelation> {
    find_all_integer_relations_custom(
        x,
        CALABI_YAU_SCALE_EXPONENT,
        CALABI_YAU_VERIFY_EXPONENT,
        CALABI_YAU_MAX_NORM,
    )
}

#[allow(dead_code)]
pub fn find_integer_relation_generic(x: &[DBig]) -> Option<IntegerRelation> {
    find_all_integer_relations_generic(x).into_iter().next()
}

/// Fast stack-allocated 4D LLL pre-filter (~2 microseconds)
#[allow(dead_code)]
pub fn quick_check_relation_4d(x: &[f64; 4]) -> Option<[i64; 4]> {
    let scale: f64 = 1e9;
    let mut basis = [[0i64; 5]; 4];
    for i in 0..4 {
        basis[i][i] = 1;
        basis[i][4] = (x[i] * scale).round() as i64;
    }

    let mut b_star = [[0.0f64; 5]; 4];
    let mut mu = [[0.0f64; 4]; 4];
    let mut b_star_sq = [0.0f64; 4];

    for i in 0..4 {
        let mut v = [0.0f64; 5];
        for col in 0..5 {
            v[col] = basis[i][col] as f64;
        }
        for j in 0..i {
            let mut dot = 0.0;
            for col in 0..5 {
                dot += v[col] * b_star[j][col];
            }
            let u = dot / b_star_sq[j];
            mu[i][j] = u;
            for col in 0..5 {
                v[col] -= u * b_star[j][col];
            }
        }
        let mut sq = 0.0;
        for col in 0..5 {
            sq += v[col] * v[col];
        }
        b_star_sq[i] = sq;
        b_star[i] = v;
    }

    let mut k = 1;
    let mut iters = 0;
    while k < 4 && iters < 200 {
        iters += 1;
        let q = mu[k][k - 1].round() as i64;
        if q != 0 {
            for col in 0..5 {
                basis[k][col] -= q * basis[k - 1][col];
            }
            let q_f = q as f64;
            for l in 0..(k - 1) {
                mu[k][l] -= q_f * mu[k - 1][l];
            }
            mu[k][k - 1] -= q_f;
        }

        let lhs = b_star_sq[k];
        let rhs = (0.99 - mu[k][k - 1] * mu[k][k - 1]) * b_star_sq[k - 1];

        if lhs >= rhs {
            for j in (0..k.saturating_sub(1)).rev() {
                let q = mu[k][j].round() as i64;
                if q != 0 {
                    let q_f = q as f64;
                    for col in 0..5 {
                        basis[k][col] -= q * basis[j][col];
                    }
                    for l in 0..j {
                        mu[k][l] -= q_f * mu[j][l];
                    }
                    mu[k][j] -= q_f;
                }
            }
            k += 1;
        } else {
            basis.swap(k, k - 1);
            let u = mu[k][k - 1];
            let b_km1_sq = b_star_sq[k - 1];
            let b_k_sq = b_star_sq[k];
            let b_prime = b_k_sq + u * u * b_km1_sq;

            if b_prime < 1e-18 {
                return None;
            }

            let new_u_k_km1 = (u * b_km1_sq) / b_prime;
            let new_b_k_sq = (b_km1_sq * b_k_sq) / b_prime;
            let new_b_km1_sq = b_prime;

            let mut new_b_star_km1 = [0.0; 5];
            let mut new_b_star_k = [0.0; 5];
            for col in 0..5 {
                new_b_star_km1[col] = b_star[k][col] + u * b_star[k - 1][col];
                new_b_star_k[col] = b_star[k - 1][col] - new_u_k_km1 * new_b_star_km1[col];
            }
            b_star[k - 1] = new_b_star_km1;
            b_star[k] = new_b_star_k;
            b_star_sq[k - 1] = new_b_km1_sq;
            b_star_sq[k] = new_b_k_sq;

            for j in 0..(k - 1) {
                let tmp = mu[k - 1][j];
                mu[k - 1][j] = mu[k][j];
                mu[k][j] = tmp;
            }
            mu[k][k - 1] = new_u_k_km1;

            for i in (k + 1)..4 {
                let mu_i_km1 = mu[i][k - 1];
                let mu_i_k = mu[i][k];
                mu[i][k - 1] = (mu_i_k * b_k_sq + mu_i_km1 * u * b_km1_sq) / b_prime;
                mu[i][k] = mu_i_km1 - u * mu_i_k;
            }

            k = k.saturating_sub(1).max(1);
        }
    }

    for row in &basis {
        let (mut m1, mut m2, mut m3, mut m4) = (row[0], row[1], row[2], row[3]);
        if m1 == 0 && m2 == 0 && m3 == 0 && m4 == 0 {
            continue;
        }

        let g = gcd(gcd(m1.abs(), m2.abs()), gcd(m3.abs(), m4.abs()));
        if g > 1 {
            m1 /= g;
            m2 /= g;
            m3 /= g;
            m4 /= g;
        }
        if m1 < 0 {
            m1 = -m1;
            m2 = -m2;
            m3 = -m3;
            m4 = -m4;
        }

        let norm_sq = (m1 * m1 + m2 * m2 + m3 * m3 + m4 * m4) as f64;
        if norm_sq.sqrt() <= MAX_NORM && row[4].abs() < 5000 {
            let eval_f64 = (m1 as f64) * x[0] + (m2 as f64) * x[1] + (m3 as f64) * x[2] + (m4 as f64) * x[3];
            if eval_f64.abs() < 1e-6 {
                return Some([m1, m2, m3, m4]);
            }
        }
    }

    None
}

/// Fast stack-allocated 5D LLL pre-filter (~3 microseconds)
#[allow(dead_code)]
pub fn quick_check_relation_5d(x: &[f64; 5]) -> Option<[i64; 5]> {
    let scale: f64 = 1e9;
    let mut basis = [[0i64; 6]; 5];
    for i in 0..5 {
        basis[i][i] = 1;
        basis[i][5] = (x[i] * scale).round() as i64;
    }

    let mut b_star = [[0.0f64; 6]; 5];
    let mut mu = [[0.0f64; 5]; 5];
    let mut b_star_sq = [0.0f64; 5];

    for i in 0..5 {
        let mut v = [0.0f64; 6];
        for col in 0..6 {
            v[col] = basis[i][col] as f64;
        }
        for j in 0..i {
            let mut dot = 0.0;
            for col in 0..6 {
                dot += v[col] * b_star[j][col];
            }
            let u = dot / b_star_sq[j];
            mu[i][j] = u;
            for col in 0..6 {
                v[col] -= u * b_star[j][col];
            }
        }
        let mut sq = 0.0;
        for col in 0..6 {
            sq += v[col] * v[col];
        }
        b_star_sq[i] = sq;
        b_star[i] = v;
    }

    let mut k = 1;
    let mut iters = 0;
    while k < 5 && iters < 250 {
        iters += 1;
        let q = mu[k][k - 1].round() as i64;
        if q != 0 {
            for col in 0..6 {
                basis[k][col] -= q * basis[k - 1][col];
            }
            let q_f = q as f64;
            for l in 0..(k - 1) {
                mu[k][l] -= q_f * mu[k - 1][l];
            }
            mu[k][k - 1] -= q_f;
        }

        let lhs = b_star_sq[k];
        let rhs = (0.99 - mu[k][k - 1] * mu[k][k - 1]) * b_star_sq[k - 1];

        if lhs >= rhs {
            for j in (0..k.saturating_sub(1)).rev() {
                let q = mu[k][j].round() as i64;
                if q != 0 {
                    let q_f = q as f64;
                    for col in 0..6 {
                        basis[k][col] -= q * basis[j][col];
                    }
                    for l in 0..j {
                        mu[k][l] -= q_f * mu[j][l];
                    }
                    mu[k][j] -= q_f;
                }
            }
            k += 1;
        } else {
            basis.swap(k, k - 1);
            let u = mu[k][k - 1];
            let b_km1_sq = b_star_sq[k - 1];
            let b_k_sq = b_star_sq[k];
            let b_prime = b_k_sq + u * u * b_km1_sq;

            if b_prime < 1e-18 {
                return None;
            }

            let new_u_k_km1 = (u * b_km1_sq) / b_prime;
            let new_b_k_sq = (b_km1_sq * b_k_sq) / b_prime;
            let new_b_km1_sq = b_prime;

            let mut new_b_star_km1 = [0.0; 6];
            let mut new_b_star_k = [0.0; 6];
            for col in 0..6 {
                new_b_star_km1[col] = b_star[k][col] + u * b_star[k - 1][col];
                new_b_star_k[col] = b_star[k - 1][col] - new_u_k_km1 * new_b_star_km1[col];
            }
            b_star[k - 1] = new_b_star_km1;
            b_star[k] = new_b_star_k;
            b_star_sq[k - 1] = new_b_km1_sq;
            b_star_sq[k] = new_b_k_sq;

            for j in 0..(k - 1) {
                let tmp = mu[k - 1][j];
                mu[k - 1][j] = mu[k][j];
                mu[k][j] = tmp;
            }
            mu[k][k - 1] = new_u_k_km1;

            for i in (k + 1)..5 {
                let mu_i_km1 = mu[i][k - 1];
                let mu_i_k = mu[i][k];
                mu[i][k - 1] = (mu_i_k * b_k_sq + mu_i_km1 * u * b_km1_sq) / b_prime;
                mu[i][k] = mu_i_km1 - u * mu_i_k;
            }

            k = k.saturating_sub(1).max(1);
        }
    }

    for row in &basis {
        let (mut m1, mut m2, mut m3, mut m4, mut m5) = (row[0], row[1], row[2], row[3], row[4]);
        if m1 == 0 && m2 == 0 && m3 == 0 && m4 == 0 && m5 == 0 {
            continue;
        }

        let g = vec_gcd(&[m1, m2, m3, m4, m5]);
        if g > 1 {
            m1 /= g;
            m2 /= g;
            m3 /= g;
            m4 /= g;
            m5 /= g;
        }
        if m1 < 0 || (m1 == 0 && m2 < 0) {
            m1 = -m1;
            m2 = -m2;
            m3 = -m3;
            m4 = -m4;
            m5 = -m5;
        }

        let norm_sq = (m1 * m1 + m2 * m2 + m3 * m3 + m4 * m4 + m5 * m5) as f64;
        if norm_sq.sqrt() <= MAX_NORM && row[5].abs() < 5000 {
            let eval_f64 = (m1 as f64) * x[0]
                + (m2 as f64) * x[1]
                + (m3 as f64) * x[2]
                + (m4 as f64) * x[3]
                + (m5 as f64) * x[4];
            if eval_f64.abs() < 1e-6 {
                return Some([m1, m2, m3, m4, m5]);
            }
        }
    }

    None
}

/// Fast stack-allocated 6D LLL pre-filter (~4 microseconds)
#[allow(dead_code)]
pub fn quick_check_relation_6d(x: &[f64; 6]) -> Option<[i64; 6]> {
    let scale: f64 = 1e9;
    let mut basis = [[0i64; 7]; 6];
    for i in 0..6 {
        basis[i][i] = 1;
        basis[i][6] = (x[i] * scale).round() as i64;
    }

    let mut b_star = [[0.0f64; 7]; 6];
    let mut mu = [[0.0f64; 6]; 6];
    let mut b_star_sq = [0.0f64; 6];

    for i in 0..6 {
        let mut v = [0.0f64; 7];
        for col in 0..7 {
            v[col] = basis[i][col] as f64;
        }
        for j in 0..i {
            let mut dot = 0.0;
            for col in 0..7 {
                dot += v[col] * b_star[j][col];
            }
            let u = dot / b_star_sq[j];
            mu[i][j] = u;
            for col in 0..7 {
                v[col] -= u * b_star[j][col];
            }
        }
        let mut sq = 0.0;
        for col in 0..7 {
            sq += v[col] * v[col];
        }
        b_star_sq[i] = sq;
        b_star[i] = v;
    }

    let mut k = 1;
    let mut iters = 0;
    while k < 6 && iters < 300 {
        iters += 1;
        let q = mu[k][k - 1].round() as i64;
        if q != 0 {
            for col in 0..7 {
                basis[k][col] -= q * basis[k - 1][col];
            }
            let q_f = q as f64;
            for l in 0..(k - 1) {
                mu[k][l] -= q_f * mu[k - 1][l];
            }
            mu[k][k - 1] -= q_f;
        }

        let lhs = b_star_sq[k];
        let rhs = (0.99 - mu[k][k - 1] * mu[k][k - 1]) * b_star_sq[k - 1];

        if lhs >= rhs {
            for j in (0..k.saturating_sub(1)).rev() {
                let q = mu[k][j].round() as i64;
                if q != 0 {
                    let q_f = q as f64;
                    for col in 0..7 {
                        basis[k][col] -= q * basis[j][col];
                    }
                    for l in 0..j {
                        mu[k][l] -= q_f * mu[j][l];
                    }
                    mu[k][j] -= q_f;
                }
            }
            k += 1;
        } else {
            basis.swap(k, k - 1);
            let u = mu[k][k - 1];
            let b_km1_sq = b_star_sq[k - 1];
            let b_k_sq = b_star_sq[k];
            let b_prime = b_k_sq + u * u * b_km1_sq;

            if b_prime < 1e-18 {
                return None;
            }

            let new_u_k_km1 = (u * b_km1_sq) / b_prime;
            let new_b_k_sq = (b_km1_sq * b_k_sq) / b_prime;
            let new_b_km1_sq = b_prime;

            let mut new_b_star_km1 = [0.0; 7];
            let mut new_b_star_k = [0.0; 7];
            for col in 0..7 {
                new_b_star_km1[col] = b_star[k][col] + u * b_star[k - 1][col];
                new_b_star_k[col] = b_star[k - 1][col] - new_u_k_km1 * new_b_star_km1[col];
            }
            b_star[k - 1] = new_b_star_km1;
            b_star[k] = new_b_star_k;
            b_star_sq[k - 1] = new_b_km1_sq;
            b_star_sq[k] = new_b_k_sq;

            for j in 0..(k - 1) {
                let tmp = mu[k - 1][j];
                mu[k - 1][j] = mu[k][j];
                mu[k][j] = tmp;
            }
            mu[k][k - 1] = new_u_k_km1;

            for i in (k + 1)..6 {
                let mu_i_km1 = mu[i][k - 1];
                let mu_i_k = mu[i][k];
                mu[i][k - 1] = (mu_i_k * b_k_sq + mu_i_km1 * u * b_km1_sq) / b_prime;
                mu[i][k] = mu_i_km1 - u * mu_i_k;
            }

            k = k.saturating_sub(1).max(1);
        }
    }

    for row in &basis {
        let (mut m1, mut m2, mut m3, mut m4, mut m5, mut m6) = (row[0], row[1], row[2], row[3], row[4], row[5]);
        if m1 == 0 && m2 == 0 && m3 == 0 && m4 == 0 && m5 == 0 && m6 == 0 {
            continue;
        }

        let g = vec_gcd(&[m1, m2, m3, m4, m5, m6]);
        if g > 1 {
            m1 /= g;
            m2 /= g;
            m3 /= g;
            m4 /= g;
            m5 /= g;
            m6 /= g;
        }
        if m1 < 0 || (m1 == 0 && m2 < 0) {
            m1 = -m1;
            m2 = -m2;
            m3 = -m3;
            m4 = -m4;
            m5 = -m5;
            m6 = -m6;
        }

        let norm_sq = (m1 * m1 + m2 * m2 + m3 * m3 + m4 * m4 + m5 * m5 + m6 * m6) as f64;
        if norm_sq.sqrt() <= MAX_NORM && row[6].abs() < 5000 {
            let eval_f64 = (m1 as f64) * x[0]
                + (m2 as f64) * x[1]
                + (m3 as f64) * x[2]
                + (m4 as f64) * x[3]
                + (m5 as f64) * x[4]
                + (m6 as f64) * x[5];
            if eval_f64.abs() < 1e-6 {
                return Some([m1, m2, m3, m4, m5, m6]);
            }
        }
    }

    None
}

/// Fast stack-allocated 7D LLL pre-filter (~5 microseconds)
#[allow(dead_code)]
pub fn quick_check_relation_7d(x: &[f64; 7]) -> Option<[i64; 7]> {
    let scale: f64 = 1e9;
    let mut basis = [[0i64; 8]; 7];
    for i in 0..7 {
        basis[i][i] = 1;
        basis[i][7] = (x[i] * scale).round() as i64;
    }

    let mut b_star = [[0.0f64; 8]; 7];
    let mut mu = [[0.0f64; 7]; 7];
    let mut b_star_sq = [0.0f64; 7];

    for i in 0..7 {
        let mut v = [0.0f64; 8];
        for col in 0..8 {
            v[col] = basis[i][col] as f64;
        }
        for j in 0..i {
            let mut dot = 0.0;
            for col in 0..8 {
                dot += v[col] * b_star[j][col];
            }
            let u = dot / b_star_sq[j];
            mu[i][j] = u;
            for col in 0..8 {
                v[col] -= u * b_star[j][col];
            }
        }
        let mut sq = 0.0;
        for col in 0..8 {
            sq += v[col] * v[col];
        }
        b_star_sq[i] = sq;
        b_star[i] = v;
    }

    let mut k = 1;
    let mut iters = 0;
    while k < 7 && iters < 250 {
        iters += 1;
        let q = mu[k][k - 1].round() as i64;
        if q != 0 {
            for col in 0..8 {
                basis[k][col] -= q * basis[k - 1][col];
            }
            let q_f = q as f64;
            for l in 0..(k - 1) {
                mu[k][l] -= q_f * mu[k - 1][l];
            }
            mu[k][k - 1] -= q_f;
        }

        let lhs = b_star_sq[k];
        let rhs = (0.99 - mu[k][k - 1] * mu[k][k - 1]) * b_star_sq[k - 1];

        if lhs >= rhs {
            for j in (0..k.saturating_sub(1)).rev() {
                let q = mu[k][j].round() as i64;
                if q != 0 {
                    let q_f = q as f64;
                    for col in 0..8 {
                        basis[k][col] -= q * basis[j][col];
                    }
                    for l in 0..j {
                        mu[k][l] -= q_f * mu[j][l];
                    }
                    mu[k][j] -= q_f;
                }
            }
            k += 1;
        } else {
            basis.swap(k, k - 1);
            let u = mu[k][k - 1];
            let b_km1_sq = b_star_sq[k - 1];
            let b_k_sq = b_star_sq[k];
            let b_prime = b_k_sq + u * u * b_km1_sq;

            if b_prime < 1e-18 {
                return None;
            }

            let new_u_k_km1 = (u * b_km1_sq) / b_prime;
            let new_b_k_sq = (b_km1_sq * b_k_sq) / b_prime;
            let new_b_km1_sq = b_prime;

            let mut new_b_star_km1 = [0.0; 8];
            let mut new_b_star_k = [0.0; 8];
            for col in 0..8 {
                new_b_star_km1[col] = b_star[k][col] + u * b_star[k - 1][col];
                new_b_star_k[col] = b_star[k - 1][col] - new_u_k_km1 * new_b_star_km1[col];
            }
            b_star[k - 1] = new_b_star_km1;
            b_star[k] = new_b_star_k;
            b_star_sq[k - 1] = new_b_km1_sq;
            b_star_sq[k] = new_b_k_sq;

            for j in 0..(k - 1) {
                let tmp = mu[k - 1][j];
                mu[k - 1][j] = mu[k][j];
                mu[k][j] = tmp;
            }
            mu[k][k - 1] = new_u_k_km1;

            for i in (k + 1)..7 {
                let mu_i_km1 = mu[i][k - 1];
                let mu_i_k = mu[i][k];
                mu[i][k - 1] = (mu_i_k * b_k_sq + mu_i_km1 * u * b_km1_sq) / b_prime;
                mu[i][k] = mu_i_km1 - u * mu_i_k;
            }

            k = k.saturating_sub(1).max(1);
        }
    }

    for row in &basis {
        let mut m = [row[0], row[1], row[2], row[3], row[4], row[5], row[6]];
        if m.iter().all(|&v| v == 0) {
            continue;
        }

        let g = vec_gcd(&m);
        if g > 1 {
            for v in &mut m {
                *v /= g;
            }
        }
        if let Some(&first_non_zero) = m.iter().find(|&&v| v != 0) {
            if first_non_zero < 0 {
                for v in &mut m {
                    *v = -*v;
                }
            }
        }

        let norm_sq: f64 = m.iter().map(|&v| (v as f64) * (v as f64)).sum();
        if norm_sq.sqrt() <= MAX_NORM && row[7].abs() < 5000 {
            let mut eval_f64 = 0.0;
            for i in 0..7 {
                eval_f64 += (m[i] as f64) * x[i];
            }
            if eval_f64.abs() < 1e-6 {
                return Some(m);
            }
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    const PI_STR: &str = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196";
    const G_STR: &str = "0.915965594177219015054603514932384110774149374281672134266498119621763019776254769479356512926115106248574422619196199579035898803325859059431594737481158406912972923635417054297";
    const LN2_STR: &str = "0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202357581305570326707516350759619307275087260465830172";

    #[test]
    fn test_known_artificial_relation_4d() {
        let prec = 160;
        let s1 = DBig::from_str(PI_STR).unwrap().with_precision(prec).value();
        let s2 = DBig::from_str(G_STR).unwrap().with_precision(prec).value();

        let two = DBig::from_parts(IBig::from(2), 0).with_precision(prec).value();
        let five = DBig::from_parts(IBig::from(5), 0).with_precision(prec).value();
        let seven = DBig::from_parts(IBig::from(7), 0).with_precision(prec).value();
        let three = DBig::from_parts(IBig::from(3), 0).with_precision(prec).value();

        // 3*K - 2*S1 + 5*S2 - 7 = 0
        let num = &(&two * &s1) - &(&five * &s2) + &seven;
        let k = &num / &three;
        let one = DBig::ONE.with_precision(prec).value();

        let x = vec![k, s1, s2, one];
        let rel = find_integer_relation_generic(&x).expect("4D relation must be found");

        assert_eq!(rel.m, vec![3, -2, 5, -7]);
        assert!(rel.norm <= MAX_NORM);
        let thresh = DBig::from_parts(IBig::ONE, VERIFY_EXPONENT).with_precision(prec).value();
        assert!(rel.residual_dbig < thresh);
    }

    #[test]
    fn test_known_artificial_relation_5d() {
        let prec = 160;
        let s1 = DBig::from_str(PI_STR).unwrap().with_precision(prec).value();
        let s2 = DBig::from_str(LN2_STR).unwrap().with_precision(prec).value();

        // K = (2*S1 - 3*S2 + 4) / (5*S1 + 1)
        // => 5*(K*S1) + 1*K - 2*S1 + 3*S2 - 4 = 0
        let five = DBig::from_parts(IBig::from(5), 0).with_precision(prec).value();
        let two = DBig::from_parts(IBig::from(2), 0).with_precision(prec).value();
        let three = DBig::from_parts(IBig::from(3), 0).with_precision(prec).value();
        let four = DBig::from_parts(IBig::from(4), 0).with_precision(prec).value();
        let one = DBig::ONE.with_precision(prec).value();

        let num = &(&two * &s1) - &(&three * &s2) + &four;
        let den = &(&five * &s1) + &one;
        let k = &num / &den;

        let k_s1 = &k * &s1;
        let x = vec![k_s1, k, s1, s2, one];

        let rel = find_integer_relation_generic(&x).expect("5D relation must be found");
        assert_eq!(rel.m, vec![5, 1, -2, 3, -4]);
        assert!(rel.norm <= MAX_NORM);
        let thresh = DBig::from_parts(IBig::ONE, VERIFY_EXPONENT).with_precision(prec).value();
        assert!(rel.residual_dbig < thresh);
    }

    #[test]
    fn test_known_artificial_relation_6d() {
        let prec = 160;
        let s1 = DBig::from_str(PI_STR).unwrap().with_precision(prec).value();
        let s2 = DBig::from_str(LN2_STR).unwrap().with_precision(prec).value();
        let s1_s2 = &s1 * &s2;

        // 7*K - 2*S1 + 3*S2 - 4*(S1*S2) + 5 = 0
        // K = (2*S1 - 3*S2 + 4*(S1*S2) - 5) / 7
        let seven = DBig::from_parts(IBig::from(7), 0).with_precision(prec).value();
        let two = DBig::from_parts(IBig::from(2), 0).with_precision(prec).value();
        let three = DBig::from_parts(IBig::from(3), 0).with_precision(prec).value();
        let four = DBig::from_parts(IBig::from(4), 0).with_precision(prec).value();
        let five = DBig::from_parts(IBig::from(5), 0).with_precision(prec).value();
        let one = DBig::ONE.with_precision(prec).value();

        let num = &(&(&two * &s1) - &(&three * &s2)) + &(&four * &s1_s2) - &five;
        let k = &num / &seven;

        let x = vec![k, s1, s2, s1_s2, one];

        let rel = find_integer_relation_generic(&x).expect("6D relation must be found");
        assert_eq!(rel.m, vec![7, -2, 3, -4, 5]);
        assert!(rel.norm <= MAX_NORM);
        let thresh = DBig::from_parts(IBig::ONE, VERIFY_EXPONENT).with_precision(prec).value();
        assert!(rel.residual_dbig < thresh);
    }

    #[test]
    fn test_known_artificial_relation_7d() {
        let prec = 160;
        let s1 = DBig::from_str(PI_STR).unwrap().with_precision(prec).value();
        let s2 = DBig::from_str(LN2_STR).unwrap().with_precision(prec).value();

        let b0 = DBig::from_str("1.41421356237309504880168872420969807856967187537694807317667973799").unwrap().with_precision(prec).value();
        let b1 = DBig::from_str("1.73205080756887729352744634150587236694280525381038062805580697945").unwrap().with_precision(prec).value();
        let b2 = DBig::from_str("2.23606797749978969640917366873127623544061835961152572427089724541").unwrap().with_precision(prec).value();
        let t1 = s1;
        let t2 = s2;
        let one = DBig::ONE.with_precision(prec).value();

        // solve for B3:
        // 4*B3 = 3*B0 - 2*B1 + B2 + 5*T1 - 6*T2 + 7
        let four = DBig::from_parts(IBig::from(4), 0).with_precision(prec).value();
        let three = DBig::from_parts(IBig::from(3), 0).with_precision(prec).value();
        let two = DBig::from_parts(IBig::from(2), 0).with_precision(prec).value();
        let five = DBig::from_parts(IBig::from(5), 0).with_precision(prec).value();
        let six = DBig::from_parts(IBig::from(6), 0).with_precision(prec).value();
        let seven = DBig::from_parts(IBig::from(7), 0).with_precision(prec).value();

        let num = &(&(&(&three * &b0) - &(&two * &b1)) + &b2) + &(&five * &t1) - &(&six * &t2) + &seven;
        let b3 = &num / &four;

        let x = vec![b0, b1, b2, b3, t1, t2, one];
        let rel = find_integer_relation_generic(&x).expect("7D relation must be found");
        assert_eq!(rel.m, vec![3, -2, 1, -4, 5, -6, 7]);
        assert!(rel.norm <= MAX_NORM);
        let thresh = DBig::from_parts(IBig::ONE, VERIFY_EXPONENT).with_precision(prec).value();
        assert!(rel.residual_dbig < thresh);

        let x_f64 = [
            x[0].to_f64().value(),
            x[1].to_f64().value(),
            x[2].to_f64().value(),
            x[3].to_f64().value(),
            x[4].to_f64().value(),
            x[5].to_f64().value(),
            x[6].to_f64().value(),
        ];
        let quick_m = quick_check_relation_7d(&x_f64).expect("quick check 7d must find candidate");
        assert_eq!(quick_m, [3, -2, 1, -4, 5, -6, 7]);
    }
}
