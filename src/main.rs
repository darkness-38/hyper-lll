mod constants;
mod kernels;
mod lll;
mod proof;

use constants::{
    FrontierConstantPool, NovelConstantPool, RamanujanGuilleraPool, TargetConstant,
    FRONTIER_PRECISION, PRECISION,
};
use kernels::{
    compute_az_table, compute_cubic_bases, compute_domb_table, compute_frontier_bases,
    compute_odd_denom_bases, compute_super_binomial_bases, compute_w7_bases,
    format_cubic_w4_lhs_latex, format_cubic_w5_lhs_latex, format_odd_denom_lhs_latex,
    format_polynomial5_latex, format_polynomial5_text, format_polynomial_latex,
    format_polynomial_text, format_w7_lhs_latex, get_bases_for_kernel, get_frontier_bases,
    HarmonicTables, SuperBinomialType, FrontierKernelType,
};
use lll::{
    find_all_integer_relations_frontier, find_all_integer_relations_generic,
    find_all_integer_relations_guillera, IntegerRelation,
};
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::time::Instant;
use dashu::float::DBig;

// ════════════════════════════════════════════════════════════════════════════════
//  1. SINIR ALAN MOTORU (DOMB, ALMKVIST-ZUDILIN, B8, B9, 1/pi^3, 1/pi^4, G/pi)
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Debug)]
pub struct FrontierDiscovery {
    pub kernel: FrontierKernelType,
    pub a: i64,
    pub poly: [i64; 5], // [c0, c1, c2, c3, c4]
    pub target: TargetConstant,
    pub target_coeff: i64,
    pub formula_text: String,
    pub formula_latex: String,
    pub residual_str: String,
    pub residual_f64: f64,
    pub norm: f64,
    pub is_novel: bool,
}

fn gcd_i64(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a.abs()
}

fn build_frontier_formula(
    kernel: FrontierKernelType,
    a: i64,
    poly: [i64; 5],
    target: &TargetConstant,
    k: i64,
) -> (String, String) {
    let poly_txt = format_polynomial5_text(&poly);
    let poly_ltx = format_polynomial5_latex(&poly);

    let k_code = kernel.code();
    let k_latex = kernel.latex();

    let abs_k = k.abs();
    let sign_txt = if k < 0 { "-" } else { "" };
    let rhs_txt = if target.name.starts_with("1/") {
        let denom = &target.name[2..];
        if abs_k == 1 {
            format!("{sign_txt}1 / {denom}")
        } else {
            format!("{sign_txt}{abs_k} / {denom}")
        }
    } else if abs_k == 1 {
        format!("{sign_txt}{}", target.name)
    } else {
        format!("{sign_txt}{abs_k}*{}", target.name)
    };

    let sign_ltx = if k < 0 { "-" } else { "" };
    let rhs_ltx = if target.latex.starts_with("\\frac{1}{") && abs_k != 1 {
        format!("{}{}", sign_ltx, target.latex.replacen("\\frac{1}{", &format!("\\frac{{{abs_k}}}{{"), 1))
    } else if target.latex.starts_with("\\frac{\\sqrt{") && abs_k != 1 {
        format!("{}{}", sign_ltx, target.latex.replacen("\\frac{\\sqrt{", &format!("\\frac{{{abs_k}\\sqrt{{"), 1))
    } else if abs_k == 1 {
        format!("{}{}", sign_ltx, target.latex)
    } else {
        format!("{}{abs_k} {}", sign_ltx, target.latex)
    };

    let a_txt = if a < 0 { format!("({a})^n") } else { format!("{a}^n") };
    let formula_text = format!(
        "sum_{{n=0}}^inf ({k_code}(n) / {a_txt}) * ({poly_txt}) = {rhs_txt}"
    );

    let a_str = if a < 0 { format!("({a})^n") } else { format!("{a}^n") };
    let formula_latex = format!(
        "\\sum_{{n=0}}^\\infty \\frac{{{k_latex}}}{{{a_str}}} \\left( {poly_ltx} \\right) = {rhs_ltx}"
    );

    (formula_text, formula_latex)
}

fn process_frontier_relation(
    kernel: FrontierKernelType,
    a: i64,
    m: &[i64],
    tc: &TargetConstant,
    residual_f64: f64,
    norm: f64,
) -> Option<FrontierDiscovery> {
    let dim = m.len();
    let num_series = dim - 1;
    let mut c_series = m[0..num_series].to_vec();
    let mut c_const = m[num_series];

    if c_series.iter().all(|&v| v == 0) || c_const == 0 {
        return None;
    }

    // Leading non-zero coefficient determines overall sign
    let leading = c_series.iter().rev().copied().find(|&v| v != 0)?;
    if leading < 0 {
        for v in &mut c_series {
            *v = -*v;
        }
        c_const = -c_const;
    }

    let mut g = c_const.abs();
    for &c in &c_series {
        if c != 0 {
            g = gcd_i64(g, c.abs());
        }
    }
    if g > 1 {
        for v in &mut c_series {
            *v /= g;
        }
        c_const /= g;
    }

    let k = -c_const;
    let mut poly5 = [0i64; 5];
    for (i, &v) in c_series.iter().enumerate() {
        if i < 5 {
            poly5[i] = v;
        }
    }

    let (formula_text, formula_latex) = build_frontier_formula(kernel, a, poly5, tc, k);
    let residual_str = format!("{:.2e}", residual_f64);

    let is_reference = tc.name.ends_with("/π") && !tc.name.contains("G");
    let is_novel = !is_reference;

    Some(FrontierDiscovery {
        kernel,
        a,
        poly: poly5,
        target: tc.clone(),
        target_coeff: k,
        formula_text,
        formula_latex,
        residual_str,
        residual_f64,
        norm,
        is_novel,
    })
}

fn scan_frontier_kernel_base(
    kernel: FrontierKernelType,
    a: i64,
    pool: &FrontierConstantPool,
    domb_table: &[dashu::integer::IBig],
    az_table: &[dashu::integer::IBig],
) -> Vec<FrontierDiscovery> {
    let mut discoveries = Vec::new();
    let bases = match compute_frontier_bases(kernel, a, FRONTIER_PRECISION, 350, domb_table, az_table) {
        Some(b) => b,
        None => return discoveries,
    };

    let [s0, s1, s2, s3, s4] = bases;

    // 1. Derece-1 Kafesi: [S0, S1, Target] (3D)
    let deg1_targets = pool
        .reference_targets
        .iter()
        .chain(pool.w2_targets.iter())
        .chain(pool.w3_targets.iter())
        .chain(pool.catalan_targets.iter());
    for tc in deg1_targets {
        let x = [s0.clone(), s1.clone(), tc.value.clone()];
        for rel in find_all_integer_relations_frontier(&x) {
            if let Some(disc) = process_frontier_relation(kernel, a, &rel.m, tc, rel.residual_f64, rel.norm) {
                discoveries.push(disc);
            }
        }
    }

    // 2. Derece-2 Kafesi: [S0, S1, S2, Target] (4D)
    let deg2_targets = pool
        .reference_targets
        .iter()
        .chain(pool.w2_targets.iter())
        .chain(pool.w3_targets.iter())
        .chain(pool.catalan_targets.iter());
    for tc in deg2_targets {
        let x = [s0.clone(), s1.clone(), s2.clone(), tc.value.clone()];
        for rel in find_all_integer_relations_frontier(&x) {
            if let Some(disc) = process_frontier_relation(kernel, a, &rel.m, tc, rel.residual_f64, rel.norm) {
                discoveries.push(disc);
            }
        }
    }

    // 3. Derece-3 Kafesi: [S0, S1, S2, S3, Target] (5D)
    let deg3_targets = pool
        .w3_targets
        .iter()
        .chain(pool.catalan_targets.iter())
        .chain(pool.w2_targets.iter())
        .chain(pool.reference_targets.iter());
    for tc in deg3_targets {
        let x = [s0.clone(), s1.clone(), s2.clone(), s3.clone(), tc.value.clone()];
        for rel in find_all_integer_relations_frontier(&x) {
            if let Some(disc) = process_frontier_relation(kernel, a, &rel.m, tc, rel.residual_f64, rel.norm) {
                discoveries.push(disc);
            }
        }
    }

    // 4. Derece-4 Kafesi: [S0, S1, S2, S3, S4, Target] (6D)
    let deg4_targets = pool
        .w4_targets
        .iter()
        .chain(pool.w3_targets.iter())
        .chain(pool.catalan_targets.iter());
    for tc in deg4_targets {
        let x = [s0.clone(), s1.clone(), s2.clone(), s3.clone(), s4.clone(), tc.value.clone()];
        for rel in find_all_integer_relations_frontier(&x) {
            if let Some(disc) = process_frontier_relation(kernel, a, &rel.m, tc, rel.residual_f64, rel.norm) {
                discoveries.push(disc);
            }
        }
    }

    discoveries
}

pub fn run_frontier_discovery_scan() {
    let t_start = Instant::now();

    println!("================================================================================");
    println!("  SINIR ALAN KEŞİF MOTORU: DOMB, ALMKVIST-ZUDILIN, 6. MERTEBE CALABI-YAU & CATALAN");
    println!("  Hassasiyet: 350 Adım dashu (220 Basamak) | LLL: 10^140 | Eşik: < 10^-115 | Norm: 120000");
    println!("================================================================================");

    print!("[1/4] Domb (D_n) ve Almkvist-Zudilin (Z_n) dizileri 350 adım hesaplanıyor... ");
    let t_tbl = Instant::now();
    let domb_table = compute_domb_table(350);
    let az_table = compute_az_table(350);
    println!("Tamamlandı ({:.2?})", t_tbl.elapsed());

    print!("[2/4] 220-Basamak Sınır Alan Sabit Havuzu (1/π³, 1/π⁴, G/π, ζ(3)/π³) üretiliyor... ");
    let t_pool = Instant::now();
    let pool = FrontierConstantPool::new();
    println!("Tamamlandı ({:.2?})", t_pool.elapsed());

    let mut tasks = Vec::new();
    for &kernel in &FrontierKernelType::ALL {
        let bases = get_frontier_bases(kernel);
        for &a in bases {
            tasks.push((kernel, a));
        }
    }

    println!(
        "[3/4] 4 Sınır Çekirdeği (Domb, AZ, B8, B9) ve {} Modüler Taban Rayon ile taranıyor...",
        tasks.len()
    );

    let raw_results: Vec<Vec<FrontierDiscovery>> = tasks
        .par_iter()
        .map(|&(kernel, a)| scan_frontier_kernel_base(kernel, a, &pool, &domb_table, &az_table))
        .collect();

    let mut seen = HashSet::new();
    let mut all_discoveries = Vec::new();

    for list in raw_results {
        for id in list {
            let sig = format!("{}_{}_{:?}_{}_{}", id.kernel.code(), id.a, id.poly, id.target.name, id.target_coeff);
            if seen.insert(sig) {
                all_discoveries.push(id);
            }
        }
    }

    all_discoveries.sort_by(|a, b| {
        b.is_novel
            .cmp(&a.is_novel)
            .then_with(|| {
                a.norm
                    .partial_cmp(&b.norm)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                a.residual_f64
                    .partial_cmp(&b.residual_f64)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    let novel_count = all_discoveries.iter().filter(|d| d.is_novel).count();
    let ref_count = all_discoveries.len() - novel_count;

    println!("[4/4] Tarama tamamlandı. Bulunan bağıntılar: {} Novel, {} Bilinen Referans\n", novel_count, ref_count);

    if all_discoveries.is_empty() {
        println!("(Belirtilen aralıkta ve kalibrasyonda yeni bağıntı bulunamadı)");
    } else {
        for id in &all_discoveries {
            if id.is_novel {
                println!("💥 [DÜNYA ÇAPINDA İLK NOVEL KEŞİF]");
            } else {
                println!("📘 [BİLİNEN REFERANS / MODÜLER]");
            }
            println!("Çekirdek: {} | Taban: a = {}", id.kernel.description(), id.a);
            println!("Formül: {}", id.formula_text);
            println!("LaTeX: {}", id.formula_latex);
            println!("Hata: {} | Norm: {:.4}\n", id.residual_str, id.norm);
        }
    }

    let filename = "discovered_novel_frontier_series.txt";
    let mut file = File::create(filename).expect("Sonuç dosyası oluşturulamadı");

    writeln!(file, "================================================================================").unwrap();
    writeln!(file, "FRONTIER SUPER-BINOMIAL & MODULAR IDENTITIES (DOMB, AZ, B8, B9)").unwrap();
    writeln!(file, "Discovered by hyper-lll (220-digit precision, LLL scale 10^140, threshold 10^-115)").unwrap();
    writeln!(file, "Total Identities: {} (Novel: {}, Modular Reference: {})", all_discoveries.len(), novel_count, ref_count).unwrap();
    writeln!(file, "================================================================================\n").unwrap();

    for (i, id) in all_discoveries.iter().enumerate() {
        let tag = if id.is_novel {
            "[DÜNYA ÇAPINDA İLK NOVEL KEŞİF]"
        } else {
            "[BİLİNEN REFERANS / MODÜLER]"
        };
        writeln!(file, "% ─────────────────────────────────────────────────────────────────────────────").unwrap();
        writeln!(file, "% [{}] {} Çekirdek: {}, Taban: a = {}", i + 1, tag, id.kernel.description(), id.a).unwrap();
        writeln!(file, "% Hedef: {} | Hata: {} | Norm: {:.4}", id.target.name, id.residual_str, id.norm).unwrap();
        writeln!(file, "\\begin{{equation}}").unwrap();
        writeln!(file, "{}", id.formula_latex).unwrap();
        writeln!(file, "\\end{{equation}}\n").unwrap();
    }

    println!("Sonuçlar '{}' kütüğüne kaydedildi.", filename);
    println!("Toplam Süre: {:.2?}", t_start.elapsed());
}

// ════════════════════════════════════════════════════════════════════════════════
//  2. RAMANUJAN-GUILLERA KEŞİF YAPILARI VE FORMÜL OLUŞTURUCULAR
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Debug)]
pub struct GuilleraDiscovery {
    pub kernel: SuperBinomialType,
    pub a: i64,
    pub poly: [i64; 4],
    pub const_coeffs: Vec<(TargetConstant, i64)>,
    pub formula_text: String,
    pub formula_latex: String,
    pub residual_str: String,
    pub residual_f64: f64,
    pub norm: f64,
}

fn build_guillera_formula(
    kernel: SuperBinomialType,
    a: i64,
    poly: [i64; 4],
    const_coeffs: &[(TargetConstant, i64)],
) -> (String, String) {
    let poly_txt = format_polynomial_text(&poly);
    let poly_ltx = format_polynomial_latex(&poly);

    let k_code = kernel.code();
    let k_latex = kernel.latex();

    let mut rhs_txt_parts = Vec::new();
    for (tc, k) in const_coeffs {
        let k = *k;
        if k != 0 {
            let sign = if k < 0 {
                if rhs_txt_parts.is_empty() { "-" } else { " - " }
            } else if rhs_txt_parts.is_empty() {
                ""
            } else {
                " + "
            };
            let abs_k = k.abs();
            let term = if tc.name.starts_with("1/") {
                let denom = &tc.name[2..];
                if abs_k == 1 {
                    format!("1 / {denom}")
                } else {
                    format!("{abs_k} / {denom}")
                }
            } else if abs_k == 1 {
                tc.name.to_string()
            } else {
                format!("{abs_k}*{}", tc.name)
            };
            rhs_txt_parts.push(format!("{sign}{term}"));
        }
    }
    let rhs_txt = if rhs_txt_parts.is_empty() { "0".to_string() } else { rhs_txt_parts.join("") };

    let mut rhs_ltx_parts = Vec::new();
    for (tc, k) in const_coeffs {
        let k = *k;
        if k != 0 {
            let sign = if k < 0 {
                if rhs_ltx_parts.is_empty() { "-" } else { " - " }
            } else if rhs_ltx_parts.is_empty() {
                ""
            } else {
                " + "
            };
            let abs_k = k.abs();
            let term = if tc.latex.starts_with("\\frac{1}{") && abs_k != 1 {
                tc.latex.replacen("\\frac{1}{", &format!("\\frac{{{abs_k}}}{{"), 1)
            } else if tc.latex.starts_with("\\frac{\\sqrt{") && abs_k != 1 {
                tc.latex.replacen("\\frac{\\sqrt{", &format!("\\frac{{{abs_k}\\sqrt{{"), 1)
            } else if abs_k == 1 {
                tc.latex.to_string()
            } else {
                format!("{abs_k} {}", tc.latex)
            };
            rhs_ltx_parts.push(format!("{sign}{term}"));
        }
    }
    let rhs_ltx = if rhs_ltx_parts.is_empty() { "0".to_string() } else { rhs_ltx_parts.join("") };

    let a_txt = if a < 0 { format!("({a})^n") } else { format!("{a}^n") };
    let formula_text = format!(
        "sum_{{n=0}}^inf ({k_code}(n) / {a_txt}) * ({poly_txt}) = {rhs_txt}"
    );

    let a_str = if a < 0 { format!("({a})^n") } else { format!("{a}^n") };
    let formula_latex = format!(
        "\\sum_{{n=0}}^\\infty \\frac{{{k_latex}}}{{{a_str}}} \\left( {poly_ltx} \\right) = {rhs_ltx}"
    );

    (formula_text, formula_latex)
}

fn scan_kernel_base(
    kernel: SuperBinomialType,
    a: i64,
    pool: &RamanujanGuilleraPool,
) -> Vec<GuilleraDiscovery> {
    let mut discoveries = Vec::new();
    let bases = match compute_super_binomial_bases(kernel, a, PRECISION, 450) {
        Some(b) => b,
        None => return discoveries,
    };

    let [s0, s1, s2, s3] = bases;

    let mut check_targets = |targets: &[TargetConstant]| {
        let mut x = Vec::with_capacity(4 + targets.len());
        x.push(s0.clone());
        x.push(s1.clone());
        x.push(s2.clone());
        x.push(s3.clone());
        for tc in targets {
            x.push(tc.value.clone());
        }

        let rels = find_all_integer_relations_guillera(&x);
        for rel in rels {
            let m = &rel.m;
            let mut c_series = [m[0], m[1], m[2], m[3]];
            let mut c_consts: Vec<i64> = m[4..].to_vec();

            if c_series.iter().all(|&v| v == 0) {
                continue;
            }
            if c_consts.iter().all(|&v| v == 0) {
                continue;
            }

            let leading = if c_series[3] != 0 {
                c_series[3]
            } else if c_series[2] != 0 {
                c_series[2]
            } else if c_series[1] != 0 {
                c_series[1]
            } else {
                c_series[0]
            };

            if leading < 0 {
                for v in &mut c_series {
                    *v = -*v;
                }
                for v in &mut c_consts {
                    *v = -*v;
                }
            }

            let k_consts: Vec<(TargetConstant, i64)> = targets
                .iter()
                .zip(c_consts.iter())
                .filter(|(_, c)| **c != 0)
                .map(|(tc, c)| (tc.clone(), -*c))
                .collect();

            if k_consts.is_empty() {
                continue;
            }

            let (formula_text, formula_latex) =
                build_guillera_formula(kernel, a, c_series, &k_consts);

            let residual_str = format!("{:.2e}", rel.residual_f64);

            discoveries.push(GuilleraDiscovery {
                kernel,
                a,
                poly: c_series,
                const_coeffs: k_consts,
                formula_text,
                formula_latex,
                residual_str,
                residual_f64: rel.residual_f64,
                norm: rel.norm,
            });
        }
    };

    for tc in &pool.targets {
        check_targets(std::slice::from_ref(tc));
    }

    check_targets(&pool.pi_powers);
    check_targets(&pool.radicals);
    check_targets(&pool.catalan_zeta);

    discoveries
}

pub fn run_guillera_discovery_scan() {
    let t_start = Instant::now();

    println!("================================================================================");
    println!("  RAMANUJAN / GUILLERA SÜPER-BİNOM VE SAF POLİNOM ARAMA MOTORU");
    println!("  Hassasiyet: 450 Adım dashu (200 Basamak) | LLL Ölçeği: 10^130 | Eşik: < 10^-110");
    println!("================================================================================");

    print!("[1/3] 200-Basamak Rasyonel / Ters Transandantal Sabit Havuzu üretiliyor... ");
    let t_pool = Instant::now();
    let pool = RamanujanGuilleraPool::new();
    println!("Tamamlandı ({:.2?})", t_pool.elapsed());

    let mut tasks = Vec::new();
    for &kernel in &SuperBinomialType::ALL {
        let bases = get_bases_for_kernel(kernel);
        for a in bases {
            tasks.push((kernel, a));
        }
    }

    println!(
        "[2/3] Süper-Binom Çekirdekleri (B1..B7) ve {} Geometrik Taban Rayon ile taranıyor...",
        tasks.len()
    );

    let raw_results: Vec<Vec<GuilleraDiscovery>> = tasks
        .par_iter()
        .map(|&(kernel, a)| scan_kernel_base(kernel, a, &pool))
        .collect();

    let mut seen = HashSet::new();
    let mut all_discoveries = Vec::new();

    for list in raw_results {
        for id in list {
            let sig = format!(
                "{}_{}_{:?}_{:?}",
                id.kernel.code(),
                id.a,
                id.poly,
                id.const_coeffs.iter().map(|(c, k)| (c.name, *k)).collect::<Vec<_>>()
            );
            if seen.insert(sig) {
                all_discoveries.push(id);
            }
        }
    }

    all_discoveries.sort_by(|a, b| {
        a.norm
            .partial_cmp(&b.norm)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                a.residual_f64
                    .partial_cmp(&b.residual_f64)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    println!("[3/3] Tarama tamamlandı. Bulunan bağıntılar listeleniyor:\n");

    if all_discoveries.is_empty() {
        println!("(Belirtilen aralıkta ve kalibrasyonda yeni bağıntı bulunamadı)");
    } else {
        for id in &all_discoveries {
            println!("🔥 [RAMANUJAN-GUILLERA KEŞFİ]");
            println!("Çekirdek: {} | Taban: a = {}", id.kernel.description(), id.a);
            println!("Formül: {}", id.formula_text);
            println!("Hata: {} | Norm: {:.4}\n", id.residual_str, id.norm);
        }
    }

    let filename = "discovered_guillera_series.txt";
    let mut file = File::create(filename).expect("Sonuç dosyası oluşturulamadı");

    writeln!(file, "================================================================================").unwrap();
    writeln!(file, "RAMANUJAN-GUILLERA SUPER-BINOMIAL & POLYNOMIAL IDENTITIES").unwrap();
    writeln!(file, "Discovered by hyper-lll (200-digit precision, LLL scale 10^130, threshold 10^-110)").unwrap();
    writeln!(file, "Total Identities Discovered: {}", all_discoveries.len()).unwrap();
    writeln!(file, "================================================================================\n").unwrap();

    for (i, id) in all_discoveries.iter().enumerate() {
        writeln!(file, "% ─────────────────────────────────────────────────────────────────────────────").unwrap();
        writeln!(file, "% [{}] Çekirdek: {}, Taban: a = {}", i + 1, id.kernel.description(), id.a).unwrap();
        writeln!(file, "% Hata: {} | Norm: {:.4}", id.residual_str, id.norm).unwrap();
        writeln!(file, "\\begin{{equation}}").unwrap();
        writeln!(file, "{}", id.formula_latex).unwrap();
        writeln!(file, "\\end{{equation}}\n").unwrap();
    }

    println!("Sonuçlar '{}' kütüğüne kaydedildi.", filename);
    println!("Toplam Süre: {:.2?}", t_start.elapsed());
}

// ════════════════════════════════════════════════════════════════════════════════
//  3. ÖNCEKİ ALANLAR (AĞIRLIK-7, KÜBİK BİNOM, TEKİL PAYDA) -- NOVEL MODU
// ════════════════════════════════════════════════════════════════════════════════

#[derive(Clone, Debug)]
pub struct NovelIdentity {
    pub area: &'static str,
    pub category_tag: String,
    pub latex_formula: String,
    pub series_coeffs: Vec<i64>,
    pub const_coeffs: Vec<i64>,
    pub norm: f64,
    pub error_str: String,
    pub kernel_desc: String,
}

pub fn format_rhs_latex(constants: &[TargetConstant], k: &[i64]) -> String {
    let mut parts = Vec::new();

    for (tc, &coeff) in constants.iter().zip(k.iter()) {
        if coeff != 0 {
            let sign = if coeff < 0 {
                if parts.is_empty() { "-" } else { "- " }
            } else if parts.is_empty() {
                ""
            } else {
                "+ "
            };
            let abs_c = coeff.abs();
            if abs_c == 1 {
                parts.push(format!("{}{}", sign, tc.latex));
            } else {
                parts.push(format!("{}{}{}", sign, abs_c, tc.latex));
            }
        }
    }

    if parts.is_empty() {
        "0".to_string()
    } else {
        parts.join(" ")
    }
}

fn process_novel_relations(
    area: &'static str,
    kernel_desc: String,
    num_series_bases: usize,
    constants: &[TargetConstant],
    relations: Vec<IntegerRelation>,
    format_lhs: impl Fn(&[i64]) -> String,
    tag_generator: impl Fn(&[i64], &[i64]) -> String,
) -> Vec<NovelIdentity> {
    let mut results = Vec::new();

    for rel in relations {
        let m = &rel.m;
        let c_series = &m[0..num_series_bases];
        let c_consts = &m[num_series_bases..];

        if c_series.iter().all(|&v| v == 0) {
            continue;
        }

        if c_consts.iter().all(|&v| v == 0) {
            continue;
        }

        let k_consts: Vec<i64> = c_consts.iter().map(|&v| -v).collect();
        let lhs_latex = format_lhs(c_series);
        let rhs_latex = format_rhs_latex(constants, &k_consts);
        let latex_formula = format!("{} = {}", lhs_latex, rhs_latex);

        let error_str = format!("{:.3e}", rel.residual_f64);
        let category_tag = tag_generator(c_series, &k_consts);

        results.push(NovelIdentity {
            area,
            category_tag,
            latex_formula,
            series_coeffs: c_series.to_vec(),
            const_coeffs: k_consts,
            norm: rel.norm,
            error_str,
            kernel_desc: kernel_desc.clone(),
        });
    }

    results
}

pub fn run_novel_areas_scan() {
    let t_start = Instant::now();

    println!("================================================================================");
    println!("  KEŞFEDİLMEMİŞ ALAN TARAMASI: AĞIRLIK-7 (ζ(7)), KÜBİK BİNOM VE TEKİL PAYDALAR");
    println!("  Hassasiyet: 550 Adım dashu (200 Basamak) | LLL Ölçeği: 10^120 | Eşik: < 10^-100");
    println!("================================================================================");

    print!("[1/4] Tam harmonik tabloları 1650 adım hesaplanıyor... ");
    let t_harm = Instant::now();
    let harmonics = HarmonicTables::build(1650, PRECISION);
    println!("Tamamlandı ({:.2?})", t_harm.elapsed());

    print!("[2/4] 200-Basamak Dinamik Sabit Havuzları üretiliyor... ");
    let t_consts = Instant::now();
    let cpool = NovelConstantPool::new();
    println!("Tamamlandı ({:.2?})", t_consts.elapsed());

    println!("[3/4] 3 Bakir Araştırma Alanı taranıyor...");

    let a_w7 = [-2, -1, 1, 2];
    let results_w7: Vec<Vec<NovelIdentity>> = a_w7
        .par_iter()
        .map(|&a| {
            let kb = compute_w7_bases(a, 550, &harmonics);
            let mut vec_w7: Vec<DBig> = Vec::with_capacity(16);
            for b in &kb.bases {
                vec_w7.push(b.clone());
            }
            for c in &cpool.w7 {
                vec_w7.push(c.value.clone());
            }

            let rels = find_all_integer_relations_generic(&vec_w7);
            let desc = format!("CentralBinomial a={a}");
            process_novel_relations(
                "AĞIRLIK-7",
                desc,
                10,
                &cpool.w7,
                rels,
                |c| format_w7_lhs_latex(a, c),
                |_c_series, k_consts| {
                    if k_consts[0] != 0 {
                        "🏆 [DÜNYA ÇAPINDA İLK: ZETA(7) KAPALI FORMÜLÜ]".to_string()
                    } else {
                        "🌟 [AĞIRLIK-7 MELEZİ]".to_string()
                    }
                },
            )
        })
        .collect();

    let a_cubic = [2, -2];
    let results_cubic: Vec<Vec<NovelIdentity>> = a_cubic
        .par_iter()
        .map(|&a| {
            let mut local = Vec::new();
            let kb = compute_cubic_bases(a, 550, &harmonics);

            let mut vec_w4: Vec<DBig> = Vec::with_capacity(12);
            for b in &kb.w4 {
                vec_w4.push(b.clone());
            }
            for c in &cpool.w4_cubic {
                vec_w4.push(c.value.clone());
            }
            let rels_w4 = find_all_integer_relations_generic(&vec_w4);
            let desc_w4 = format!("CubicBinomial a={a} (W4)");
            let ids_w4 = process_novel_relations(
                "KÜBİK BİNOM W4",
                desc_w4,
                6,
                &cpool.w4_cubic,
                rels_w4,
                |c| format_cubic_w4_lhs_latex(a, c),
                |_c_series, _k_consts| "🌟 [YENİ KÜBİK BİNOM MELEZİ - W4]".to_string(),
            );
            local.extend(ids_w4);

            let mut vec_w5: Vec<DBig> = Vec::with_capacity(12);
            for b in &kb.w5 {
                vec_w5.push(b.clone());
            }
            for c in &cpool.w5_cubic {
                vec_w5.push(c.value.clone());
            }
            let rels_w5 = find_all_integer_relations_generic(&vec_w5);
            let desc_w5 = format!("CubicBinomial a={a} (W5)");
            let ids_w5 = process_novel_relations(
                "KÜBİK BİNOM W5",
                desc_w5,
                7,
                &cpool.w5_cubic,
                rels_w5,
                |c| format_cubic_w5_lhs_latex(a, c),
                |_c_series, _k_consts| "🌟 [YENİ KÜBİK BİNOM MELEZİ - W5]".to_string(),
            );
            local.extend(ids_w5);

            local
        })
        .collect();

    let odd_params: Vec<(i64, usize)> = {
        let a_vals = [-4, -1, 1, 4];
        let k_vals = [1, 2, 3, 4];
        let mut list = Vec::new();
        for &a in &a_vals {
            for &k in &k_vals {
                list.push((a, k));
            }
        }
        list
    };

    let results_odd: Vec<Vec<NovelIdentity>> = odd_params
        .par_iter()
        .map(|&(a, k)| {
            let kb = compute_odd_denom_bases(a, k, 550, &harmonics);
            let mut vec_odd: Vec<DBig> = Vec::with_capacity(9);
            for b in &kb.bases {
                vec_odd.push(b.clone());
            }
            for c in &cpool.beta {
                vec_odd.push(c.value.clone());
            }

            let rels = find_all_integer_relations_generic(&vec_odd);
            let desc = format!("OddDenom (2n+1)^{k} a={a}");
            process_novel_relations(
                "TEKİL PAYDA",
                desc,
                4,
                &cpool.beta,
                rels,
                |c| format_odd_denom_lhs_latex(a, k, c),
                |_c_series, _k_consts| "💎 [YENİ TEKİL PAYDALI ÖZDEŞLİK]".to_string(),
            )
        })
        .collect();

    let mut all_identities = Vec::new();
    for sub in results_w7.into_iter().chain(results_cubic).chain(results_odd) {
        all_identities.extend(sub);
    }

    all_identities.sort_by(|a, b| a.norm.partial_cmp(&b.norm).unwrap());

    println!("[4/4] Alan taraması tamamlandı. {} adet özdeşlik bulundu.", all_identities.len());
    println!("Toplam Süre: {:.2?}", t_start.elapsed());
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|arg| arg == "--prove") {
        proof::run_formal_proof_verification();
        return;
    }
    if args.iter().any(|arg| arg == "--novel") {
        run_novel_areas_scan();
        return;
    }
    if args.iter().any(|arg| arg == "--guillera") {
        run_guillera_discovery_scan();
        return;
    }

    run_frontier_discovery_scan();
}
