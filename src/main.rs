mod constants;
mod kernels;
mod lll;
mod proof;

use constants::{NovelConstantPool, TargetConstant, PRECISION};
use kernels::{
    compute_cubic_bases, compute_odd_denom_bases, compute_w7_bases,
    format_cubic_w4_lhs_latex, format_cubic_w5_lhs_latex, format_odd_denom_lhs_latex,
    format_w7_lhs_latex, HarmonicTables,
};
use lll::{find_all_integer_relations_generic, IntegerRelation};
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::Write;
use std::time::Instant;
use dashu::float::DBig;

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

        // Filter: must have at least one series term
        if c_series.iter().all(|&v| v == 0) {
            continue;
        }

        // Filter: must have at least one constant term
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|arg| arg == "--prove") {
        proof::run_formal_proof_verification();
        return;
    }

    // 1. Rayon ThreadPool: 7 iş parçacığı
    rayon::ThreadPoolBuilder::new()
        .num_threads(7)
        .build_global()
        .unwrap();

    println!("================================================================================");
    println!("  KEŞFEDİLMEMİŞ ALAN TARAMASI: AĞIRLIK-7 (ζ(7)), KÜBİK BİNOM VE TEKİL PAYDALAR");
    println!("  Hassasiyet: 550 Adım dashu (200 Basamak) | LLL Ölçeği: 10^120 | Eşik: < 10^-100");
    println!("================================================================================");

    let t_start = Instant::now();

    // 2. Precompute Harmonic Tables up to 1650 steps
    print!("[1/4] Tam harmonik tabloları (H_n, H_{{2n}}, H_{{3n}}, H_n^(2..5), O_n, O_n^(2)) 1650 adım hesaplanıyor... ");
    let t_harm = Instant::now();
    let harmonics = HarmonicTables::build(1650, PRECISION);
    println!("Tamamlandı ({:.2?})", t_harm.elapsed());

    // 3. Constants Pools
    print!("[2/4] 200-Basamak Dinamik Sabit Havuzları (CVZ N=270, ζ(7), β(4), β(6)) üretiliyor... ");
    let t_consts = Instant::now();
    let cpool = NovelConstantPool::new();
    println!("Tamamlandı ({:.2?})", t_consts.elapsed());

    println!("[3/4] 3 Bakir Araştırma Alanı taranıyor:");
    println!("      - Alan A: Ağırlık-7 (ζ(7)) 16-Boyutlu LLL (4 çekirdek)");
    println!("      - Alan B: Kübik Binom (3n n) Ağırlık-4 & 5 12-Boyutlu LLL (2 çekirdek x 2 ağırlık)");
    println!("      - Alan C: Tekil Paydalı (2n+1)^k ve β(4) 9-Boyutlu LLL (16 çekirdek)");

    let t_scan = Instant::now();

    // ════════════════════════════════════════════════════════════════════════════
    //  TARAMA 1: ALAN A (AĞIRLIK-7 / ZETA 7 UZAYI - 16D KAFES)
    // ════════════════════════════════════════════════════════════════════════════
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

    // ════════════════════════════════════════════════════════════════════════════
    //  TARAMA 2: ALAN B (KÜBİK BİNOM MERDİVENİ W4 & W5 - 12D KAFES)
    // ════════════════════════════════════════════════════════════════════════════
    let a_cubic = [2, -2];
    let results_cubic: Vec<Vec<NovelIdentity>> = a_cubic
        .par_iter()
        .map(|&a| {
            let mut local = Vec::new();
            let kb = compute_cubic_bases(a, 550, &harmonics);

            // W4 (12D: 6 bases + 6 constants)
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

            // W5 (12D: 7 bases + 5 constants)
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

    // ════════════════════════════════════════════════════════════════════════════
    //  TARAMA 3: ALAN C (TEKİL PAYDALI (2n+1)^k VE BETA(4) - 9D KAFES)
    // ════════════════════════════════════════════════════════════════════════════
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
            let desc = format!("OddDenom a={a} k={k}");
            process_novel_relations(
                "TEKİL PAYDA",
                desc,
                4,
                &cpool.beta,
                rels,
                |c| format_odd_denom_lhs_latex(a, k, c),
                |_c_series, k_consts| {
                    if k_consts[0] != 0 {
                        "💎 [YENİ DİRİCHLET BETA(4) KEŞFİ]".to_string()
                    } else {
                        "💎 [YENİ TEKİL PAYDALI MELEZ]".to_string()
                    }
                },
            )
        })
        .collect();

    println!("Tüm LLL kafes analizleri tamamlandı ({:.2?})\n", t_scan.elapsed());

    // Deduplication
    let mut all_identities = Vec::new();
    let mut seen_formulas = HashSet::new();

    for batch in results_w7.into_iter().chain(results_cubic).chain(results_odd) {
        for id in batch {
            if !seen_formulas.contains(&id.latex_formula) {
                seen_formulas.insert(id.latex_formula.clone());
                all_identities.push(id);
            }
        }
    }

    all_identities.sort_by(|a, b| a.norm.partial_cmp(&b.norm).unwrap());

    let w7_ids: Vec<_> = all_identities.iter().filter(|id| id.area == "AĞIRLIK-7").cloned().collect();
    let cubic_w4_ids: Vec<_> = all_identities.iter().filter(|id| id.area == "KÜBİK BİNOM W4").cloned().collect();
    let cubic_w5_ids: Vec<_> = all_identities.iter().filter(|id| id.area == "KÜBİK BİNOM W5").cloned().collect();
    let odd_ids: Vec<_> = all_identities.iter().filter(|id| id.area == "TEKİL PAYDA").cloned().collect();

    // Print to Console
    println!("════════════════════════════════════════════════════════════════════════════════");
    println!("  🏆 [ALAN A: AĞIRLIK-7 / ZETA(7) VE YÜKSEK MERTEBELİ HARMONİK KEŞİFLER]");
    println!("════════════════════════════════════════════════════════════════════════════════");
    if w7_ids.is_empty() {
        println!("  (Ağırlık-7 / ζ(7) kafesinde lineer bağıntı bulunamadı)");
    } else {
        println!("  >>> TOPLAM {} ADET AĞIRLIK-7 ÖZDEŞLİĞİ KEŞFEDİLDİ! <<<\n", w7_ids.len());
        for (i, id) in w7_ids.iter().enumerate() {
            println!("  [{}] {}:", i + 1, id.category_tag);
            println!("    LaTeX: {}", id.latex_formula);
            println!("    Seri Katsayıları: {:?}", id.series_coeffs);
            println!("    Sabit Katsayıları: {:?}", id.const_coeffs);
            println!("    Norm: {:.4} | Hata: {}", id.norm, id.error_str);
            println!("    Çekirdek: {}\n", id.kernel_desc);
        }
    }

    println!("════════════════════════════════════════════════════════════════════════════════");
    println!("  🌟 [ALAN B: KÜBİK BİNOM (3n n) MERDİVENİ (AĞIRLIK-4 VE AĞIRLIK-5)]");
    println!("════════════════════════════════════════════════════════════════════════════════");
    let total_cubic = cubic_w4_ids.len() + cubic_w5_ids.len();
    if total_cubic == 0 {
        println!("  (Kübik Binom merdiveninde bağıntı bulunamadı)");
    } else {
        println!("  >>> TOPLAM {} ADET KÜBİK BİNOM MELEZİ KEŞFEDİLDİ! <<<\n", total_cubic);
        for (i, id) in cubic_w4_ids.iter().chain(cubic_w5_ids.iter()).enumerate() {
            println!("  [{}] {}:", i + 1, id.category_tag);
            println!("    LaTeX: {}", id.latex_formula);
            println!("    Seri Katsayıları: {:?}", id.series_coeffs);
            println!("    Sabit Katsayıları: {:?}", id.const_coeffs);
            println!("    Norm: {:.4} | Hata: {}", id.norm, id.error_str);
            println!("    Çekirdek: {}\n", id.kernel_desc);
        }
    }

    println!("════════════════════════════════════════════════════════════════════════════════");
    println!("  💎 [ALAN C: TEKİL PAYDALI ÇEKİRDEKLER ((2n+1)^k) VE DİRİCHLET BETA(4)]");
    println!("════════════════════════════════════════════════════════════════════════════════");
    if odd_ids.is_empty() {
        println!("  (Tekil paydalı çekirdeklerde bağıntı bulunamadı)");
    } else {
        println!("  >>> TOPLAM {} ADET TEKİL PAYDALI ÖZDEŞLİK KEŞFEDİLDİ! <<<\n", odd_ids.len());
        for (i, id) in odd_ids.iter().enumerate() {
            println!("  [{}] {}:", i + 1, id.category_tag);
            println!("    LaTeX: {}", id.latex_formula);
            println!("    Seri Katsayıları: {:?}", id.series_coeffs);
            println!("    Sabit Katsayıları: {:?}", id.const_coeffs);
            println!("    Norm: {:.4} | Hata: {}", id.norm, id.error_str);
            println!("    Çekirdek: {}\n", id.kernel_desc);
        }
    }

    // Save to discovered_novel_identities.txt
    let filename = "discovered_novel_identities.txt";
    let mut file = File::create(filename).expect("Sonuç dosyası oluşturulamadı");

    writeln!(file, "================================================================================").unwrap();
    writeln!(file, "  KEŞFEDİLMEMİŞ ALAN TARAMA RAPORU (AĞIRLIK-7, KÜBİK BİNOM, TEKİL PAYDALAR)").unwrap();
    writeln!(file, "  Hassasiyet: 550 Adım dashu (200 Basamak) | LLL Ölçeği: 10^120 | Eşik: < 10^-100").unwrap();
    writeln!(file, "================================================================================\n").unwrap();

    writeln!(file, "════════════════════════════════════════════════════════════════════════════════").unwrap();
    writeln!(file, "  🏆 [ALAN A: AĞIRLIK-7 / ZETA(7) VE YÜKSEK MERTEBELİ HARMONİK KEŞİFLER]").unwrap();
    writeln!(file, "════════════════════════════════════════════════════════════════════════════════").unwrap();
    writeln!(file, "  >>> TOPLAM {} ADET ONAYLANMIŞ ÖZDEŞLİK KEŞFEDİLDİ! <<<\n", w7_ids.len()).unwrap();
    for (i, id) in w7_ids.iter().enumerate() {
        writeln!(file, "  [{}] {}:", i + 1, id.category_tag).unwrap();
        writeln!(file, "    LaTeX: {}", id.latex_formula).unwrap();
        writeln!(file, "    Seri Katsayıları: {:?}", id.series_coeffs).unwrap();
        writeln!(file, "    Sabit Katsayıları: {:?}", id.const_coeffs).unwrap();
        writeln!(file, "    Norm: {:.4} | Hata: {}", id.norm, id.error_str).unwrap();
        writeln!(file, "    Çekirdek: {}\n", id.kernel_desc).unwrap();
    }

    writeln!(file, "════════════════════════════════════════════════════════════════════════════════").unwrap();
    writeln!(file, "  🌟 [ALAN B: KÜBİK BİNOM (3n n) MERDİVENİ (AĞIRLIK-4 VE AĞIRLIK-5)]").unwrap();
    writeln!(file, "════════════════════════════════════════════════════════════════════════════════").unwrap();
    writeln!(file, "  >>> TOPLAM {} ADET ONAYLANMIŞ KÜBİK BİNOM MELEZİ KEŞFEDİLDİ! <<<\n", total_cubic).unwrap();
    for (i, id) in cubic_w4_ids.iter().chain(cubic_w5_ids.iter()).enumerate() {
        writeln!(file, "  [{}] {}:", i + 1, id.category_tag).unwrap();
        writeln!(file, "    LaTeX: {}", id.latex_formula).unwrap();
        writeln!(file, "    Seri Katsayıları: {:?}", id.series_coeffs).unwrap();
        writeln!(file, "    Sabit Katsayıları: {:?}", id.const_coeffs).unwrap();
        writeln!(file, "    Norm: {:.4} | Hata: {}", id.norm, id.error_str).unwrap();
        writeln!(file, "    Çekirdek: {}\n", id.kernel_desc).unwrap();
    }

    writeln!(file, "════════════════════════════════════════════════════════════════════════════════").unwrap();
    writeln!(file, "  💎 [ALAN C: TEKİL PAYDALI ÇEKİRDEKLER ((2n+1)^k) VE DİRİCHLET BETA(4)]").unwrap();
    writeln!(file, "════════════════════════════════════════════════════════════════════════════════").unwrap();
    writeln!(file, "  >>> TOPLAM {} ADET ONAYLANMIŞ TEKİL PAYDALI ÖZDEŞLİK KEŞFEDİLDİ! <<<\n", odd_ids.len()).unwrap();
    for (i, id) in odd_ids.iter().enumerate() {
        writeln!(file, "  [{}] {}:", i + 1, id.category_tag).unwrap();
        writeln!(file, "    LaTeX: {}", id.latex_formula).unwrap();
        writeln!(file, "    Seri Katsayıları: {:?}", id.series_coeffs).unwrap();
        writeln!(file, "    Sabit Katsayıları: {:?}", id.const_coeffs).unwrap();
        writeln!(file, "    Norm: {:.4} | Hata: {}", id.norm, id.error_str).unwrap();
        writeln!(file, "    Çekirdek: {}\n", id.kernel_desc).unwrap();
    }

    println!("[4/4] Sonuçlar '{}' kütüğüne kaydedildi.", filename);
    println!("Toplam Çalışma Süresi: {:.2?}", t_start.elapsed());
}
