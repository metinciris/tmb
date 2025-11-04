#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Basit TMB Raporlayıcı (Qiagen/CLC VCF'ler için)
- Kullanıcı ayarı gerektirmez; sağlam varsayılanlarla çalışır.
- VCF başlığından örnek adını otomatik algılar (son sütun).
- QUAL/DP/ALT okuma/VAF eşikleri sabittir (aşağıda CONFIG içinde).
- STR/artefakt elemesi: homopolimer ve di-nükleotid tekrarlarını kaba olarak filtreler.
- Çıktı:
  * Hasta başına TXT rapor
  * Tek Excel dosyası: summary (hasta bazlı), variants (detay)
Not: VCF INFO alanında fonksiyonel anotasyon yoksa, klinik/raporlanabilir TMB'yi (non-syn coding) birebir taklit edemez.
Burada panel tabanlı "ham filtreli TMB" hesaplanır.
"""

import os, sys, re, math, argparse, datetime, textwrap
from collections import defaultdict
from dataclasses import dataclass, asdict

# Pandas opsiyonel ama Excel için gerekli
try:
    import pandas as pd
except Exception as e:
    pd = None

# ---------- Sabit Varsayılanlar ----------
CONFIG = {
    "panel_name": "Qiagen/CLC Panel",
    "panel_mb": 1.20,        # Panelin raporlanabilir/hedef alanı (Mb). Kurumunuza göre güncelleyiniz.
    "reference": "hg19",
    "qual_min": 50.0,
    "dp_min": 100,
    "alt_min": 5,
    "vaf_min": 0.05,
    "require_pass": True,    # FILTER alanı PASS ya da '.' olsun
    "drop_str_artifacts": True
}

# ---------- Yardımcılar ----------

def is_homopolymer(seq: str) -> bool:
    """AAA..., CCCC... gibi homopolimer?"""
    if not seq:
        return False
    return len(set(seq)) == 1 and len(seq) >= 3

def is_dinuc_repeat(seq: str) -> bool:
    """ACACAC..., GTGT... gibi di-nükleotid tekrar?"""
    if not seq or len(seq) < 4:
        return False
    if len(seq) % 2 != 0:
        # tek sayılarda da bir derece tekrar olabilir ama kaba filtre
        pass
    unit = seq[:2]
    return seq == (unit * (len(seq) // 2)) + unit[:(len(seq) % 2)]

def looks_like_str_artifact(ref: str, alt: str) -> bool:
    """Basit STR artefakt filtresi: homopolimer ya da di-nükleotid tekrar içeren alleller"""
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    # uzun ALT listeyi virgüle göre ayırıp herhangi biri tekrar ise artefakt say
    alts = alt.split(",") if "," in alt else [alt]
    for allele in [ref] + alts:
        if is_homopolymer(allele) or is_dinuc_repeat(allele):
            return True
    return False

def safe_float(x, default=math.nan):
    try:
        return float(x)
    except Exception:
        return default

def parse_format(fmt: str, sample: str):
    """FORMAT ve örnek sütunundan GT, DP ve CLCAD2'yi çıkar.
    CLCAD2: 'refDepth,alt1Depth,alt2Depth,...'
    DP: toplam filtreli okuma
    VAF hesaplamak için en yüksek alt derinliği seçilir.
    """
    result = {"GT": None, "DP": None, "CLCAD2": None, "ALT_AD_MAX": None, "VAF": None}
    try:
        keys = fmt.split(":") if fmt else []
        vals = sample.split(":") if sample else []
        m = dict(zip(keys, vals))
        gt = m.get("GT")
        dp = m.get("DP")
        ad = m.get("CLCAD2")
        result["GT"] = gt
        result["DP"] = int(dp) if dp is not None and dp != "." else None
        result["CLCAD2"] = ad
        if ad and ad != ".":
            parts = [p for p in ad.split(",") if p != ""]
            # parts[0] ref, geri kalan alt alleller
            if len(parts) >= 2:
                alt_depths = [int(x) for x in parts[1:] if re.fullmatch(r"\d+", x)]
                if alt_depths:
                    alt_max = max(alt_depths)
                    result["ALT_AD_MAX"] = alt_max
                    if result["DP"] and result["DP"] > 0:
                        result["VAF"] = alt_max / result["DP"]
        return result
    except Exception:
        return result

@dataclass
class VarRow:
    sample: str
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    filt: str
    dp: int
    ad_alt_max: int
    vaf: float
    kept: bool
    reason: str

# ---------- VCF Okuyucu ----------

def read_vcf_variants(vcf_path: str):
    """VCF'den varyantları (satırları) üretir. Header'ı çöz, sample sütununu bul."""
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(vcf_path)

    with open(vcf_path, "r", encoding="utf-8", errors="ignore") as f:
        header_cols = None
        sample_name = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_cols = line.strip().lstrip("#").split("\t")
                # beklenen minimum sütunlar
                required = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
                if not all(c in header_cols for c in required):
                    raise ValueError("Beklenen VCF sütunları bulunamadı (#CHROM .. QUAL .. FILTER).")
                # FORMAT ve örnek isimleri varsa tespit et
                if "FORMAT" in header_cols and len(header_cols) > header_cols.index("FORMAT") + 1:
                    sample_name = header_cols[-1]
                break

        if header_cols is None:
            raise ValueError("VCF başlığı (#CHROM ...) bulunamadı.")

    # Tekrar aç ve veri satırlarını işle
    variants = []
    with open(vcf_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            # kolon indexleri
            try:
                idx = {name: i for i, name in enumerate(header_cols)}
                chrom = parts[idx["CHROM"]]
                pos = int(parts[idx["POS"]])
                ref = parts[idx["REF"]]
                alt = parts[idx["ALT"]]
                qual = safe_float(parts[idx["QUAL"]])
                filt = parts[idx["FILTER"]] if "FILTER" in idx else "."
                fmt = parts[idx["FORMAT"]] if "FORMAT" in idx else None
                sample_col = parts[idx[sample_name]] if sample_name else None
            except Exception:
                # hatalı satırları atla
                continue

            parsed = parse_format(fmt, sample_col) if sample_name else {"DP": None, "ALT_AD_MAX": None, "VAF": None, "GT": None}
            dp = parsed.get("DP")
            ad_alt_max = parsed.get("ALT_AD_MAX")
            vaf = parsed.get("VAF")

            # Çoklu ALT'ları string olarak koru (detay sayfasında görünür)
            variants.append({
                "sample": sample_name or os.path.basename(vcf_path),
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "qual": qual, "filt": filt,
                "dp": dp, "ad_alt_max": ad_alt_max, "vaf": vaf
            })

    return sample_name or os.path.basename(vcf_path), variants

# ---------- Filtreleme ve TMB ----------

def keep_variant(v, cfg):
    """Varyantı rapor TMB hesabına dahil etme kriterleri."""
    reasons = []
    # FILTER PASS
    if cfg["require_pass"] and v["filt"] not in (".", "PASS"):
        return False, "FILTER!=PASS"

    # Temel alanlar mevcut mu?
    if v["qual"] is None or math.isnan(v["qual"]):
        return False, "QUAL_missing"
    if v["dp"] is None:
        return False, "DP_missing"
    if v["ad_alt_max"] is None:
        return False, "AD_alt_missing"
    if v["vaf"] is None:
        return False, "VAF_missing"

    # Eşikler
    if v["qual"] < cfg["qual_min"]:
        return False, f"QUAL<{cfg['qual_min']}"
    if v["dp"] < cfg["dp_min"]:
        return False, f"DP<{cfg['dp_min']}"
    if v["ad_alt_max"] < cfg["alt_min"]:
        return False, f"ALT<{cfg['alt_min']}"
    if v["vaf"] < cfg["vaf_min"]:
        return False, f"VAF<{cfg['vaf_min']}"

    # STR artefakt
    if cfg["drop_str_artifacts"] and looks_like_str_artifact(v["ref"], v["alt"]):
        return False, "STR_artifact"

    return True, "OK"

def compute_tmb(kept_count: int, panel_mb: float):
    if panel_mb and panel_mb > 0:
        return kept_count / panel_mb
    return float("nan")

# ---------- Raporlama ----------

def make_text_report(sample_name, kept, dropped, outdir, cfg, vcf_path):
    """Hasta bazlı TXT rapor"""
    tmb = compute_tmb(len(kept), cfg["panel_mb"])
    today = datetime.date.today().isoformat()

    lines = []
    lines.append(f"*** TMB RAPORU ***\n")
    lines.append(f"Tarih: {today}")
    lines.append(f"Örnek: {sample_name}")
    lines.append(f"VCF: {vcf_path}")
    lines.append(f"Panel: {cfg['panel_name']} ({cfg['reference']})")
    lines.append(f"Panel boyutu (Mb): {cfg['panel_mb']:.3f}")
    lines.append("Filtre Eşikleri: QUAL≥{qual_min}, DP≥{dp_min}, ALT≥{alt_min}, VAF≥{vaf_min}".format(**cfg))
    lines.append(f"Filter PASS şartı: {'Evet' if cfg['require_pass'] else 'Hayır'}")
    lines.append(f"STR artefakt elemesi: {'Evet' if cfg['drop_str_artifacts'] else 'Hayır'}")
    lines.append("")
    lines.append(f"Nitelikli varyant sayısı: {len(kept)}")
    lines.append(f"TMB (varyant/Mb): {tmb:.2f}" if not math.isnan(tmb) else "TMB: hesaplanamadı (panel_mb eksik)")
    lines.append("")
    lines.append("Notlar:")
    lines.append("- Bu değer panel tabanlı, filtrelenmiş ham TMB'dir. Fonksiyonel anotasyon yoksa benign/silent dışlama yapılamaz.")
    lines.append("- Klinik raporlama için laboratuvarınızın doğrulanmış eşiği ve panel kapsamı esas alınmalıdır.")
    lines.append("- STR/tekrar bölgelerindeki çağrılar basit sezgisel kuralla elenmiştir.")
    lines.append("")

    # İlk 20 varyant
    lines.append("İlk 20 nitelikli varyant (CHROM:POS REF>ALT | QUAL | DP | ALT_AD_MAX | VAF):")
    for r in kept[:20]:
        vaf_str = f"{r['vaf']:.3f}" if r.get("vaf") is not None and not math.isnan(r["vaf"]) else ""
        lines.append(f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} | {r['qual']:.2f} | {r['dp']} | {r['ad_alt_max']} | {vaf_str}")

    report_path = os.path.join(outdir, f"{sample_name}_TMB_report.txt")
    with open(report_path, "w", encoding="utf-8") as w:
        w.write("\n".join(lines))
    return report_path, tmb

def write_excel(all_rows, summary_rows, outdir):
    if pd is None:
        return None
    df_all = pd.DataFrame(all_rows)
    df_sum = pd.DataFrame(summary_rows)
    xlsx_path = os.path.join(outdir, "tmb_summary_and_variants.xlsx")
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
        df_sum.to_excel(xw, index=False, sheet_name="summary")
        df_all.to_excel(xw, index=False, sheet_name="variants")
    return xlsx_path

# ---------- Çalıştırıcı ----------

def process_vcf(vcf_path: str, outdir: str, cfg: dict):
    sample, variants = read_vcf_variants(vcf_path)
    kept_rows, drop_rows = [], []
    for v in variants:
        keep, reason = keep_variant(v, cfg)
        row = {
            "sample": v["sample"],
            "chrom": v["chrom"],
            "pos": v["pos"],
            "ref": v["ref"],
            "alt": v["alt"],
            "qual": v["qual"],
            "filter": v["filt"],
            "dp": v["dp"],
            "ad_alt_max": v["ad_alt_max"],
            "vaf": v["vaf"],
            "kept": keep,
            "reason": reason
        }
        if keep:
            kept_rows.append(row)
        else:
            drop_rows.append(row)
    txt_path, tmb = make_text_report(sample, kept_rows, drop_rows, outdir, cfg, vcf_path)
    return sample, kept_rows, drop_rows, txt_path, tmb

def parse_args(argv):
    p = argparse.ArgumentParser(description="Basit TMB raporlayıcı (Qiagen/CLC VCF). Varsayılanlarla çalışır.")
    p.add_argument("--vcf", nargs="+", help="Bir veya daha fazla VCF dosyası yolu")
    p.add_argument("--outdir", default=".", help="Çıktı klasörü (varsayılan: .)")
    # gelişmiş ayarlar gerekmeden sabit; ama istenirse override
    p.add_argument("--panel-mb", type=float, default=CONFIG["panel_mb"])
    p.add_argument("--panel-name", default=CONFIG["panel_name"])
    p.add_argument("--reference", default=CONFIG["reference"])
    p.add_argument("--no-pass", action="store_true", help="FILTER=PASS şartını KAPAT")
    p.add_argument("--no-str", action="store_true", help="STR artefakt elemesini KAPAT")
    return p.parse_args(argv)

def main(argv=None):
    # Tk ile dosya seçimi (argüman yoksa)
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        try:
            import tkinter as tk
            from tkinter import filedialog, messagebox
            root = tk.Tk()
            root.withdraw()
            vcf_paths = filedialog.askopenfilenames(title="VCF dosyalarını seçin", filetypes=[("VCF", "*.vcf *.vcf.gz"), ("Tümü","*.*")])
            if not vcf_paths:
                return 1
            outdir = filedialog.askdirectory(title="Çıktı klasörünü seçin") or os.getcwd()
            args = argparse.Namespace(
                vcf=list(vcf_paths),
                outdir=outdir,
                panel_mb=CONFIG["panel_mb"],
                panel_name=CONFIG["panel_name"],
                reference=CONFIG["reference"],
                no_pass=False,
                no_str=False
            )
        except Exception:
            # GUI açılamıyorsa CLI'ye düş
            args = parse_args(argv)
    else:
        args = parse_args(argv)

    cfg = CONFIG.copy()
    cfg["panel_mb"] = args.panel_mb
    cfg["panel_name"] = args.panel_name
    cfg["reference"] = args.reference
    if args.no_pass:
        cfg["require_pass"] = False
    if args.no_str:
        cfg["drop_str_artifacts"] = False

    os.makedirs(args.outdir, exist_ok=True)

    all_rows = []
    summary_rows = []
    txt_paths = []
    for vcf_path in args.vcf or []:
        try:
            sample, kept, dropped, txt_path, tmb_val = process_vcf(vcf_path, args.outdir, cfg)
            txt_paths.append(txt_path)
            all_rows.extend(kept + dropped)
            summary_rows.append({
                "sample": sample,
                "vcf": vcf_path,
                "panel_name": cfg["panel_name"],
                "reference": cfg["reference"],
                "panel_mb": cfg["panel_mb"],
                "kept_variants": len(kept),
                "dropped_variants": len(dropped),
                "tmb_variants_per_mb": round(tmb_val, 3) if not math.isnan(tmb_val) else None,
                "qual_min": cfg["qual_min"],
                "dp_min": cfg["dp_min"],
                "alt_min": cfg["alt_min"],
                "vaf_min": cfg["vaf_min"],
                "require_pass": cfg["require_pass"],
                "drop_str_artifacts": cfg["drop_str_artifacts"]
            })
        except Exception as e:
            # Hatalı dosyayı summary'ye not düş
            summary_rows.append({
                "sample": os.path.basename(vcf_path),
                "vcf": vcf_path,
                "panel_name": cfg["panel_name"],
                "reference": cfg["reference"],
                "panel_mb": cfg["panel_mb"],
                "kept_variants": None,
                "dropped_variants": None,
                "tmb_variants_per_mb": None,
                "error": str(e)
            })

    xlsx_path = write_excel(all_rows, summary_rows, args.outdir) if (pd is not None) else None

    # Konsol özeti
    print("Tamamlandı.")
    if txt_paths:
        for pth in txt_paths:
            print("TXT:", pth)
    if xlsx_path:
        print("Excel:", xlsx_path)

    return 0

if __name__ == "__main__":
    sys.exit(main())
