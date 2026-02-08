"""Run statistical tests and generate summary tables and figures."""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.contingency_tables import Table2x2
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from pathlib import Path

INP = Path("results/tables/master.tsv")


def to_binary_label(cls):
    if cls in ("Pathogenic", "Likely_pathogenic"):
        return 1
    if cls in ("Benign", "Likely_benign"):
        return 0
    return None


def parse_score(x):
    if pd.isna(x):
        return np.nan
    s = str(x)
    if "(" in s and ")" in s:
        try:
            return float(s.split("(")[1].split(")")[0])
        except:
            return np.nan
    try:
        return float(s)
    except:
        return np.nan


def fisher_block(sub):
    sub = sub.dropna(subset=["plddt", "y"]).copy()
    sub["hi_conf"] = sub["plddt"].apply(lambda x: x >= 70)
    P_hi = int(((sub.y == 1) & (sub.hi_conf == True)).sum())
    P_lo = int(((sub.y == 1) & (sub.hi_conf == False)).sum())
    B_hi = int(((sub.y == 0) & (sub.hi_conf == True)).sum())
    B_lo = int(((sub.y == 0) & (sub.hi_conf == False)).sum())
    odds, p = fisher_exact([[P_hi, P_lo], [B_hi, B_lo]])
    ci_low = np.nan
    ci_high = np.nan
    try:
        t = Table2x2([[P_hi, P_lo], [B_hi, B_lo]])
        ci_low, ci_high = t.oddsratio_confint()
    except:
        pass
    return {
        "P_hi": P_hi,
        "P_lo": P_lo,
        "B_hi": B_hi,
        "B_lo": B_lo,
        "odds_ratio": odds,
        "p_value": p,
        "or_ci_low": ci_low,
        "or_ci_high": ci_high
    }


def domain_filter(df, domain):
    if domain == "None":
        return df[df["uniprot_domain_group"] == "None"]
    return df[df["uniprot_domain_group"].astype(str).str.contains(domain, na=False)]

def bootstrap_auc(y, score, n_boot=2000, seed=1):
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    score = np.asarray(score)
    if len(np.unique(y)) < 2:
        return np.nan, np.nan, np.nan
    aucs = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        ys = y[idx]
        ss = score[idx]
        if len(np.unique(ys)) < 2:
            continue
        aucs.append(roc_auc_score(ys, ss))
    if len(aucs) == 0:
        return np.nan, np.nan, np.nan
    auc = roc_auc_score(y, score)
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return auc, lo, hi

def revstat_filter(df, mode):
    s = df["revstat"].astype(str)
    if mode == "expert_panel":
        return df[s.str.contains("reviewed_by_expert_panel", na=False)]
    if mode == "multiple_submitters":
        return df[s.str.contains("multiple_submitters", na=False) & s.str.contains("no_conflicts", na=False)]
    return df


def main():
    Path("results/tables").mkdir(parents=True, exist_ok=True)
    Path("results/figures").mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INP, sep="\t")

    df["y"] = df["class"].apply(to_binary_label)
    bin_df = df[df["y"].notnull()].copy()
    bin_df["y"] = bin_df["y"].astype(int)

    table1 = df.groupby(["class", "revstat"]).size().reset_index(name="n")
    table1.to_csv("results/tables/table1_counts.tsv", sep="\t", index=False)

    uniq = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    uniq_counts = uniq.groupby("class").size().reset_index(name="n_unique")
    uniq_plddt = uniq[uniq["plddt"].notnull()].groupby("class").size().reset_index(name="n_unique_with_plddt")
    uniq_counts = uniq_counts.merge(uniq_plddt, on="class", how="left").fillna(0)
    uniq_counts.to_csv("results/tables/table_unique_counts.tsv", sep="\t", index=False)

    fisher_rows = []
    overall = fisher_block(bin_df)
    overall["group"] = "All"
    fisher_rows.append(overall)
    for d in ["RING", "BRCT", "None"]:
        sub = domain_filter(bin_df, d)
        if len(sub) == 0:
            continue
        row = fisher_block(sub)
        row["group"] = d
        fisher_rows.append(row)
    pd.DataFrame(fisher_rows).to_csv("results/tables/table2_fisher.tsv", sep="\t", index=False)

    rev_rows = []
    for mode in ["expert_panel", "multiple_submitters"]:
        sub = revstat_filter(bin_df, mode)
        if len(sub) == 0:
            continue
        row = fisher_block(sub)
        row["group"] = mode
        rev_rows.append(row)
    pd.DataFrame(rev_rows).to_csv("results/tables/table2_fisher_revstat.tsv", sep="\t", index=False)

    plddt_P = bin_df.loc[bin_df.y == 1, "plddt"].dropna().values
    plddt_B = bin_df.loc[bin_df.y == 0, "plddt"].dropna().values
    if len(plddt_P) > 0 and len(plddt_B) > 0:
        mw_p = mannwhitneyu(plddt_P, plddt_B, alternative="two-sided").pvalue
    else:
        mw_p = np.nan
    pd.DataFrame([{"test": "Mann-Whitney", "p_value": mw_p}]).to_csv(
        "results/tables/plddt_mannwhitney.tsv", sep="\t", index=False
    )

    if "SIFT" in bin_df.columns:
        bin_df["SIFT_score"] = bin_df["SIFT"].apply(parse_score)
    else:
        bin_df["SIFT_score"] = np.nan
    if "PolyPhen" in bin_df.columns:
        bin_df["PPH2_score"] = bin_df["PolyPhen"].apply(parse_score)
    else:
        bin_df["PPH2_score"] = np.nan

    ok = bin_df.dropna(subset=["SIFT_score", "PPH2_score", "plddt", "y"]).copy()
    auc_rows = []
    if len(ok) > 0 and ok["y"].nunique() == 2:
        auc_sift, lo_sift, hi_sift = bootstrap_auc(ok["y"], -ok["SIFT_score"])
        auc_pph, lo_pph, hi_pph = bootstrap_auc(ok["y"], ok["PPH2_score"])
        auc_rows.append({
            "group": "All",
            "auc_sift": auc_sift,
            "auc_sift_ci_low": lo_sift,
            "auc_sift_ci_high": hi_sift,
            "auc_polyphen": auc_pph,
            "auc_polyphen_ci_low": lo_pph,
            "auc_polyphen_ci_high": hi_pph,
            "n": len(ok)
        })
        ok["hi_conf"] = ok["plddt"] >= 70
        for b in ["<70", ">=70"]:
            d = ok[ok["hi_conf"] == (b == ">=70")]
            if len(d) < 20 or d["y"].nunique() != 2:
                continue
            auc_sift, lo_sift, hi_sift = bootstrap_auc(d["y"], -d["SIFT_score"])
            auc_pph, lo_pph, hi_pph = bootstrap_auc(d["y"], d["PPH2_score"])
            auc_rows.append({
                "group": b,
                "auc_sift": auc_sift,
                "auc_sift_ci_low": lo_sift,
                "auc_sift_ci_high": hi_sift,
                "auc_polyphen": auc_pph,
                "auc_polyphen_ci_low": lo_pph,
                "auc_polyphen_ci_high": hi_pph,
                "n": len(d)
            })
        for b in [">=90", "70-89.999", "50-69.999", "<50"]:
            d = ok[ok["plddt_bin"] == b]
            if len(d) < 20 or d["y"].nunique() != 2:
                continue
            auc_sift, lo_sift, hi_sift = bootstrap_auc(d["y"], -d["SIFT_score"])
            auc_pph, lo_pph, hi_pph = bootstrap_auc(d["y"], d["PPH2_score"])
            auc_rows.append({
                "group": b,
                "auc_sift": auc_sift,
                "auc_sift_ci_low": lo_sift,
                "auc_sift_ci_high": hi_sift,
                "auc_polyphen": auc_pph,
                "auc_polyphen_ci_low": lo_pph,
                "auc_polyphen_ci_high": hi_pph,
                "n": len(d)
            })
    pd.DataFrame(auc_rows).to_csv("results/tables/auc_summary.tsv", sep="\t", index=False)

    auc_rev_rows = []
    for mode in ["expert_panel", "multiple_submitters"]:
        sub = revstat_filter(ok, mode)
        if len(sub) < 20 or sub["y"].nunique() != 2:
            continue
        auc_sift, lo_sift, hi_sift = bootstrap_auc(sub["y"], -sub["SIFT_score"])
        auc_pph, lo_pph, hi_pph = bootstrap_auc(sub["y"], sub["PPH2_score"])
        auc_rev_rows.append({
            "group": mode,
            "auc_sift": auc_sift,
            "auc_sift_ci_low": lo_sift,
            "auc_sift_ci_high": hi_sift,
            "auc_polyphen": auc_pph,
            "auc_polyphen_ci_low": lo_pph,
            "auc_polyphen_ci_high": hi_pph,
            "n": len(sub)
        })
    pd.DataFrame(auc_rev_rows).to_csv("results/tables/auc_summary_revstat.tsv", sep="\t", index=False)

    vus = df[df["class"] == "Uncertain_significance"].copy()
    if "SIFT" in vus.columns:
        vus["SIFT_score"] = vus["SIFT"].apply(parse_score)
    else:
        vus["SIFT_score"] = np.nan
    if "PolyPhen" in vus.columns:
        vus["PPH2_score"] = vus["PolyPhen"].apply(parse_score)
    else:
        vus["PPH2_score"] = np.nan
    cand = vus[(vus["SIFT_score"] < 0.05) & (vus["PPH2_score"] > 0.9) & (vus["plddt"] >= 70)]
    cand = cand.sort_values(["PPH2_score", "SIFT_score"], ascending=[False, True]).head(15)
    cand.to_csv("results/tables/vus_candidates_top15.tsv", sep="\t", index=False)

    order = ["Benign", "Likely_benign", "Uncertain_significance", "Likely_pathogenic", "Pathogenic"]
    data = [df.loc[df["class"] == c, "plddt"].dropna().values for c in order]
    plt.figure()
    plt.boxplot(data, tick_labels=order, showfliers=False)
    plt.ylabel("pLDDT (AlphaFold, per-residue)")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig("results/figures/plddt_by_class.png", dpi=300)

    bin_counts = df.groupby(["class", "plddt_bin"]).size().reset_index(name="n")
    pivot = bin_counts.pivot(index="class", columns="plddt_bin", values="n").fillna(0)
    pivot = pivot.reindex(order)
    bin_order = [">=90", "70-89.999", "50-69.999", "<50", "NA"]
    pivot = pivot.reindex(columns=bin_order, fill_value=0)
    fig, ax = plt.subplots()
    x = np.arange(len(pivot.index))
    bottom = np.zeros(len(pivot))
    for b in bin_order:
        vals = pivot[b].values
        ax.bar(x, vals, bottom=bottom, label=b)
        bottom += vals
    ax.set_ylabel("Variant count")
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, rotation=30, ha="right")
    ax.legend(title="pLDDT bin", fontsize=8)
    fig.tight_layout()
    fig.savefig("results/figures/plddt_bin_stacked.png", dpi=300)


if __name__ == "__main__":
    main()
