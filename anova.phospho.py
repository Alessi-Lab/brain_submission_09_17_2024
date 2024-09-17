import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests

if __name__ == '__main__':
    df = pd.read_csv(r"H:\20230713_143813_22_phospho_no-norm_Report.output.B.normalized.tsv", sep="\t")
    columns_map = {}
    print(df.columns)
    for i in df.columns:
        if i not in "PTM_collapse_key	EG.ModifiedPeptide	PEP.StrippedSequence	PG.UniProtIds	PTM_localization".split("\t"):
            spl = i.split(".")
            columns_map[i] = {"condition": ".".join(spl[0:-1]), "replicate": spl[-1]}
    print(columns_map)
    # perform anova analysis on each protein for any group without R1441G in name
    for i, r in df.iterrows():
        groups = {}
        for c in columns_map:
            if "R1441G" not in columns_map[c]["condition"] and "iPD" not in columns_map[c]["condition"] and "B" in columns_map[c]["replicate"]:
                if columns_map[c]["condition"] not in groups:
                    groups[columns_map[c]["condition"]] = []
                if not np.isinf(r[c]):
                    groups[columns_map[c]["condition"]].append(r[c])
        F, p = f_oneway(*groups.values())
        df.at[i, "p-value"] = p

    # perform multiple testing correction
    df["fdr_bh"] = multipletests(df["p-value"], method="fdr_bh")[1]
    df.to_csv(r"H:\discover.ANOVA.all.G2019S.and.controls.phospho.tsv", sep="\t", index=False)

    df = pd.read_csv(r"H:\20230713_143813_22_phospho_no-norm_Report.output.S.normalized.tsv", sep="\t")
    columns_map = {}
    print(df.columns)
    for i in df.columns:
        if i not in "PTM_collapse_key	EG.ModifiedPeptide	PEP.StrippedSequence	PG.UniProtIds	PTM_localization".split("\t"):
            spl = i.split(".")
            columns_map[i] = {"condition": ".".join(spl[0:-1]), "replicate": spl[-1]}
    print(columns_map)
    # perform anova analysis on each protein for any group without R1441G in name
    for i, r in df.iterrows():
        groups = {}
        for c in columns_map:
            if "R1441G" not in columns_map[c]["condition"] and "iPD" not in columns_map[c]["condition"] and "S" in \
                    columns_map[c]["replicate"]:
                if columns_map[c]["condition"] not in groups:
                    groups[columns_map[c]["condition"]] = []
                if not np.isinf(r[c]):
                    groups[columns_map[c]["condition"]].append(r[c])
        F, p = f_oneway(*groups.values())
        df.at[i, "p-value"] = p

    # perform multiple testing correction
    df["fdr_bh"] = multipletests(df["p-value"], method="fdr_bh")[1]
    df.to_csv(r"H:\validation.ANOVA.all.G2019S.and.controls.phospho.tsv", sep="\t", index=False)

    df = pd.read_csv(r"H:\20230713_143813_22_phospho_no-norm_Report.output.D.normalized.tsv", sep="\t")
    columns_map = {}
    print(df.columns)
    for i in df.columns:
        if i not in "PTM_collapse_key	EG.ModifiedPeptide	PEP.StrippedSequence	PG.UniProtIds	PTM_localization".split("\t"):
            spl = i.split(".")
            columns_map[i] = {"condition": ".".join(spl[0:-1]), "replicate": spl[-1]}
    print(columns_map)
    # perform anova analysis on each protein for any group without R1441G in name
    for i, r in df.iterrows():
        groups = {}
        for c in columns_map:
            if "R1441G" in columns_map[c]["condition"] and "iPD" not in columns_map[c]["condition"]:
                if columns_map[c]["condition"] not in groups:
                    groups[columns_map[c]["condition"]] = []
                if not np.isinf(r[c]):
                    groups[columns_map[c]["condition"]].append(r[c])
        F, p = f_oneway(*groups.values())
        df.at[i, "p-value"] = p

    # perform multiple testing correction
    df["fdr_bh"] = multipletests(df["p-value"], method="fdr_bh")[1]
    df.to_csv(r"H:\expand.ANOVA.all.R1441G.and.controls.phospho.tsv", sep="\t", index=False)

    df = pd.read_csv(r"H:\20230713_143813_22_phospho_no-norm_Report.output.meta.normalized.tsv", sep="\t")
    columns_map = {}
    print(df.columns)
    for i in df.columns:
        if i not in "PTM_collapse_key	EG.ModifiedPeptide	PEP.StrippedSequence	PG.UniProtIds	PTM_localization".split(
                "\t"):
            spl = i.split(".")
            columns_map[i] = {"condition": ".".join(spl[0:-1]), "replicate": spl[-1]}
    print(columns_map)
    # perform anova analysis on each protein for any group without R1441G in name
    for i, r in df.iterrows():
        groups = {}
        for c in columns_map:
            if "R1441G" not in columns_map[c]["condition"] and "iPD" not in columns_map[c]["condition"]:
                if columns_map[c]["condition"] not in groups:
                    groups[columns_map[c]["condition"]] = []
                if not np.isinf(r[c]):
                    groups[columns_map[c]["condition"]].append(r[c])
        F, p = f_oneway(*groups.values())
        df.at[i, "p-value"] = p

    # perform multiple testing correction
    df["fdr_bh"] = multipletests(df["p-value"], method="fdr_bh")[1]
    df.to_csv(r"H:\meta.ANOVA.ALL-G2019S.and.C.phospho.tsv", sep="\t", index=False)

    df = pd.read_csv(r"H:\20230713_143813_22_phospho_no-norm_Report.output.meta.normalized.tsv", sep="\t")
    columns_map = {}
    print(df.columns)
    for i in df.columns:
        if i not in "PTM_collapse_key	EG.ModifiedPeptide	PEP.StrippedSequence	PG.UniProtIds	PTM_localization".split(
                "\t"):
            spl = i.split(".")
            columns_map[i] = {"condition": ".".join(spl[0:-1]), "replicate": spl[-1]}
    print(columns_map)
    # perform anova analysis on each protein for any group without R1441G in name
    for i, r in df.iterrows():
        groups = {}
        for c in columns_map:
            if "iPD" not in columns_map[c]["condition"]:
                if columns_map[c]["condition"] not in groups:
                    groups[columns_map[c]["condition"]] = []
                if not np.isinf(r[c]):
                    groups[columns_map[c]["condition"]].append(r[c])
        F, p = f_oneway(*groups.values())
        df.at[i, "p-value"] = p

    # perform multiple testing correction
    df["fdr_bh"] = multipletests(df["p-value"], method="fdr_bh")[1]
    df.to_csv(r"H:\meta.ANOVA.ALL-G2019S.all.R1441G.and.C.phospho.tsv", sep="\t", index=False)




