# -*- coding: utf-8 -*-
"""make_figures_ml.py — 生成 ML 部分的稿件级图（ROC 曲线 + 性能总览）"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import roc_curve, roc_auc_score

BASE = _ROOT
PROC = f"{BASE}/01_数据/处理后数据"
OUT = f"{BASE}/03_结果/ML"
FIG = f"{BASE}/04_图表"
os.makedirs(FIG, exist_ok=True)
SEEDS = list(range(1, 11))


def load_cohort(gse):
    expr = pd.read_csv(f"{PROC}/{gse}_log2_expr.csv", index_col=0)
    ann = pd.read_csv(f"{PROC}/{gse}_annotation.csv", index_col=0)
    return expr, ann


def make_xy(expr, ann, genes, label_col, pos, neg):
    mask = ann[label_col].isin([pos, neg])
    ann = ann[mask]
    used = [g for g in genes if g in expr.index]
    X = expr.loc[used, ann.index].T.values
    y = (ann[label_col] == pos).astype(int).values
    return X, y, used


def fit_lasso(X, y, seed):
    pipe = Pipeline([("sc", StandardScaler()),
                     ("m", LogisticRegression(penalty="l1", solver="liblinear",
                                              max_iter=5000, random_state=seed))])
    gs = GridSearchCV(pipe, {"m__C": [0.01, 0.1, 1.0, 10.0]},
                      cv=StratifiedKFold(3, shuffle=True, random_state=seed),
                      scoring="roc_auc")
    gs.fit(X, y)
    coef = gs.best_estimator_.named_steps["m"].coef_.ravel()
    sel = np.where(np.abs(coef) > 1e-8)[0]
    if len(sel) == 0:
        sel = np.arange(X.shape[1])
    model = Pipeline([("sc", StandardScaler()),
                      ("m", LogisticRegression(penalty="l1", solver="liblinear",
                                               C=gs.best_params_["m__C"], max_iter=5000,
                                               random_state=seed))])
    model.fit(X[:, sel], y)
    return sel, model


def fit_rf(X, y, seed):
    m = RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=1)
    m.fit(X, y)
    sel = np.argsort(-m.feature_importances_)[:max(2, int(0.3 * X.shape[1]))]
    model = RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=1)
    model.fit(X[:, sel], y)
    return sel, model


def fit_svmrfe(X, y, seed):
    k = max(2, int(0.3 * X.shape[1]))
    r = RFE(SVC(kernel="linear", random_state=seed), n_features_to_select=k, step=1)
    r.fit(X, y)
    sel = np.where(r.support_)[0]
    model = SVC(kernel="linear", probability=True, random_state=seed)
    model.fit(X[:, sel], y)
    return sel, model


def fit_xgb(X, y, seed):
    from xgboost import XGBClassifier
    m = XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05,
                      eval_metric="logloss", random_state=seed, n_jobs=1,
                      use_label_encoder=False)
    m.fit(X, y)
    sel = np.argsort(-m.feature_importances_)[:max(2, int(0.3 * X.shape[1]))]
    model = XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05,
                          eval_metric="logloss", random_state=seed, n_jobs=1,
                          use_label_encoder=False)
    model.fit(X[:, sel], y)
    return sel, model


FITTERS = {"LASSO": fit_lasso, "RF": fit_rf, "SVM-RFE": fit_svmrfe, "XGBoost": fit_xgb}


def roc_panel(ax, y_true, scores: dict, title, colors):
    for name, s in scores.items():
        fpr, tpr, _ = roc_curve(y_true, s)
        auc = roc_auc_score(y_true, s)
        ax.plot(fpr, tpr, lw=2, label=f"{name} (AUC={auc:.3f})", color=colors.get(name, "grey"))
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(title, fontsize=11)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)


def main():
    pool = pd.read_csv(f"{BASE}/03_结果/WGCNA/ML_candidate_pool.csv")
    genes = pool["gene"].tolist()
    e431, a431 = load_cohort("GSE134431")
    e999, a999 = load_cohort("GSE199939")

    # Task A 训练
    m = a431["group2"].isin(["DFU", "DFS"])
    Xtr = e431.loc[genes, a431.index[m]].T.values
    ytr = (a431.loc[m, "group2"] == "DFU").astype(int).values

    # 各算法最终模型 + 共识
    per_algo = {}
    for name, fn in FITTERS.items():
        sel, model = fn(Xtr, ytr, 1)
        per_algo[name] = ([genes[i] for i in sel], model)

    # 共识基因（>=3/4 算法跨 10 种子）
    sel_counts = {g: 0 for g in genes}
    for name, fn in FITTERS.items():
        for seed in SEEDS:
            sel, _ = fn(Xtr, ytr, seed)
            for i in sel:
                sel_counts[genes[i]] += 1
    consensus = [g for g, c in sel_counts.items() if c >= 3]

    # 外部分数
    m9 = a999["group2"].isin(["DFU", "control"])
    ye = (a999.loc[m9, "group2"] == "DFU").astype(int).values
    scores = {}
    for name, (sel, model) in per_algo.items():
        common = [g for g in sel if g in e999.index]
        if len(common) < len(sel):
            continue
        Xe = e999.loc[common, a999.index[m9]].T.values
        scores[name] = model.predict_proba(Xe)[:, 1] if hasattr(model, "predict_proba") \
            else model.decision_function(Xe)
    lr = Pipeline([("sc", StandardScaler()),
                   ("m", LogisticRegression(max_iter=5000, random_state=1))])
    sel_idx = [genes.index(g) for g in consensus]
    lr.fit(Xtr[:, sel_idx], ytr)
    Xe_c = e999.loc[consensus, a999.index[m9]].T.values
    scores["Consensus-LR"] = lr.predict_proba(Xe_c)[:, 1]

    # 图 1：外部验证 ROC（GSE199939）
    colors = {"LASSO": "#C0392B", "RF": "#2166AC", "SVM-RFE": "#1B7837",
              "XGBoost": "#7B3294", "Consensus-LR": "#E08214"}
    fig, ax = plt.subplots(figsize=(6.5, 6))
    roc_panel(ax, ye, scores,
              "External validation (GSE199939): DFU vs control", colors)
    fig.tight_layout()
    fig.savefig(f"{FIG}/ML_ROC_external_GSE199939.png", dpi=300)
    plt.close(fig)

    # 图 2：共识基因表达箱线图（GSE134431 + GSE199939）
    n = len(consensus)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.8 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for i, g in enumerate(consensus):
        ax = axes[i]
        for gse, e, a, pos, neg, lab in [
            ("GSE134431", e431, a431, "DFU", "DFS", "DFU"),
            ("GSE199939", e999, a999, "DFU", "control", "Control")]:
            mask = a["group2"].isin([pos, neg])
            vals = e.loc[g, a.index[mask]].values
            grp = a.loc[mask, "group2"].values
            data = [vals[grp == pos], vals[grp == neg]]
            bp = ax.boxplot(data, positions=[0, 1] if gse == "GSE134431" else [0.25, 1.25],
                            widths=0.22, patch_artist=True, showfliers=False)
            for patch, c in zip(bp["boxes"], ["#C0392B", "#2166AC"]):
                patch.set_facecolor(c); patch.set_alpha(0.6)
        ax.set_title(g, fontsize=10)
        ax.set_xticks([0.125, 1.125]); ax.set_xticklabels(["DFU", "Ctrl"], fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
    for j in range(n, len(axes)):
        axes[j].axis("off")
    fig.suptitle("Consensus KpyRN genes: expression in DFU vs control (two cohorts)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{FIG}/ML_consensus_gene_boxplots.png", dpi=300)
    plt.close(fig)

    # 保存 ROC 曲线数据（审计）
    roc_data = {"sample": a999.index[m9].tolist(), "y": ye.tolist()}
    for k, v in scores.items():
        roc_data[k] = np.round(v, 6).tolist()
    pd.DataFrame(roc_data).to_csv(f"{OUT}/TaskA_external_ROC_scores.csv", index=False)

    # 图 3：算法性能总览（CV + 外部）
    perf = pd.read_csv(f"{OUT}/TaskA_diagnosis_performance.csv")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    cv = perf[perf["setting"] == "train_repeated5foldCV"]
    ext = perf[perf["setting"] == "ext_GSE199939"]
    algos = ["LASSO", "RF", "SVM-RFE", "XGBoost", "Consensus-LR"]
    cv_vals = [cv[cv["algorithm"] == a]["auc_mean"].values[0] if len(cv[cv["algorithm"] == a]) else np.nan for a in algos]
    ext_vals = [ext[ext["algorithm"] == a]["auc_mean"].values[0] if len(ext[ext["algorithm"] == a]) else np.nan for a in algos]
    x = np.arange(len(algos))
    ax.bar(x - 0.2, cv_vals, width=0.38, label="Training (repeated 5-fold CV)", color="#2166AC", alpha=0.85)
    ax.bar(x + 0.2, ext_vals, width=0.38, label="External (GSE199939)", color="#C0392B", alpha=0.85)
    ax.axhline(0.5, color="grey", ls="--", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(algos, fontsize=9)
    ax.set_ylabel("AUC"); ax.set_ylim(0, 1.08)
    ax.legend(fontsize=9)
    ax.set_title("Task A (DFU diagnosis): cross-validated vs external performance", fontsize=11)
    for i, (c, e) in enumerate(zip(cv_vals, ext_vals)):
        if np.isfinite(c):
            ax.text(i - 0.2, c + 0.02, f"{c:.2f}", ha="center", fontsize=8)
        if np.isfinite(e):
            ax.text(i + 0.2, e + 0.02, f"{e:.2f}", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{FIG}/ML_performance_overview.png", dpi=300)
    plt.close(fig)
    print("figures saved:", os.listdir(FIG))


if __name__ == "__main__":
    main()
