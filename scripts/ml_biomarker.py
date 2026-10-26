# -*- coding: utf-8 -*-
"""
ml_biomarker.py — 四算法共识机器学习（任务 #6，v2 修复版）
项目：丙酮酰化（Kpy）与糖尿病足（DFU）
v2 修复：
  1) LASSO 增加折内 C 网格嵌套选择（修复 C=0.1 过度正则化导致的 CV nan）
  2) 外部验证做特征对齐：所选基因缺失时记录 skipped，不静默错配
  3) XGBoost 保持标准配置；其单基因支配行为如实记录（多算法共识的意义正在于此）
设计（预注册）：
  Task A（诊断）：训练 GSE134431 (DFU vs DFS) → 外部验证 GSE199939 (DFU vs 对照) 等
  Task B（结局）：训练 GSE134431 (Healer vs Nonhealer) → 外部验证 GSE143735
  算法：LASSO(L1-logistic, 嵌套C) / RF / SVM-RFE / XGBoost，10 个随机种子
  评估：训练集内重复分层 5 折 out-of-fold AUC；外部队列 AUC + DeLong 95%CI
  共识：跨种子被 >=3/4 算法选中的基因
"""
import os

import os as _os
_ROOT = _os.environ.get("KPYRN_PROJECT_ROOT") or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))

import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score

BASE = _ROOT
PROC = f"{BASE}/01_数据/处理后数据"
OUT = f"{BASE}/03_结果/ML"
os.makedirs(OUT, exist_ok=True)
SEEDS = list(range(1, 11))  # 预注册：种子 1-10
LASSO_GRID = {"m__C": [0.01, 0.1, 1.0, 10.0]}


# ---------------- DeLong 95% CI ----------------
def compute_midrank(z):
    J = np.argsort(z)
    Zr = z[J]
    N = len(z)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Zr[j] == Zr[i]:
            j += 1
        T[i:j] = (i + j - 1) / 2.0 + 1
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T
    return T2


def delong_ci(y_true, y_score):
    y_true = np.asarray(y_true)
    order = np.argsort(-y_score)
    y_true, y_score = y_true[order], np.asarray(y_score)[order]
    pos, neg = y_true == 1, ~(y_true == 1)
    n_pos, n_neg = int(pos.sum()), int(neg.sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan, np.nan, np.nan
    auc = roc_auc_score(y_true, y_score)
    px, nx = y_score[pos], y_score[neg]
    sorted_n, sorted_p = np.sort(nx), np.sort(px)
    v10 = np.empty(n_pos); v01 = np.empty(n_neg)
    for i in range(n_pos):
        less = np.searchsorted(sorted_n, px[i], side="left")
        equal = np.searchsorted(sorted_n, px[i], side="right") - less
        v10[i] = (less + 0.5 * equal) / n_neg
    for j in range(n_neg):
        greater = n_pos - np.searchsorted(sorted_p, nx[j], side="right")
        equal = np.searchsorted(sorted_p, nx[j], side="right") - np.searchsorted(sorted_p, nx[j], side="left")
        v01[j] = (greater + 0.5 * equal) / n_pos
    var_auc = np.var(v10) / n_pos + np.var(v01) / n_neg
    se = np.sqrt(max(var_auc, 0.0))
    return auc, max(auc - 1.96 * se, 0.0), min(auc + 1.96 * se, 1.0)


# ---------------- 数据 ----------------
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
    return X, y, used, ann.index.tolist()


# ---------------- 算法：返回 (所选特征索引, 在所选特征上拟合的 model, 参数) ----------------
def fit_lasso(X, y, seed):
    pipe = Pipeline([("sc", StandardScaler()),
                     ("m", LogisticRegression(penalty="l1", solver="liblinear",
                                              max_iter=5000, random_state=seed))])
    gs = GridSearchCV(pipe, LASSO_GRID,
                      cv=StratifiedKFold(3, shuffle=True, random_state=seed),
                      scoring="roc_auc")
    gs.fit(X, y)
    coef = gs.best_estimator_.named_steps["m"].coef_.ravel()
    sel = np.where(np.abs(coef) > 1e-8)[0]
    if len(sel) == 0:
        sel = np.arange(X.shape[1])
    best_C = gs.best_params_["m__C"]
    model = Pipeline([("sc", StandardScaler()),
                      ("m", LogisticRegression(penalty="l1", solver="liblinear",
                                               C=best_C, max_iter=5000, random_state=seed))])
    model.fit(X[:, sel], y)
    return sel, model, {"m__C": best_C}


def fit_rf(X, y, seed):
    m = RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=1)
    m.fit(X, y)
    k = max(2, int(0.3 * X.shape[1]))
    sel = np.argsort(-m.feature_importances_)[:k]
    model = RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=1)
    model.fit(X[:, sel], y)
    return sel, model, None


def fit_svmrfe(X, y, seed):
    k = max(2, int(0.3 * X.shape[1]))
    est = SVC(kernel="linear", random_state=seed)
    sel = RFE(est, n_features_to_select=k, step=1)
    sel.fit(X, y)
    idx = np.where(sel.support_)[0]
    model = SVC(kernel="linear", probability=True, random_state=seed)
    model.fit(X[:, idx], y)
    return idx, model, None


def fit_xgb(X, y, seed):
    from xgboost import XGBClassifier
    m = XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05,
                      eval_metric="logloss", random_state=seed, n_jobs=1,
                      use_label_encoder=False)
    m.fit(X, y)
    k = max(2, int(0.3 * X.shape[1]))
    sel = np.argsort(-m.feature_importances_)[:k]
    model = XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05,
                          eval_metric="logloss", random_state=seed, n_jobs=1,
                          use_label_encoder=False)
    model.fit(X[:, sel], y)
    return sel, model, None


FITTERS = {"LASSO": fit_lasso, "RF": fit_rf, "SVM-RFE": fit_svmrfe, "XGBoost": fit_xgb}


def score_of(model, X):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X)[:, 1]
    return model.decision_function(X)


def run_task(task_name, train_gse, train_label_col, train_pos, train_neg, ext_cohorts, genes):
    print(f"\n===== {task_name} =====", flush=True)
    tr_expr, tr_ann = load_cohort(train_gse)
    Xtr, ytr, tr_genes, tr_samples = make_xy(tr_expr, tr_ann, genes, train_label_col,
                                             train_pos, train_neg)
    print(f"训练: {train_gse} n={len(ytr)} (pos={int(ytr.sum())}, neg={int((ytr==0).sum())}), 特征={len(tr_genes)}", flush=True)

    results = []
    # ---- 1) 训练集内重复分层 5 折 out-of-fold（特征选择在折内） ----
    for algo_name, fitter in FITTERS.items():
        oof_aucs, best_Cs = [], []
        for seed in SEEDS:
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
            oof = np.full(len(ytr), np.nan)
            for tr_idx, te_idx in skf.split(Xtr, ytr):
                Xa, ya = Xtr[tr_idx], ytr[tr_idx]
                Xb = Xtr[te_idx]
                if len(np.unique(ya)) < 2:
                    continue
                sel, model, params = fitter(Xa, ya, seed)
                if params:
                    best_Cs.append(params.get("m__C"))
                if len(sel) == 0:
                    sel = np.arange(Xa.shape[1])
                # 用折内重训的模型预测（fitter 已返回在 Xa[:,sel] 上拟合的模型）
                oof[te_idx] = score_of(model, Xb[:, sel])
            if np.isfinite(oof).all() and np.unique(oof).size > 1:
                oof_aucs.append(roc_auc_score(ytr, oof))
        rec = {"algorithm": algo_name, "setting": "train_repeated5foldCV",
               "n_seeds": len(oof_aucs),
               "auc_mean": float(np.mean(oof_aucs)) if oof_aucs else np.nan,
               "auc_sd": float(np.std(oof_aucs)) if oof_aucs else np.nan}
        if best_Cs:
            rec["lasso_best_C_median"] = float(np.median(best_Cs))
        results.append(rec)
        msg = f"  [CV] {algo_name}: OOF AUC = {np.mean(oof_aucs):.3f} ± {np.std(oof_aucs):.3f} (n={len(oof_aucs)})" if oof_aucs else f"  [CV] {algo_name}: FAILED"
        print(msg, flush=True)

    # ---- 2) 全训练集选特征 -> 外部验证（特征对齐） ----
    sel_counts = {g: 0 for g in tr_genes}
    per_algo_sel, per_algo_model = {}, {}
    for algo_name, fitter in FITTERS.items():
        genes_by_seed = set()
        for seed in SEEDS:
            sel, _, _ = fitter(Xtr, ytr, seed)
            for i in sel:
                genes_by_seed.add(tr_genes[i])
        per_algo_sel[algo_name] = sorted(genes_by_seed)
        for g in genes_by_seed:
            sel_counts[g] += 1
        # 最终模型：seed=1 在全训练集上（所选基因 + 在该子集上拟合的模型）
        sel1, model1, _ = fitter(Xtr, ytr, 1)
        per_algo_model[algo_name] = ([tr_genes[i] for i in sel1], model1)

    consensus = sorted([g for g, c in sel_counts.items() if c >= 3])
    print(f"  共识基因（>=3/4 算法，跨 10 种子）: {len(consensus)} -> {consensus}", flush=True)

    for algo_name, (sel, model) in per_algo_model.items():
        for gse, lcol, pos, neg in ext_cohorts:
            try:
                e_expr, e_ann = load_cohort(gse)
                common = [g for g in sel if g in e_expr.index]
                if len(common) < len(sel):
                    print(f"  [EXT] {algo_name} -> {gse}: SKIPPED（模型基因 {[g for g in sel if g not in e_expr.index]} 在外部队列缺失）", flush=True)
                    results.append({"algorithm": algo_name, "setting": f"ext_{gse}",
                                    "n_seeds": 1, "auc_mean": np.nan, "note": "skipped_missing_features"})
                    continue
                Xe, ye, _, _ = make_xy(e_expr, e_ann, common, lcol, pos, neg)
                if len(np.unique(ye)) < 2:
                    continue
                score = score_of(model, Xe)
                auc, lo, hi = delong_ci(ye, score)
                results.append({"algorithm": algo_name, "setting": f"ext_{gse}", "n_seeds": 1,
                                "auc_mean": auc, "auc_sd": np.nan, "ci_low": lo, "ci_high": hi,
                                "n": len(ye), "pos": int(ye.sum()), "neg": int((ye == 0).sum()),
                                "n_features_used": len(common)})
                print(f"  [EXT] {algo_name} -> {gse}: AUC={auc:.3f} (95%CI {lo:.3f}-{hi:.3f}), n={len(ye)}, feats={len(common)}", flush=True)
            except Exception as e:
                print(f"  [EXT] {algo_name} -> {gse}: FAILED {type(e).__name__}: {e}", flush=True)

    # ---- 3) 共识基因 logistic 回归（最终模型） ----
    if consensus:
        for gse, lcol, pos, neg in ext_cohorts:
            try:
                e_expr, e_ann = load_cohort(gse)
                common = [g for g in consensus if g in e_expr.index]
                if len(common) < 2:
                    print(f"  [EXT] Consensus-LR -> {gse}: SKIPPED（共识基因在外部队列缺失）", flush=True)
                    continue
                sel_idx = [tr_genes.index(g) for g in common]
                model = Pipeline([("sc", StandardScaler()),
                                  ("m", LogisticRegression(max_iter=5000, random_state=1))])
                model.fit(Xtr[:, sel_idx], ytr)
                Xe, ye, _, e_samples = make_xy(e_expr, e_ann, common, lcol, pos, neg)
                if len(np.unique(ye)) < 2:
                    continue
                score = model.predict_proba(Xe)[:, 1]
                auc, lo, hi = delong_ci(ye, score)
                results.append({"algorithm": "Consensus-LR", "setting": f"ext_{gse}", "n_seeds": 1,
                                "auc_mean": auc, "auc_sd": np.nan, "ci_low": lo, "ci_high": hi,
                                "n": len(ye), "pos": int(ye.sum()), "neg": int((ye == 0).sum()),
                                "n_features_used": len(common)})
                print(f"  [EXT] Consensus-LR({len(common)} genes) -> {gse}: AUC={auc:.3f} (95%CI {lo:.3f}-{hi:.3f})", flush=True)
                pd.DataFrame({"sample": e_samples, "score": score, "y": ye}).to_csv(
                    f"{OUT}/{task_name}_ConsensusLR_{gse}_scores.csv", index=False)
            except Exception as e:
                print(f"  [EXT] Consensus-LR -> {gse}: FAILED {e}", flush=True)

    pd.DataFrame(results).to_csv(f"{OUT}/{task_name}_performance.csv", index=False)
    pd.DataFrame({"gene": list(sel_counts.keys()),
                  "n_algorithms": list(sel_counts.values()),
                  "consensus": [g in consensus for g in sel_counts.keys()]}).to_csv(
        f"{OUT}/{task_name}_gene_selection.csv", index=False)
    return consensus


if __name__ == "__main__":
    pool = pd.read_csv(f"{BASE}/03_结果/WGCNA/ML_candidate_pool.csv")
    genes = pool["gene"].tolist()
    print(f"候选池: {len(genes)} 基因", flush=True)

    consA = run_task(
        "TaskA_diagnosis", "GSE134431", "group2", "DFU", "DFS",
        [("GSE199939", "group2", "DFU", "control"),
         ("GSE143735", "group2", "DFU_healer", "DM_no_ulcer"),
         ("GSE143735", "group2", "DFU_nonhealer", "DM_no_ulcer")],
        genes)

    consB = run_task(
        "TaskB_healing", "GSE134431", "healing", "Healer", "Nonhealer",
        [("GSE143735", "group2", "DFU_healer", "DFU_nonhealer")],
        genes)

    print("\nML DONE", flush=True)
