# GitHub 发布操作手册（中文）

本手册说明如何把 `07_复现包/` 发布为 GitHub 仓库 `yyx-4113/kpy-dfu-transcriptome`，并完成 Zenodo 归档与版本 tag。

> 目标读者：本项目作者本人。全部操作可在本机 Windows + Git Bash 下完成。

---

## 0. 发布前检查清单

| # | 检查项 | 命令 / 方法 | 通过标准 |
|---|---|---|---|
| 1 | 复现包组装完成 | `python 02_脚本/build_repro_package.py` | 无报错；文件数与上轮一致 |
| 2 | 审计通过 | `python 02_脚本/verify_manuscript_audit.py` | `未溯源 token：0 处` |
| 3 | 脚本语法 | `python -m py_compile 07_复现包/scripts/*.py` | 无输出 |
| 4 | 无残余硬编码路径 | `grep -rl "D:/2026" 07_复现包/scripts/` | 无匹配 |
| 5 | MANIFEST 校验和正确 | 见 §2.3 | `0 mismatches` |
| 6 | 无大文件误入 | `find 07_复现包 -size +50M` | 无输出 |
| 7 | 署名与学位 | 人工核对 `CITATION.cff` | 仅 `given-names`/`family-names`/ORCID，**无 MD/PhD/MS 后缀** |

> ⚠ **第 7 项特别提醒**：本项目作者仅本科学历（B.M.）。`CITATION.cff`、投稿系统 author 头衔、cover letter 中**一律不得出现 MD / PhD / MS**。若投稿系统自动带出 "Yongxin Yang, MD"，必须删除 MD。

---

## 1. 本地 Git 仓库初始化

```bash
cd "D:/2026.10/丙酮酰化与神经病理性疼痛/丙酮酰化与糖尿病足/07_复现包"

git init
git config user.name  "Yongxin Yang"
git config user.email "960856791@qq.com"
git add .
git status                    # 逐项确认没有意外文件
```

**提交前务必确认**：
- `07_复现包/.github/` 存在（workflow 文件）
- `07_复现包/.gitignore` 生效（没有 `__pycache__`、没有 `.npz`）
- 没有 `01_数据/原始数据/` 等大目录（构建脚本本就不复制它们）

---

## 2. 提交与校验

```bash
git commit -m "Reproducibility package v1.0.0 for Kpy-DFU transcriptomic footprint study"
git log --stat -1
```

### 2.3 校验 MANIFEST（发布后也要跑一次）

```bash
cd 07_复现包
python - <<'PY'
import csv, hashlib, pathlib
bad = 0
for r in csv.DictReader(open('MANIFEST.csv', encoding='utf-8')):
    p = pathlib.Path(r['path'])
    if not p.exists():
        print('MISSING', r['path']); bad += 1; continue
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    if h != r['sha256']:
        print('MISMATCH', r['path']); bad += 1
print(f'{bad} problems')
PY
```

---

## 3. 创建 GitHub 仓库并推送

```bash
# 在 GitHub 网页创建空仓库（不要勾选 README/LICENSE/.gitignore）
# 仓库名：kpy-dfu-transcriptome
# 可见性：Public（投稿需要可公开访问的 DOI 链接）

cd 07_复现包
git remote add origin https://github.com/yyx-4113/kpy-dfu-transcriptome.git
git branch -M main
git push -u origin main
```

**认证**：使用 Personal Access Token（Fine-grained，权限 `Contents: Read and write`），不要用密码。若已配置 `gh` CLI：

```bash
gh auth status
gh repo create yyx-4113/kpy-dfu-transcriptome --public --source=. --push
```

---

## 4. 版本 tag 与 release

```bash
cd 07_复现包
git tag -a v1.0.0 -m "Kpy-DFU transcriptomic footprint: reproduction package v1.0.0"
git push origin v1.0.0

gh release create v1.0.0 \
  --title "Kpy × DFU reproduction package v1.0.0" \
  --notes "Reproduction package for the Kpy-DFU transcriptomic footprint manuscript. Data availability statement should cite this release tag and the associated Zenodo DOI." \
  || git push origin v1.0.0
```

> **数据可用性声明必须引用具体版本**，不要写"master 分支"。推荐写法：
> "Analysis code, derived result tables and the audit ledger are available at
> https://github.com/yyx-4113/kpy-dfu-transcriptome (release v1.0.0,
> commit <sha>, archived at Zenodo DOI 10.5281/zenodo.<record>)."

---

## 5. Zenodo 归档（大体量文件 + 整体快照）

本仓库只含代码与小体积产物；下列内容需 Zenodo 归档：

| 内容 | 清单 |
|---|---|
| 单细胞 count 矩阵、处理后表达矩阵 | `ZENODO_ARCHIVE_LIST.csv` |
| 大型元数据参考（GENCODE gtf、MSigDB gmt、Finngen 表） | 同上 |

### 5.1 创建 Zenodo 记录

1. 登录 <https://zenodo.org> → **New upload**
2. 上传 `07_复现包/` 的压缩包（不含 `.git`）：

```bash
cd "D:/2026.10/丙酮酰化与神经病理性疼痛/丙酮酰化与糖尿病足"
tar -czf kpy_dfu_repro_v1.0.0.tar.gz --exclude='.git' 07_复现包/
```

3. 元数据填写：
   - **Title**: Reproducibility package — Transcriptomic footprint of the lysine-pyruvylation regulatory network in diabetic foot ulcer tissue
   - **Upload type**: Software
   - **Keywords**: lysine-pyruvylation; diabetic foot ulcer; transcriptomics; meta-analysis; single-cell RNA-seq; Mendelian randomisation; molecular docking
   - **Creators**: Yang, Yongxin (ORCID 0009-0004-9698-6552) — **不要填写学位后缀**
   - **License**: MIT
   - **Related identifiers**: 填入 GitHub 仓库 URL（`Related identifier type: Software`）
4. 发布后记录 DOI（形如 `10.5281/zenodo.XXXXXXX`）与 **concept DOI** 与 **version DOI**
5. 把 version DOI 写入稿件的 Data Availability 段落

### 5.2 校验 Zenodo 上传完整性

```bash
# 对照 Zenodo 记录中的文件列表与本地清单
python - <<'PY'
import csv
rows = list(csv.DictReader(open('07_复现包/ZENODO_ARCHIVE_LIST.csv', encoding='utf-8')))
print(f'{len(rows)} files to archive, total '
      f'{sum(int(r["bytes"]) for r in rows)/1024/1024:.1f} MB')
PY
```

---

## 6. 投稿时需要同步更新的三处

| 位置 | 改什么 |
|---|---|
| 稿件 Data Availability 段 | GitHub URL + release tag + commit sha + Zenodo version DOI |
| `CITATION.cff` 的 `repository-code` / `url` | 确认指向真实仓库地址；`version` 与 release tag 一致 |
| `07_复现包/README.md` §7 | 填入 Zenodo DOI 与版本号 |

> 纪律：**不要在 Zenodo 记录创建之前就把 DOI 写进稿件并投稿**。必须先取得真实 DOI 并回读确认可访问，再写入稿件。

---

## 7. 修改后如何同步（增量）

```bash
cd 07_复现包
# 1) 重新组装（会覆盖旧文件并刷新 MANIFEST）
python ../02_脚本/build_repro_package.py
# 2) 复跑审计
python ../02_脚本/verify_manuscript_audit.py
# 3) 提交
git add -A
git commit -m "Update results, audit ledger and supplementary tables"
git push origin main
```

若正文数字有改动，务必按**修上游 bug → 全下游一致性扫描**的流程走：grep 全部旧数字指纹 → 重跑 `build_audit_facts.py` → 重跑 `verify_manuscript_audit.py` → 确认「未溯源 token = 0」→ 再提交。

---

## 8. 常见问题

| 现象 | 处理 |
|---|---|
| `git push` 报认证失败 | 用 Fine-grained PAT；确认 token 有 `Contents: Read and write` |
| 仓库 > 100 MB 被拒 | 检查 `.gitignore` 是否生效；`07_复现包` 本身约 32 MB，不应超限 |
| Zenodo 上传 502 / 超时 | 分卷上传或用 Zenodo 命令行客户端 `zenodo upload` |
| 中文路径导致脚本写入失败 | R 脚本开头已含 `Sys.setlocale("LC_ALL", "Chinese_China.65001")`；若仍失败，改为在 ASCII 路径下运行并管道输出 |
| BLAST v4 查询失败 | NCBI 限流；脚本已内置 `time.sleep(4)` 与失败重试，等待后重跑（幂等） |
