# 预注册笔记 — 基因集与分析规则锁定

锁定时间：2026-10-29（正式分析运行之前）。锁定后不得更换基因集、阈值或分组；如需偏离，必须作为单独的敏感性分析并明示。

## 1. 预定义基因集（分析前锁定）

| 基因集 | 基因 | 依据 |
|---|---|---|
| SMC_CONTRACTILE | MYH11, ACTA2, TAGLN, CNN1, MYL9, MYLK, LMOD1, RGS5 | 经典收缩表型标记（试点已用前 5 个） |
| ECM_REMODELING | MMP9, TIMP1, COL1A1, POSTN, MMP2, COL3A1, FN1, LOX | 原文报告的 ECM 重构轴 [105, 108] |
| INFLAMMATION | IL6, CXCL2, CCL2, VCAM1, ICAM1, CXCL8, CCL20 | 原文报告的炎症轴 [105, 108] |
| MS_CHANNEL | PIEZO1, PIEZO2, TMC4, TMC5, TMC6, TMC8, TMEM63A, TMEM63B | GO:0008381，第一轮审计锁定版本，原样沿用 |
| OSTEogenic | RUNX2, SOX9, SPP1, BGLAP, BMP2 | CHUV 原文的成骨转换轴 [108]（作为其结论的复制检验） |

阳性对照：SMC_CONTRACTILE 必须在两队列显著下调（原文结论 [105, 108]）；若不成立，先停止并检查流程。

## 2. 分析规则（锁定）

1. 主模型：limma（log2CPM + limma-trend 或 voom），设计 `~ patient + time`，配对。
2. 敏感性：edgeR QL、DESeq2，方向一致性为硬标准。
3. 显著性：BH FDR < 0.05；效应量报告 logFC + 95% CI。
4. 跨队列一致性：两队列 FDR<0.05 的基因并集中方向一致率；二项检验（H0: p=0.5）；预注册通过阈值 >80%，止损阈值 ≤60%。
5. 基因集检验：camera（配对设计存于 model 矩阵）+ fgsea（按 logFC 排序），双方法均报告。
6. 空间：Visium 斑点聚合到患者×切片；GeoMx 患者阻断；只做定位，不做 spot 级 p 值包装。
7. GSE241205：FPKM log2 后配对 limma，仅预定义签项方向评分（探索性）。
8. PIEZO1 次要分析：单基因跨队列方向 + MS_CHANNEL 集 + Visium 斑点定位，一次性执行。

## 3. 数据版本

- GSE310013_raw_counts.txt.gz（下载于 2026-10-29）
- GSE326710_Raw_counts_all_samples.csv.gz（同上）
- GSE299685_Normalized_Counts.xlsx（同上）
- GSE310351：6/6 样本 matrix/features/barcodes/positions（2026-10-29）
- GSE241205_gene_fpkm.txt.gz（同上）
