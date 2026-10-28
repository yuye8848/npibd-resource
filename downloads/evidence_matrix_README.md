# NP×IBD 逐基因证据矩阵 — 说明

构建脚本：`build_evidence_matrix.R`（base R，无包依赖；输入 md5 见 `evidence_matrix_qc.txt`）。
构建日期：2026-10-28。行数：590 基因（唯一 SYMBOL，无重复行；脚本内置 `stopifnot(nrow==590)` 防合并膨胀）。

## 基因宇宙

- 590 ＝ NP×UC 551 ∪ NP×CD 226（方向一致交集，稿件发现集）。
- 核心 187 ＝ 551 ∩ 226（已验证：核心文件与交集完全一致）。
- 551 文件名中的"921"为历史遗留命名，其内容实为 187 核心基因（188 行含表头），不要按文件名理解。

## 列定义

| 列 | 含义 | 来源文件 |
|---|---|---|
| SYMBOL | 基因符号 | 各来源 join 键 |
| in_NPxUC551 / in_NPxCD226 / in_core187 | 程序归属 0/1 | pilot/shared_NPxUC.csv、shared_NPxCD.csv、shared_NPxIBD921.csv |
| dir_NP / dir_UC / dir_CD | 方向（1 上调、-1 下调） | 同上共享文件（dir_IBD 在 UC 文件中为 UC 方向，在 CD 文件中为 CD 方向） |
| NP_logFC / NP_padj | 鼻息肉发现集差异（GSE338426，14 NP vs 8 对照） | pilot/NP_DEGs.csv |
| GSE198950_logFC / GSE198950_padj | 鼻侧独立验证队列差异（仅 DEG 者有值；23/590） | analysis/01_validation/GSE198950_DEGs.csv |
| UC_logFC / UC_adjP / UC_isDEG | UC 结肠差异（GSE75214；多探针对应同一基因时取 \|logFC\| 最大探针，与稿件 pilot 方法一致） | sci-ar-ibd/pilot/IBD_allgenes.csv、IBD_DEGs.csv（稿件 pilot 本人所用输入） |
| CD_logFC / CD_adjP / CD_isDEG | CD 回肠差异（同上方法） | sci-ar-ibd/pilot/CDileum_allgenes.csv、CDileum_DEGs.csv |
| sc_epi_diff / sc_epi_p | 鼻侧单细胞上皮供体级 pseudobulk 差异（GSE235711） | analysis/06_revision/pseudobulk_gene_tests_epithelial.csv |
| sc_imm_diff / sc_imm_p | 鼻侧单细胞免疫供体级 pseudobulk 差异 | analysis/06_revision/pseudobulk_gene_tests_immune.csv |
| gut_sc_epi_diff / gut_sc_epi_p | 肠侧单细胞上皮供体级差异（GSE116222，UC n=6 vs 健康 n=3，即稿件肠侧单细胞比较） | analysis/06_revision/GSE116222_epithelial_gene_tests.csv |
| gut_sc_epi_infl_diff / gut_sc_epi_infl_p | 同上文件的炎症亚组敏感性比较（UC 炎症 n=3 vs 健康 n=3） | 同上 |
| spatial_CD45_logFC / spatial_CD45_adjP | GeoMx CD45 区室疾病差异 | analysis/03_spatial/geomx_disease_DEGs_CD45.csv |
| spatial_PanCK_logFC / spatial_PanCK_adjP | GeoMx PanCK 区室疾病差异 | analysis/03_spatial/geomx_disease_DEGs_PanCK.csv |
| spatial_compartment_logFC / spatial_compartment_FDR | GeoMx CD45 vs PanCK 区室对比 | analysis/03_spatial/geomx_compartment_DEGs.csv |
| protein_logFC / protein_p / protein_padj | 肠侧蛋白组差异（PXD050335 结肠黏膜，炎症 vs 非炎症） | analysis/08_proteomics/PXD050335_mucosa_IvsN_protein_stats.csv |
| dupi_logFC / dupi_FDR | dupilumab 治疗后 vs 治疗前逐基因差异（GSE267331） | sci-crswnp-resource/analysis/step2_bulk/de/GSE267331_dupilumab_post_vs_pre.csv（与稿件同一队列／同一对比的逐基因 companion 计算；稿件正文报的是程序水平配对检验） |
| placebo_logFC / placebo_FDR | 安慰剂后 vs 前逐基因差异（阴性对照） | sci-crswnp-resource/analysis/step2_bulk/de/GSE267331_placebo_post_vs_pre.csv |
| badge_bulk / badge_sc / badge_spatial / badge_protein | 四层证据徽章（见下） | 由各层数据存在性导出 |
| n_layers_4 | 四层徽章之和（1–4） | 导出 |
| dupi_sig | dupilumab FDR<0.05 标志 | 导出 |

## 徽章规则（按 SPEC 第 1 节：徽章＝该层有数据，不是显著性）

- badge_bulk：恒为 1（矩阵宇宙本身就是方向一致共享 DEG，属构造性事实）。
- badge_sc：上皮或免疫任一区室有可检验的供体级数据（p 非 NA）。注意：两区室各有 7／9 个基因因两组表达全为 0 而 p＝NA，属"无数据"，不是"不显著"，在矩阵中以 NA 呈现。
- badge_spatial：CD45 或 PanCK 任一疾病对比有数据（169 基因在 GeoMx panel 内）。
- badge_protein：PXD050335 蛋白统计表中有该基因（112 基因）。
- 显著性一律看各层 p／FDR 列，不由徽章代替。

## QC 实测（2026-10-28 构建运行）

- 行数 590；UC 551、CD 226、核心 187；核心＝交集（TRUE）。
- 共享文件 NP logFC 与 NP_DEGs 最大绝对差 0（n＝551，CD 独有基因无共享文件 logFC，以 NP_DEGs 为准）。
- 方向一致性：NP 符号＝dir_NP 590／590；UC 成员 551／551；CD 成员 226／226。
- dupilumab 逐基因均值：UP 臂 −0.0244（n＝480）vs DOWN 臂 ＋0.1508（n＝98），与稿件程序水平逆转结论方向一致；安慰剂呈反向（UP ＋0.1169、DOWN −0.0738），符合阴性对照。
- 覆盖（tested／significant＜0.05）：鼻上皮 sc 93／34；鼻免疫 sc 469／38；肠上皮 sc 30／5；肠炎症亚组 30／0；GeoMx CD45 169／0；GeoMx PanCK 169／0；GeoMx 区室对比 169／154；蛋白 112／13；dupilumab 578／88；GSE198950 23／23；UC 590／566；CD 590／301。
- 四层分布：1 层 28 基因、2 层 315、3 层 217、4 层 30。

## 队列级（非逐基因）证据——不在矩阵列中，勿误读为逐基因

- GSE193677（IBD n＝2,490）：稿件验证为程序分数水平（12／12 显著），无逐基因表，故矩阵无此列。
- GSE198950：程序水平 GSEA 名义显著（NES 1.26／−1.42），逐基因仅 23／590 为 DEG，与稿件"inconclusive"定性一致。
- dupilumab 程序水平配对检验（DOWN 臂逆转 p＝0.0023／0.0041）见稿件正文与 paired_tests_GSE267331.csv；矩阵 dupi 列为逐基因 companion 计算。

## 诚实说明与限制

- GeoMx 疾病对比逐基因在 FDR＜0.05 下为 0／169（40 ROI 逐基因功效低）；稿件的空间结论是程序水平（区室差 1.066，p＝1.2×10⁻⁴）＋区室对比 154／169 显著，不要把矩阵的 0 理解为"空间无证据"。
- 肠侧单细胞层很薄：上皮 30／590、炎症亚组 0 显著，n＝6 vs 3；这是 GSE116222 本身的局限，矩阵如实呈现，不放大。
- GSE198950 仅 23／590 为 DEG（16 CRS vs 4 对照，欠功效），与稿件 inconclusive 一致。
- 构建中曾出现一次合并膨胀（肠侧文件含两个比较致行数 590→622），被 QC 发现后修复（两个比较拆为独立列对＋stopifnot 守卫），最终版行数 590 已验证；中间错误版未用于任何下游。
