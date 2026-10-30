// Shared helpers: nav, bilingual toggle (EN default), CSV parsing, hypergeometric test.
const I18N = {
  nav_home: {en: "Home", zh: "首页"},
  nav_genes: {en: "Gene cards", zh: "基因证据卡"},
  nav_checker: {en: "List checker", zh: "列表检查器"},
  nav_regs: {en: "Regulators", zh: "调控子"},
  nav_dl: {en: "Downloads", zh: "下载"},
  dl_file: {en: "File", zh: "文件"},
  dl_contents: {en: "Contents", zh: "内容说明"},
  dl_r1: {en: "Per-gene evidence matrix, 590 genes × all cohorts (bulk, single-cell, spatial, protein, dupilumab)", zh: "逐基因证据矩阵，590 个基因 × 全部队列（bulk、单细胞、空间、蛋白、dupilumab）"},
  dl_r2: {en: "Program gene sets: NP×UC UP 473 / DOWN 78, NP×CD UP 174 / DOWN 52, core UP/DOWN (187 total)", zh: "程序基因集：NP×UC 上调 473 / 下调 78，NP×CD 上调 174 / 下调 52，核心上调/下调（共 187）"},
  dl_r3: {en: "SCENIC regulons, 278 regulators with target lists", zh: "SCENIC 调控子，278 个调控因子及靶基因列表"},
  dl_r4: {en: "Regulon × program-arm overlap statistics (1,112 rows)", zh: "调控子 × 程序臂重叠统计（1,112 行）"},
  dl_r5: {en: "R scoring function (locked manuscript formula; zmean default, ssGSEA option; reproduces Fig. 12 and the dupilumab paired pattern — see acceptance tests in the repository)", zh: "R 打分函数（稿件固定公式；默认 zmean，可选 ssGSEA；复现图 12 与 dupilumab 配对模式——验收测试见仓库）"},
  dl_r6: {en: "Matrix documentation, sources, and honest limitations", zh: "矩阵文档、数据来源与如实标注的局限"},
  dl_r7: {en: "Build QC log with input file md5s", zh: "构建质控日志（含输入文件 md5）"},
  dl_zenodo_h: {en: "Zenodo record:", zh: "Zenodo 记录："},
  dl_zenodo_a: {en: " the evidence matrix is published at ", zh: " 证据矩阵已发布于 "},
  dl_zenodo_b: {en: " (CC-BY-4.0) and cited in the manuscript Data Availability section. Gene sets and the scoring function are in the table above and in the GitHub repository.", zh: "（CC-BY-4.0），并在稿件 Data Availability 部分引用。基因集与打分函数见上表及 GitHub 仓库。"},
  dl_small: {en: "All files are derived from the published analysis; no new statistics were computed for this site. Program scores use the locked formula (UP mean z minus DOWN mean z for zmean; GSVA ssGSEA for bulk validation).", zh: "本站全部文件均源自稿件分析，未为建站计算新统计量。程序分数使用固定公式（zmean 为上调臂均值 z 减下调臂均值 z；bulk 验证用 GSVA ssGSEA）。"},
  lang_btn: {en: "中文", zh: "English"},
  site_sub: {en: "A shared epithelial–immune transcriptional program in nasal polyps and inflammatory bowel disease — queryable evidence resource",
             zh: "鼻息肉与炎症性肠病共享上皮–免疫转录程序 — 可查询证据资源"},
  footer: {en: "Source data: NP×IBD multi-omics analysis (GEO + PRIDE PXD050335). DEG rule: FDR < 0.05 and |log2FC| > 0.585. Scores: zmean default (reproduces Fig. 12 and the dupilumab paired analysis); ssGSEA available in score_program.R (reproduces bulk validation).",
           zh: "数据来源：NP×IBD 多组学分析（GEO＋PRIDE PXD050335）。DEG 标准：FDR＜0.05 且 |log2FC|＞0.585。打分：默认 zmean（复现 Fig.10 与 dupilumab 配对分析）；ssGSEA 见 score_program.R（复现 bulk 验证）。"}
};
function lang() { return localStorage.getItem("npibd_lang") || "en"; }
function setLang(l) { localStorage.setItem("npibd_lang", l); applyLang(); }
function applyLang() {
  const l = lang();
  document.querySelectorAll("[data-i18n]").forEach(el => {
    const k = el.getAttribute("data-i18n");
    if (I18N[k]) el.textContent = I18N[k][l];
  });
  const b = document.getElementById("langBtn");
  if (b) b.textContent = I18N.lang_btn[l];
}
function nav() {
  const l = lang();
  const links = [["index.html", "nav_home"], ["genes.html", "nav_genes"], ["checker.html", "nav_checker"], ["regulators.html", "nav_regs"], ["download.html", "nav_dl"]];
  const n = document.createElement("nav");
  links.forEach(([href, k]) => {
    const a = document.createElement("a"); a.href = href; a.setAttribute("data-i18n", k); a.textContent = I18N[k][l]; n.appendChild(a);
  });
  const b = document.createElement("button"); b.id = "langBtn";
  b.onclick = () => setLang(lang() === "en" ? "zh" : "en");
  b.textContent = I18N.lang_btn[l]; n.appendChild(b);
  document.body.insertBefore(n, document.body.firstChild.nextSibling);
}
function parseCSV(text) {
  const rows = []; let row = [], cur = "", q = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (q) { if (c === '"') { if (text[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
    else if (c === '"') q = true;
    else if (c === ",") { row.push(cur); cur = ""; }
    else if (c === "\n") { row.push(cur); rows.push(row); row = []; cur = ""; }
    else if (c !== "\r") cur += c;
  }
  if (cur !== "" || row.length) { row.push(cur); rows.push(row); }
  const head = rows[0]; const out = [];
  for (let r = 1; r < rows.length; r++) {
    if (rows[r].length < head.length) continue;
    const o = {}; head.forEach((h, j) => o[h] = rows[r][j]); out.push(o);
  }
  return out;
}
// log C(n,k) via summed logs; hypergeometric survival P(X >= k)
function logC(n, k) {
  if (k < 0 || k > n) return -Infinity;
  k = Math.min(k, n - k); let s = 0;
  for (let i = 1; i <= k; i++) s += Math.log(n - k + i) - Math.log(i);
  return s;
}
function hyperP(N, K, n, k) {
  const denom = logC(N, n); let mx = -Infinity; const terms = [];
  for (let i = k; i <= Math.min(K, n); i++) {
    const t = logC(K, i) + logC(N - K, n - i) - denom;
    terms.push(t); if (t > mx) mx = t;
  }
  let s = 0; terms.forEach(t => s += Math.exp(t - mx));
  return Math.min(1, Math.exp(mx) * s);
}
function fmtP(p) { return p < 1e-300 ? "<1e-300" : p.toExponential(2); }
document.addEventListener("DOMContentLoaded", () => { nav(); applyLang(); });
