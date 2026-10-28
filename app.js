// Shared helpers: nav, bilingual toggle (EN default), CSV parsing, hypergeometric test.
const I18N = {
  nav_home: {en: "Home", zh: "首页"},
  nav_genes: {en: "Gene cards", zh: "基因证据卡"},
  nav_checker: {en: "List checker", zh: "列表检查器"},
  nav_regs: {en: "Regulators", zh: "调控子"},
  nav_dl: {en: "Downloads", zh: "下载"},
  lang_btn: {en: "中文", zh: "English"},
  site_sub: {en: "A shared epithelial–immune transcriptional program in nasal polyps and inflammatory bowel disease — queryable evidence resource",
             zh: "鼻息肉与炎症性肠病共享上皮–免疫转录程序 — 可查询证据资源"},
  footer: {en: "Source data: NP×IBD multi-omics analysis (GEO + PRIDE PXD050335). DEG rule: FDR < 0.05 and |log2FC| > 0.585. Scores: zmean default (reproduces Fig. 10 and the dupilumab paired analysis); ssGSEA available in score_program.R (reproduces bulk validation).",
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
