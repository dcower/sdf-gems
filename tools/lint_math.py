"""
Enforce GitHub-safe math in Markdown: inline math as $`...`$ and display math as ```math blocks.

Why: GitHub processes Markdown BEFORE handing `$...$` / `$$...$$` math to MathJax (verified against the
GitHub Markdown API):
  - every backslash before ASCII punctuation is consumed: \\# -> # (MathJax error), \\% -> % (a TeX comment that
    silently deletes the rest of the formula), \\{ \\} -> { }, \\| -> |, \\, \\; \\! -> literal punctuation,
    \\\\ -> \\ (matrix/cases rows collapse);
  - emphasis and link syntax still run across math: two `y^*` in one paragraph pair into <em>, `[$x$]`
    parses as link text, so many spans are never recognised as math and show up as raw `$...$` text.
GitHub's documented fix is its Markdown-proof syntax, $`...`$ (inline) and ```math (display): the TeX inside
is taken verbatim. Two extra rules found by testing against the GitHub Markdown API:
  - inside a GFM table row, a `|` still splits the cell (and `\\|` loses its backslash) even within code spans,
    so math in table rows must not contain `|` (use \\Vert, \\vert, \\lVert...);
  - ```math blocks nested inside list items are NOT rendered when the document also contains inline math, so
    display math inside a list item is written as a display-style inline formula on its own line:
    $`\\displaystyle ...`$.

Run via UV:
  uv run tools/lint_math.py              # fail (exit 1) on any $...$ / $$...$$ math or `|` in table-row math
  uv run tools/lint_math.py --fix        # convert $...$ -> $`...`$ and $$...$$ -> ```math, fix table-row pipes
  uv run tools/lint_math.py --render     # also render every formula with MathJax 3 (GitHub's renderer) and
                                         #   KaTeX via Node.js (packages cached in ~/.cache/sdf-gems-mathcheck)
  uv run tools/lint_math.py --github     # also render all files through the GitHub Markdown API and verify every
                                         #   formula is recognised as math (unauthenticated: 60 requests/hour;
                                         #   set GITHUB_TOKEN for more)
  uv run tools/lint_math.py docs/foo.md  # limit to specific files
"""

import argparse
import glob
import html
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request
from collections import Counter, namedtuple

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# kind: "dollar" ($...$), "display" ($$...$$), "code" ($`...`$), "fence" (```math block)
Span = namedtuple("Span", "start end bstart bend kind")


def find_math(text):
    """All math in a Markdown string (ordinary code spans and non-math fenced blocks are skipped)."""
    masked = list(text)

    def mask(a, b):
        for i in range(a, b):
            if masked[i] != "\n":
                masked[i] = " "

    spans = []
    for m in re.finditer(r"^([ \t]*)(```|~~~)([^\n]*)\n(.*?)^[ \t]*\2[ \t]*$", text, flags=re.S | re.M):
        if m.group(3).strip() == "math":
            spans.append(Span(m.start(), m.end(), m.start(4), m.end(4), "fence"))
        mask(m.start(), m.end())
    for m in re.finditer(r"\$`([^`\n]+)`\$", "".join(masked)):
        spans.append(Span(m.start(), m.end(), m.start(1), m.end(1), "code"))
        mask(m.start(), m.end())
    for m in re.finditer(r"(?<!`)(`+)(?!`).*?(?<!`)\1(?!`)", "".join(masked), flags=re.S):
        mask(m.start(), m.end())
    t, i = "".join(masked), 0
    n = len(t)
    while i < n:
        if t.startswith("$$", i):
            j = t.find("$$", i + 2)
            if j < 0:
                break
            spans.append(Span(i, j + 2, i + 2, j, "display"))
            i = j + 2
            continue
        if t[i] == "$" and (i == 0 or t[i - 1] != "\\"):
            j = i + 1
            while j < n and t[j] != "\n" and not (t[j] == "$" and t[j - 1] != "\\"):
                j += 1
            if j < n and t[j] == "$" and j > i + 1:
                spans.append(Span(i, j + 1, i + 1, j, "dollar"))
                i = j + 1
                continue
        i += 1
    return sorted(spans)


def line_of(text, pos):
    return text.count("\n", 0, pos) + 1


def in_table_row(text, pos):
    ls = text.rfind("\n", 0, pos) + 1
    return text[ls:pos].lstrip().startswith("|")


def pipe_free(tex):
    """Remove every `|` from TeX (for GFM table cells): \\| -> \\Vert, | -> \\vert, \\left| -> \\left\\vert."""
    tex = re.sub(r"\\\|", r"\\Vert ", tex)
    tex = re.sub(r"(?<!\\)\|", r"\\vert ", tex)
    return re.sub(r"(\\(?:Vert|vert)) +(?=[\s}$`)]|$)", r"\1", tex)


MD_ONLY = "*[]()-+/<@?"  # backslash escapes that only mean something to Markdown, not inside TeX math
GITHUB_DISALLOWED_MACROS = {
    "operatorname", "DeclareMathOperator", "colorbox", "fcolorbox", "bbox",
    "href", "url", "class", "style", "cssId", "unicode", "htmlData",
    "newcommand", "renewcommand", "def", "let", "gdef", "edef", "xdef", "mathchoice",
}


def sanitize_github_tex(tex):
    """Replace GitHub-disallowed macros (e.g. \\operatorname -> \\mathrm) and non-ASCII characters inside TeX."""
    tex = re.sub(r"\\operatorname\*?\s*\{", r"\\mathrm{", tex)
    tex = re.sub(r"\\operatorname\*?\s+([A-Za-z]+)", r"\\mathrm{\1}", tex)
    for k, v in [("ü", "u"), ("ö", "o"), ("ä", "a"), ("é", "e"), ("è", "e"), ("ê", "e"), ("à", "a"), ("â", "a"), ("ô", "o"), ("î", "i"), ("ç", "c"), ("–", "-"), ("—", "-"), ("−", "-")]:
        tex = tex.replace(k, v)
    return tex


def strip_markdown_escapes(tex):
    """`y^\\*` / `\\[-2, 2]` were Markdown escapes; inside verbatim $`...`$ they would be invalid TeX."""
    return sanitize_github_tex(re.sub(r"\\([A-Za-z]+|.)", lambda m: m.group(1) if m.group(1) in MD_ONLY else m.group(0), tex))


def _row_code_spans(row):
    """(start, end, is_math) for code spans in one line, CommonMark-style: a backtick run is closed by the next run
    of the same length; $`...`$ is GitHub math."""
    out, i, n = [], 0, len(row)
    while i < n:
        if row[i] == "\\" and i + 1 < n:
            i += 2
            continue
        if row[i] != "`":
            i += 1
            continue
        j = i
        while j < n and row[j] == "`":
            j += 1
        run = j - i
        k = j
        close = -1
        while k < n:
            if row[k] == "`":
                m = k
                while m < n and row[m] == "`":
                    m += 1
                if m - k == run:
                    close = k
                    break
                k = m
            else:
                k += 1
        if close < 0:
            i = j
            continue
        is_math = run == 1 and i > 0 and row[i - 1] == "$" and close + 1 < n and row[close + 1] == "$"
        out.append((i, close + run, is_math))
        i = close + run
    return out


def table_code_pipes(text):
    """[(start, end, fixed)] for ordinary code spans in GFM table rows that contain an unescaped `|`
    (it splits the cell even inside backticks; GFM's documented escape is `\\|`)."""
    out = []
    for m in re.finditer(r"^[ \t]*\|.*$", text, flags=re.M):
        row, base = m.group(0), m.start()
        for a, b, is_math in _row_code_spans(row):
            span = row[a:b]
            if not is_math and re.search(r"(?<!\\)\|", span):
                out.append((base + a, base + b, re.sub(r"(?<!\\)\|", r"\\|", span)))
    return out


def convert(text):
    """Return (new_text, n_converted, problems) with all $...$ / $$...$$ math in GitHub-safe syntax."""
    # Fix typo pattern $`...$ (missing closing backtick before $) on a single line
    text, n_typo = re.subn(r"\$`([^`\n$]+)\$(?!`)", r"$`\1`$", text)
    problems, edits = [], []
    for sp in find_math(text):
        body = text[sp.bstart:sp.bend]
        ln = line_of(text, sp.start)
        if sp.kind in ("dollar", "code"):
            tex = sanitize_github_tex(body.strip()) if sp.kind == "code" else strip_markdown_escapes(body.strip())
            if in_table_row(text, sp.start) and "|" in tex:
                tex = pipe_free(tex)
            if "`" in tex:
                problems.append((ln, "math contains a backtick; rewrite manually", body))
                continue
            new = f"$`{tex}`$"
            if new != text[sp.start:sp.end]:
                edits.append((sp.start, sp.end, new))
        elif sp.kind in ("display", "fence"):
            ls = text.rfind("\n", 0, sp.start) + 1
            le = text.find("\n", sp.end)
            le = len(text) if le < 0 else le
            prefix, suffix = text[ls:sp.start], text[sp.end:le]
            clean_body = sanitize_github_tex(body.strip()) if sp.kind == "fence" else strip_markdown_escapes(body.strip())
            if sp.kind == "display" and (prefix.strip() or suffix.strip()):
                # display math sharing its line with text: GitHub renders it inline-style anyway
                tex = clean_body
                if "`" in tex:
                    problems.append((ln, "math contains a backtick; rewrite manually", body))
                    continue
                edits.append((sp.start, sp.end, f"$`{tex}`$"))
            else:
                lines = [ln_.strip() for ln_ in clean_body.split("\n") if ln_.strip()]
                if prefix:
                    # Indented (list item): GitHub does not render ```math blocks nested in lists when the
                    # document also has inline math, so use a display-style inline formula on its own line.
                    if "`" in clean_body:
                        problems.append((ln, "math contains a backtick; rewrite manually", body))
                        continue
                    rep = prefix + "$`\\displaystyle " + " ".join(lines) + "`$"
                    if rep != text[ls:le]:
                        edits.append((ls, le, rep))
                else:
                    # Join multi-line ```math blocks where lines end with \\ so GitHub never turns \\\n into \\\
                    joined = " ".join(lines) if any(l.endswith("\\\\") for l in lines) else "\n".join(lines)
                    rep = f"```math\n{joined}\n```"
                    if rep != text[ls:le]:
                        edits.append((ls, le, rep))
    for s, e, r in sorted(edits, reverse=True):
        text = text[:s] + r + text[e:]
    pipes = table_code_pipes(text)
    for s, e, r in sorted(pipes, reverse=True):
        text = text[:s] + r + text[e:]
    return text, len(edits) + len(pipes) + n_typo, problems


def lint(text):
    """Problems in already-converted text."""
    out = []
    for m in re.finditer(r"\$`([^`\n$]+)\$(?!`)", text):
        out.append((line_of(text, m.start()), "unclosed $`...$ math span (missing closing backtick before $)", m.group(0)))
    for sp in find_math(text):
        body = text[sp.bstart:sp.bend]
        ln = line_of(text, sp.start)
        if sp.kind in ("dollar", "display"):
            out.append((ln, "use $`...`$ (inline) or a ```math block (display) instead of $-delimited math", body))
        elif sp.kind == "code" and in_table_row(text, sp.start) and "|" in body:
            out.append((ln, "`|` inside math in a table row splits the cell; use \\Vert / \\vert", body))
        elif sp.kind == "fence":
            ls = text.rfind("\n", 0, sp.start) + 1
            if text[ls:sp.start].strip() != "" or len(text[ls:sp.start]) > 0:
                out.append((ln, "indented ```math block inside list item is not rendered by GitHub; use $`\\displaystyle ...`$", body[:80]))
            if re.search(r"\\\\[ \t]*\n", body):
                out.append((ln, "line ending with \\\\ inside ```math is corrupted into \\\\\\ by GitHub; join rows onto one line", body[:80]))
        for mac in re.findall(r"\\([A-Za-z]+)", body):
            if mac in GITHUB_DISALLOWED_MACROS:
                out.append((ln, f"GitHub disallows macro \\{mac}; use \\mathrm{{...}} instead", body[:80]))
        if any(ord(ch) > 127 for ch in body):
            out.append((ln, "non-ASCII character inside math; use ASCII/LaTeX accents", body[:80]))
    for s, e, r in table_code_pipes(text):
        out.append((line_of(text, s), "unescaped `|` in a code span in a table row splits the cell; write \\|", text[s:e]))
    lines = text.split("\n")
    for i in range(len(lines) - 1):
        if re.match(r"^\s*\|.*\|\s*$", lines[i]) and re.match(r"^\s*\|[\s:\-|]+\|\s*$", lines[i + 1]):
            ncol = len(re.findall(r"(?<!\\)\|", lines[i + 1]))
            j = i
            while j < len(lines) and lines[j].strip().startswith("|"):
                if j != i + 1 and len(re.findall(r"(?<!\\)\|", lines[j])) != ncol:
                    out.append((j + 1, "table row has a different number of cells than its header", lines[j][:80]))
                j += 1
    return out


# ----------------------------------------------------------------------------- GitHub oracle
def _github_render(text):
    headers = {"Accept": "application/vnd.github+json", "Content-Type": "application/json"}
    tok = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if tok:
        headers["Authorization"] = f"Bearer {tok}"
    req = urllib.request.Request("https://api.github.com/markdown",
                                 data=json.dumps({"text": text, "mode": "gfm"}).encode(), headers=headers)
    return urllib.request.urlopen(req, timeout=60).read().decode()


def github_check(files, batch_bytes=300_000):
    """{file: [(line, formula), ...]} for formulas GitHub does not render as math (files batched per request)."""
    sep, result, batches, cur, size = "SDFGEMSFILEBREAK", {}, [], [], 0
    for f in files:
        if cur and size + os.path.getsize(f) > batch_bytes:
            batches.append(cur)
            cur, size = [], 0
        cur.append(f)
        size += os.path.getsize(f)
    if cur:
        batches.append(cur)
    for b in batches:
        srcs = [open(f, encoding="utf-8").read() for f in b]
        out = _github_render("".join(f"\n\n{sep}{i}\n\n{s}" for i, s in enumerate(srcs)))
        for f, src, part in zip(b, srcs, re.split(sep + r"\d+", out)[1:]):
            got = Counter(re.sub(r"\s+", " ", html.unescape(html.unescape(m)).strip().strip("$").strip())
                          for m in re.findall(r"<math-renderer[^>]*>(.*?)</math-renderer>", part, flags=re.S))
            miss = []
            for sp in find_math(src):
                key = re.sub(r"\s+", " ", src[sp.bstart:sp.bend].strip())
                if sp.kind in ("dollar", "display"):
                    key = re.sub(r"\\([!-/:-@\[-`{-~])", r"\1", key)  # GitHub's escape stripping
                if got[key] > 0:
                    got[key] -= 1
                else:
                    miss.append((line_of(src, sp.start), src[sp.start:sp.end]))
            if miss:
                result[f] = miss
    return result


# ----------------------------------------------------------------------------- MathJax / KaTeX
CHECK_MJS = r"""
import katex from "katex";
import { mathjax } from "mathjax-full/js/mathjax.js";
import { TeX } from "mathjax-full/js/input/tex.js";
import { SVG } from "mathjax-full/js/output/svg.js";
import { liteAdaptor } from "mathjax-full/js/adaptors/liteAdaptor.js";
import { RegisterHTMLHandler } from "mathjax-full/js/handlers/html.js";
import { AllPackages } from "mathjax-full/js/input/tex/AllPackages.js";
import fs from "fs";
RegisterHTMLHandler(liteAdaptor());
const tex = new TeX({ packages: AllPackages, formatError: (j, e) => { throw e; } });
const doc = mathjax.document("", { InputJax: tex, OutputJax: new SVG() });
const bad = [];
for (const it of JSON.parse(fs.readFileSync(process.argv[2], "utf8"))) {
  for (const [name, f] of [["MathJax", () => doc.convert(it.tex, { display: it.display })],
      ["KaTeX", () => katex.renderToString(it.tex, { throwOnError: true, displayMode: it.display, strict: "ignore" })]]) {
    try { f(); } catch (e) { bad.push({ ...it, renderer: name, error: String(e.message).split("\n")[0].slice(0, 120) }); }
  }
}
console.log(JSON.stringify(bad));
"""


def render_check(files):
    cache = os.path.expanduser("~/.cache/sdf-gems-mathcheck")
    if not os.path.isdir(os.path.join(cache, "node_modules", "mathjax-full")):
        os.makedirs(cache, exist_ok=True)
        subprocess.run(["npm", "init", "-y"], cwd=cache, check=True, capture_output=True)
        subprocess.run(["npm", "install", "--silent", "mathjax-full@3", "katex"], cwd=cache, check=True)
    open(os.path.join(cache, "check.mjs"), "w").write(CHECK_MJS)
    items = []
    for f in files:
        text = open(f, encoding="utf-8").read()
        for sp in find_math(text):
            if sp.kind in ("code", "fence"):
                items.append({"file": f, "line": line_of(text, sp.start), "tex": text[sp.bstart:sp.bend],
                              "display": sp.kind == "fence"})
    json.dump(items, open(os.path.join(cache, "items.json"), "w"))
    res = subprocess.run(["node", "check.mjs", "items.json"], cwd=cache, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stderr)
        return None, len(items)
    return json.loads(res.stdout), len(items)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*")
    ap.add_argument("--fix", action="store_true", help="convert to GitHub-safe math syntax in place")
    ap.add_argument("--render", action="store_true", help="render all math with MathJax 3 + KaTeX (needs node)")
    ap.add_argument("--github", action="store_true", help="verify rendering through the GitHub Markdown API")
    a = ap.parse_args()
    os.chdir(REPO)
    files = a.files or (["README.md"] + sorted(glob.glob("docs/**/*.md", recursive=True)))
    total = 0
    if a.fix:
        nconv = 0
        for f in files:
            text = open(f, encoding="utf-8").read()
            new, n, probs = convert(text)
            nconv += n
            if new != text:
                open(f, "w", encoding="utf-8").write(new)
            for ln, msg, body in probs:
                total += 1
                print(f"{f}:{ln}: {msg}: {body[:100]}")
        print(f"Converted {nconv} formulas.")
    for f in files:
        for ln, msg, body in lint(open(f, encoding="utf-8").read()):
            total += 1
            print(f"{f}:{ln}: {msg}: {body[:100]}")
    if a.render:
        bad, n = render_check(files)
        if bad is None:
            total += 1
        else:
            for b in bad:
                total += 1
                print(f"{b['file']}:{b['line']}: {b['renderer']} error: {b['error']}")
            print(f"Rendered {n} formulas with MathJax 3 + KaTeX: {len(bad)} errors.")
    if a.github:
        try:
            miss = github_check(files)
        except urllib.error.HTTPError as e:
            print(f"GitHub API error ({e.code} {e.reason}); unauthenticated limit is 60 requests/hour. "
                  "Set GITHUB_TOKEN or retry later.")
            sys.exit(2)
        for f, items in miss.items():
            for ln, s in items:
                total += 1
                print(f"{f}:{ln}: GitHub does not render this as math: {s[:100]}")
        print(f"GitHub check: {sum(len(v) for v in miss.values())} formulas not rendered as math.")
    print(f"{total} math problem(s).")
    sys.exit(1 if total else 0)


if __name__ == "__main__":
    main()
