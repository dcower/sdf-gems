"""
Maintenance & audit script for `sdf-gems`.

1. Regenerates `sdf_gems.glsl` with `#include` directives for every `.glsl` file under `glsl/`
   (`glsl/common/common.glsl` first, followed by all other categories in sorted order).
2. Verifies that every `.glsl` file (except `common.glsl`) has a matching `.md` documentation file
   under `docs/<category>/<name>.md`.
3. Audits all Markdown image references in `README.md` and `docs/**/*.md` to verify that every
   referenced PNG exists on disk.
4. Lints all Markdown math for GitHub-safe syntax via `tools/lint_math.py` (inline $`...`$, display ```math;
   no `|` inside table-row math or table-row code spans).

Run via UV:
  uv run tools/update_umbrella.py
"""

import glob
import os
import re
import sys


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(repo_root)

    glsl_files = sorted(glob.glob("glsl/**/*.glsl", recursive=True))
    common_file = "glsl/common/common.glsl"
    ordered_glsl = [common_file] + [f for f in glsl_files if f != common_file]

    lines = [
        "// sdf-gems: Auto-generated single-include umbrella header",
        "#ifndef SDF_GEMS_GLSL",
        "#define SDF_GEMS_GLSL",
        "",
    ]
    for path in ordered_glsl:
        lines.append(f'#include "{path}"')
    lines.extend(["", "#endif // SDF_GEMS_GLSL", ""])

    with open("sdf_gems.glsl", "w") as f:
        f.write("\n".join(lines))

    missing_docs = []
    missing_catalog_links = []
    missing_root_links = []
    with open("README.md") as f:
        readme_txt = f.read()
    cat_readmes = {}
    for cr in sorted(glob.glob("docs/*/README.md")):
        cat = cr.split(os.sep)[1]
        with open(cr) as f:
            cat_readmes[cat] = f.read()
        if f"(docs/{cat}/README.md)" not in readme_txt:
            missing_root_links.append(cr)

    for g in ordered_glsl:
        if g == common_file:
            continue
        parts = g.split(os.sep)
        cat = parts[1]
        stem = os.path.splitext(parts[2])[0]
        doc_path = re.sub(r"^glsl/", "docs/", re.sub(r"\.glsl$", ".md", g))
        if not os.path.exists(doc_path):
            missing_docs.append((g, doc_path))
        cat_txt = cat_readmes.get(cat, "")
        if f"({stem}.md)" not in cat_txt or f"(../../{g})" not in cat_txt:
            missing_catalog_links.append(g)
        if f"({doc_path})" not in readme_txt and f"({g})" not in readme_txt:
            missing_root_links.append(g)

    md_files = ["README.md"] + sorted(glob.glob("docs/**/*.md", recursive=True))
    missing_images = []
    for md in md_files:
        base = os.path.dirname(md)
        with open(md) as f:
            txt = f.read()
        for m in re.findall(r'(?:!\[.*?\]\((.*?)\)|src="(.*?)")', txt):
            rel_img = m[0] or m[1]
            full_img = os.path.normpath(os.path.join(base, rel_img))
            if not os.path.exists(full_img):
                missing_images.append((md, rel_img, full_img))

    sys.path.insert(0, os.path.join(repo_root, "tools"))
    from lint_math import lint as lint_math

    math_problems = []
    for md in md_files:
        with open(md) as f:
            math_problems += [(md, ln, msg) for ln, msg, _ in lint_math(f.read())]

    print(f"Updated sdf_gems.glsl with {len(ordered_glsl)} GLSL headers.")
    print(f"Audited {len(md_files)} Markdown files.")
    if missing_docs:
        print("ERROR: Missing documentation files:", missing_docs)
        sys.exit(1)
    if missing_catalog_links:
        print("ERROR: GLSL headers missing from docs/<category>/README.md tables:", missing_catalog_links)
        sys.exit(1)
    if missing_root_links:
        print("ERROR: Entries missing from root README.md:", missing_root_links)
        sys.exit(1)
    if missing_images:
        print("ERROR: Missing image files:", missing_images)
        sys.exit(1)
    if math_problems:
        for md, ln, msg in math_problems[:50]:
            print(f"ERROR: {md}:{ln}: {msg}")
        print(f"ERROR: {len(math_problems)} math problem(s); run `uv run tools/lint_math.py --fix`.")
        sys.exit(1)
    print("Audit passed: 0 missing docs, 0 broken image links, 0 math problems.")


if __name__ == "__main__":
    main()
