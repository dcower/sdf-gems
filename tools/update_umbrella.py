"""
Maintenance & audit script for `sdf-gems`.

1. Regenerates `sdf_gems.glsl` with `#include` directives for every `.glsl` file under `glsl/`
   (`glsl/common/common.glsl` first, followed by all other categories in sorted order).
2. Verifies that every `.glsl` file (except `common.glsl`) has a matching `.md` documentation file
   under `docs/<category>/<name>.md`.
3. Audits all Markdown image references in `README.md` and `docs/**/*.md` to verify that every
   referenced PNG exists on disk.

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
    for g in ordered_glsl:
        if g == common_file:
            continue
        doc_path = re.sub(r"^glsl/", "docs/", re.sub(r"\.glsl$", ".md", g))
        if not os.path.exists(doc_path):
            missing_docs.append((g, doc_path))

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

    print(f"Updated sdf_gems.glsl with {len(ordered_glsl)} GLSL headers.")
    print(f"Audited {len(md_files)} Markdown files.")
    if missing_docs:
        print("ERROR: Missing documentation files:", missing_docs)
        sys.exit(1)
    if missing_images:
        print("ERROR: Missing image files:", missing_images)
        sys.exit(1)
    print("Audit passed: 0 missing docs, 0 broken image links.")


if __name__ == "__main__":
    main()
