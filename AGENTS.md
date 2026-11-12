# sdf-gems Workspace Rules

## Privacy & Attribution
- **Never include personal names or local usernames** in any code comment, documentation file (`.md`), rule/skill file, image, or git commit message/author trailer.
- Attribute AI-discovered primitives, operators, and optimizations to the model used (e.g., **Gemini 4 Argon**) in their respective `docs/<category>/<name>.md` files, and avoid naming people in the top-level `README.md` tables.

## Git & Repository Hygiene
- **One logical unit per commit**: When adding or upgrading primitives/operators in batches, commit each primitive or operator separately with clear before/after metrics.
- **Keep the repository root clean**:
  - Never leave temporary test scripts, scratch files, or one-off debug harnesses in the repository root. Put temporary scripts in an external scratch directory outside the repository.
  - Place reusable rendering and verification scripts in `tools/`.
- **Single-Translation-Unit Compatibility**: Every `.glsl` header in `glsl/` must compile cleanly when included together via `sdf_gems.glsl` in any order. Guard shared helpers or cross-included functions with `#ifndef SDF_GEMS_HAS_<NAME>` macros.

## Documentation & README Formatting
- **Modular Structure**: Every SDF or operator belongs in its own `glsl/<category>/<name>.glsl` header (using the consolidated top-level categories: `primitives3d`, `polyhedra`, `curves3d`, `primitives2d`, `booleans`, `gradients3d`, `deform`, `repetition`, `volume`, `optics`, `duality`, and `common`—avoid creating new one-off singleton category directories) with a matching `docs/<category>/<name>.md` explanation and rendered previews in `docs/images/<category>/`.
- **GitHub-Safe Math Only**: Write inline math as `` $`...`$ `` and display math as a ```` ```math ```` block (inside list items, use `` $`\displaystyle ...`$ `` on its own line). Never use bare `$...`` $` / `$ ``$...$$`: GitHub eats backslash escapes (`\%`, `\{`, `\|`, `\\`, `\,` ...) and runs emphasis/link syntax across them. Inside table rows, math and code spans must not contain an unescaped `|`. `uv run tools/update_umbrella.py` enforces this; `uv run tools/lint_math.py --fix --render --github` converts, renders with MathJax + KaTeX, and verifies against GitHub's renderer.
- **Clean Root `README.md` & Per-Category Galleries (`docs/<category>/README.md`)**: Keep the root `README.md` concise (category showcase table + collapsed quick-search index), and place individual 3-column visual rows (`Preview | Function | Summary`) in their respective `docs/<category>/README.md` gallery so pages stay fast and never scroll horizontally.

## Mandatory Visual & Numerical Verification
- **Inspect Every Rendered Image**: Never commit a generated `.png` without viewing it via `read_file` first. Verify that:
  1. 3D previews (`docs/images/<category>/<name>.png`) are `360x270` 3D raymarched ground-plane renders (not 2D slices or wide multi-panel banners) with the object properly framed and unclipped.
  2. 2D cross-section slices (`*_slice.png`) clearly show interior and exterior level sets without `NaN` holes, branch discontinuities, or blobby distortions.
- **Test the Literal `.glsl` Code**: Never rely solely on a Python prototype to validate a GLSL shader. Always compile and numerically test the actual `.glsl` file (e.g., via a C++17 GLSL-compatibility test harness) against brute-force ground truth before committing.
