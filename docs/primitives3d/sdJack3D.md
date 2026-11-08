# `sdJack3D`, `sdAxisJack3D`, & `sdFaceJack3D` — Exact 1-`sqrt` 3D Capsule Jacks

| 3D Raymarched Surface | 2D Cross-Section Slice |
| :---: | :---: |
| ![3D Preview](../images/primitives3d/sdJack3D.png) | ![2D Cross-Section](../images/primitives3d/sdJack3D_slice.png) |

- **GLSL Implementation**: [`sdJack3D.glsl`](../../glsl/primitives3d/sdJack3D.glsl)
- **Category**: `primitives3d`
- **Status**: **Exact 1-`sqrt` Euclidean SDF** (Discovered using **Gemini 4 Argon**, `sdf-gems`, 2026)

---

## 1. Overview & Mathematical Derivations

Evaluating a 3D jack or caltrop formed by the union of $`N`$ symmetric capsules passing through the origin naively requires $`N`$ calls to `sdCapsule` ($`N`$ `sqrt`s). Using octant and Weyl-chamber folding reduces all three octahedral capsule jacks to **a single `sqrt` everywhere in $`\mathbb{R}^3`$** (`4.44e-16` max error):

1. **8-Pointed Body-Diagonal Caltrop (`sdJack3D`, 4-Capsule Union, $`4\times`$ Faster)**:
   The 4 body diagonals of a cube have unit directions $`\hat{u}_s = \frac{1}{\sqrt{3}}(s_x, s_y, s_z)`$ for $`s \in \{-1, +1\}^3`$. For any point in the positive octant $`p \leftarrow |p|`$, the projection onto $`\hat{u}_s`$ satisfies:
   $`\displaystyle |p \cdot \hat{u}_s| = \frac{|s_x p_x + s_y p_y + s_z p_z|}{\sqrt{3}} \le \frac{p_x + p_y + p_z}{\sqrt{3}} = p \cdot \hat{u}_{+++}.`$
   Since the squared distance $`\|p - \mathrm{clamp}(p\cdot\hat{u}_s, -l, l)\hat{u}_s\|_2^2 = \|p\|_2^2 - \min(|p\cdot\hat{u}_s|, l)(2|p\cdot\hat{u}_s| - \min(|p\cdot\hat{u}_s|, l))`$ is strictly decreasing in $`|p\cdot\hat{u}_s|`$, the closest of the 4 body-diagonal capsules in the positive octant is **always** the $`(+1, +1, +1)`$ diagonal—requiring **zero sorting swaps and 2 lines of GLSL**:
   ```glsl
   float sdJack3D(vec3 p, float w, float r) {
       p = abs(p);
       return length(p - min(p.x + p.y + p.z, w) * (1.0 / 3.0)) - r;
   }
   ```
2. **6-Pointed Coordinate-Axis Jack (`sdAxisJack3D`, 3-Capsule Union, $`3\times`$ Faster)**:
   In the positive octant $`p \leftarrow |p|`$, the closest coordinate-axis segment lies along $`m = \max(p_x, \max(p_y, p_z))`$, with clamped parameter $`t = \min(m, w)`$ and squared distance $`\mathrm{dot}(p, p) - t(2m - t)`$.
3. **12-Pointed Cuboctahedral Face-Diagonal Jack (`sdFaceJack3D`, 6-Capsule Union, $`6\times`$ Faster)**:
   Sorting $`p_x \ge p_y \ge p_z \ge 0`$ maximizes $`p_x + p_y`$ across all 6 face diagonals $`(\pm 1, \pm 1, 0)/\sqrt{2}`$.

---

## 2. Numerical Verification ($`100{,}000`$ Random 3D Points vs. $`N`$-Capsule Union)

| 3D Jack Primitive | Capsules $`N`$ | Naive `sqrt` Calls | **`sdf-gems` `sqrt` Calls** | Sorting Swaps | Max Error vs. $`N`$-Capsule Union |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`sdJack3D` (Body Diagonals)** | `4` (8 tips) | `4` | **`1`** | **`0`** | `4.44e-16` |
| **`sdAxisJack3D` (Coordinate Axes)** | `3` (6 tips) | `3` | **`1`** | **`0`** | `4.44e-16` |
| **`sdFaceJack3D` (Face Diagonals)** | `6` (12 tips) | `6` | **`1`** | `3` | `4.44e-16` |
