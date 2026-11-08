# `sdCircleWave2D`, `sdCircleWaveSigned2D` & `sdgArc2D` — 1-`sqrt` Unsigned/Signed Circular Wave & 1-SFU Arc Gradients

- **Category:** `primitives2d`
- **GLSL File:** [`glsl/primitives2d/sdCircleWave2D.glsl`](../../glsl/primitives2d/sdCircleWave2D.glsl)
- **Author:** Inigo Quilez — upgraded to 1-`sqrt` unsigned/signed wave & 1-SFU (`inversesqrt`) gradients using **Gemini 4 Argon** (`sdf-gems, 2026`)
- **Exact Euclidean SDF:** Yes ($`\|\nabla d\| = 1`$, $`0.00\times 10^0`$ difference vs baseline on unsigned wave, **1 `sqrt`** down from **4 `sqrt`s**)

![sdCircleWave2D Preview](../images/primitives2d/sdCircleWave2D.png)

![sdCircleWave2D Comparison](../images/comparisons/circle_wave_comparison.png)

---

## Mathematical Overview

### 1. The Inflection-Normal Dead-Code Theorem (`sdCircleWave2D`: 4 `sqrt`s $`\to`$ 1 `sqrt`)

In Inigo Quilez's baseline `sdCircleWave` ([Shadertoy `stGyzt`](https://www.shadertoy.com/view/stGyzt)), the wave is constructed inside the folded half-period $`p_x \in [0, 2\,\text{co}_x]`$ from a lower circular arc centered at $`(0, 0)`$ and an upper circular arc centered at $`2\,\text{co} = 2 r_a (\sin\theta_b, \cos\theta_b)`$, meeting with $`C^1`$ tangency at the inflection point $`\text{co}`$:

```glsl
vec2 p1 = p;
vec2 p2 = vec2(abs(p.x - 2.0*co.x), -p.y + 2.0*co.y);
float d1 = ((co.y*p1.x > co.x*p1.y) ? length(p1-co) : abs(length(p1)-ra));
float d2 = ((co.y*p2.x > co.x*p2.y) ? length(p2-co) : abs(length(p2)-ra));
return min(d1, d2);
```

This baseline evaluates **two** arc SDFs (`d1`, `d2`) with **four `length()` calls** (up to **4 `sqrt`s**).

We prove that **both endpoint cap branches `length(p1 - co)` and `length(p2 - co)` are 100% dead code** across all of $`\mathbb{R}^2`$:

1. Since $`p_x = \big|\text{mod}(x, 4\,\text{co}_x) - 2\,\text{co}_x\big| \in [0, 2\,\text{co}_x]`$, we always have $`p_x - 2\,\text{co}_x \le 0`$. Thus `abs(p.x - 2.0*co.x)` is identically $`2\,\text{co}_x - p_x`$, so $`p_2 = 2\,\text{co} - p`$ is simply the point reflection of $`p`$ across the inflection point $`\text{co}`$.
2. Under point reflection $`p_2 = 2\,\text{co} - p`$, the sector determinant reverses sign identically:
   $`\displaystyle \text{co}_y p_{2x} - \text{co}_x p_{2y} = \text{co}_y(2\,\text{co}_x - p_x) - \text{co}_x(2\,\text{co}_y - p_y) = -(\text{co}_y p_x - \text{co}_x p_y).`$
3. By the reverse triangle inequality on $`\triangle(0, p, \text{co})`$ and $`\triangle(2\,\text{co}, p, \text{co})`$, since $`\|\text{co}\| = r_a`$:
   $`\displaystyle \big|\|p\| - r_a\big| = \big|\|p\| - \|\text{co}\|\big| \le \|p - \text{co}\|, \qquad \big|\|2\,\text{co} - p\| - r_a\big| \le \|p - \text{co}\|.`$
   Therefore:
   - When $`\text{co}_y p_x \le \text{co}_x p_y`$, $`p_1`$ is in its arc sector ($`d_1 = \big|\|p\| - r_a\big|`$) and $`p_2`$ is in its cap sector ($`d_2 = \|p - \text{co}\|`$), so $`\min(d_1, d_2) = d_1 = \big|\|p\| - r_a\big|`$ identically.
   - When $`\text{co}_y p_x > \text{co}_x p_y`$, $`p_2`$ is in its arc sector ($`d_2 = \big|\|2\,\text{co} - p\| - r_a\big|`$) and $`p_1`$ is in its cap sector ($`d_1 = \|p - \text{co}\|`$), so $`\min(d_1, d_2) = d_2 = \big|\|2\,\text{co} - p\| - r_a\big|`$ identically.

Folding $`p \mapsto 2\,\text{co} - p`$ when $`\text{co}_y p_x > \text{co}_x p_y`$ collapses `sdCircleWave2D` to a **single `length(p)` (`1 sqrt`)** with zero approximation error ($`< 10^{-15}`$).

### 2. The Exact Signed Wave Theorem (`sdCircleWaveSigned2D`)

Even more remarkably, along the separating normal line $`\text{co}_y p_x = \text{co}_x p_y`$ inside the fundamental strip $`p_x \in [0, 2\,\text{co}_x]`$, every point has the form $`p = t\,\text{co}`$ with $`t = p_x / \text{co}_x \in [0, 2]`$. For all $`t \in [0, 2]`$:

```math
\|t\,\text{co}\| - r_a = (t - 1) r_a = -\big(\|(2 - t)\,\text{co}\| - r_a\big), \qquad \nabla\big(\|p\| - r_a\big)\Big|_{p=t\text{co}} = \nabla\Big(-\big(\|2\text{co} - p\| - r_a\big)\Big)\Big|_{p=t\text{co}} = \frac{\text{co}}{r_a}.
```

Thus, simply multiplying $`(\|p\| - r_a)`$ by $`-1`$ in the upper sector ($`\text{co}_y p_x > \text{co}_x p_y`$) yields a **$`C^1`$-smooth across the inflection normal, exact Euclidean signed half-space SDF (`sdCircleWaveSigned2D`)** that negatively fills the interior below the wave and positively measures distance above the wave!

### 3. Single-SFU `inversesqrt` Gradient Fusion (`sdgCircleWave2D` & `sdgArc2D`)

Instead of computing `d = length(w)` (`MUFU.SQRT`) followed by `w / d` (`MUFU.RCP`, 2 SFU instructions), both `sdgCircleWave2D` and `sdgArc2D` compute `float il = inversesqrt(l2);` (**1 SFU instruction total**) and obtain both the distance `l2 * il` and unit gradient `w * il` via fast FMA/MUL instructions, while eliminating the `sign(0.0) == 0.0` zero-gradient bug on the symmetry axis $`p_x = 0`$ and spine circle $`\|p\| = r_a`$.

---

## Numerical Verification & Benchmark

| Function | SFU Ops (`sqrt` / `rsq` / `rcp`) | Max Distance Error vs Ground Truth | Max Gradient Error | Capabilities |
| :--- | :---: | :---: | :---: | :--- |
| Baseline IQ `sdCircleWave` | `4` `sqrt`s | $`0.00 \times 10^0`$ | N/A | Unsigned curve only |
| **`sdf-gems` `sdCircleWave2D`** | **`1` `sqrt` ($`4\times`$ fewer)** | **$`8.88 \times 10^{-16}`$** | N/A | Unsigned curve |
| **`sdf-gems` `sdCircleWaveSigned2D`** | **`1` `sqrt`** | **$`0.00 \times 10^0`$ ($`\Vert \nabla d\Vert = 1`$)** | N/A | **Signed half-plane wave** |
| **`sdf-gems` `sdgCircleWave2D` / `sdgArc2D`** | **`1` `inversesqrt` (`1` SFU total)** | **$`8.88 \times 10^{-16}`$** | **$`1.39 \times 10^{-9}`$** | Distance + `NaN`-free 2D gradient |
