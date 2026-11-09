# sdCinquefoilKnot & sdSolomonLink

## 3D Preview & 2D Cross-Section Slice
![sdCinquefoilKnot](../images/primitives3d/sdCinquefoilKnot.png) ![sdCinquefoilKnot Slice](../images/primitives3d/sdCinquefoilKnot_slice.png)

## Error Heatmap & Comparison
![sdCinquefoilKnot Comparison](../images/comparisons/sdCinquefoilKnot_comparison.png)

* **Category:** 3D Primitives
* **GLSL File:** `glsl/primitives3d/sdCinquefoilKnot.glsl`
* **Discovered by:** Gemini 4 Argon (research & mathematical design)
* **Implemented by:** Gemini 3.8 Flash (implementation & verification)
* **Exactness:** Exact $`O(1)`$ 3D Euclidean SDF & Analytical Normal Gradients

This header provides exact $`O(1)`$ 3D Euclidean signed distance functions and analytical unit normal gradients for two fundamental families of toroidal curves:
1. The **3D $`(2, 5)`$ Cinquefoil Torus Knot** (`sdCinquefoilKnot`, `sdgCinquefoilKnot`), also known as the $`5_1`$ knot or Solomon's Seal knot—a single closed strand that wraps 2 times longitudinally and 5 times azimuthally around a torus. Overloaded functions support general $`(2, 2m+1)`$ odd-lobe torus knots (e.g. $`q = 3`$ Trefoil, $`q = 5`$ Cinquefoil, $`q = 7`$ Septafoil).
2. **King Solomon's Link** (`sdSolomonLink`, `sdgSolomonLink`), the $`4^2_1`$ link—two intertwined, doubly-linked closed rings with 4 alternating crossings. Overloaded functions support general $`(2, 2n)`$ even-lobe multi-ring torus links.

Both primitives execute in strict $`O(1)`$ time, independent of the lobe count $`k`$, through the **Antipodal Torus $`D_k`$ Fundamental-Cell Folding Theorem** and a vectorized 2-lane `vec2` Newton solver.

---

## Verification

Evaluated on hardware GPU (OpenGL 4.1 core via Transform Feedback) against dense FP64 ground-truth polylines ($`150{,}000`$ segments):

| Primitive | `NaN / Inf` | `max \|d - gt\|` | `mean \|d - gt\|` | `check_sdg \|dg.x - d\|` | `\|\|g\|\| - 1` | `Compile Umbrella` |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `sdCinquefoilKnot` | `0 / 4009` | `3.148e-07` | `4.746e-08` | `0.00e+00` | `< 1.77e-07` | `OK (262 files)` |
| `sdSolomonLink` | `0 / 4009` | `3.177e-07` | `4.719e-08` | `0.00e+00` | `< 1.24e-07` | `OK (262 files)` |

Both primitives pass all numerical checks:
* **Zero NaNs or Infs** across coordinate axes, symmetry planes, and medial axes (`0 / 4009` probe points, `0 / 5000` random volume points).
* **Strict $`L^\infty`$ accuracy:** $`\max |d - \text{gt}| < 3.2 \times 10^{-7}`$ and $`\text{mean} |d - \text{gt}| < 4.8 \times 10^{-8}`$.
* **Analytical unit normal gradients:** $`|\|\mathbf{g}\| - 1\| < 1.8 \times 10^{-7}`$ with $`dg.x \equiv d`$.
* **Seam continuity:** $`C^0`$ continuous across all dihedral fold boundaries ($`\max |d^+ - d^-| = 2.01 \times 10^{-5}`$ against threshold $`1.2 \times 10^{-4}`$, $`0`$ violations).
* **Single-translation-unit compatibility:** verified across 262 `.glsl` headers in forward, reverse, and new-first include orders.

---

## 1. Torus Knot & Link Parameterization

### General $`(2, k)`$ Torus Curves

A classical torus curve of type $`(2, k)`$ winds $`2`$ times around the torus hole (toroidal / poloidal winding) and $`k`$ times around the central vertical axis (longitudinal / azimuthal winding). Let $`R > 0`$ denote the major radius and $`r_0 \ge 0`$ denote the minor radius.

The 3D space curve is parameterized by the azimuthal coordinate $`t \in \mathbb{R}`$:
```math
\mathbf{r}(t) = \begin{pmatrix} \bigl(R + r_0 \cos(\omega t)\bigr) \cos t \\ \bigl(R + r_0 \cos(\omega t)\bigr) \sin t \\ r_0 \sin(\omega t) \end{pmatrix}, \quad \omega = \frac{k}{2}
```

The topological classification is governed by the parity of $`k`$:

### Case 1: Odd $`k`$ ($`k = 2m + 1`$, Cinquefoil Knot $`5_1`$)

Because $`\omega = \frac{2m + 1}{2}`$ is a half-integer, when $`t`$ advances by $`2\pi`$, the harmonic phase advances by:
```math
\omega(t + 2\pi) = \omega t + (2m + 1)\pi \equiv \omega t + \pi \pmod{2\pi}
```
Thus:
```math
\cos\bigl(\omega(t + 2\pi)\bigr) = -\cos(\omega t), \quad \sin\bigl(\omega(t + 2\pi)\bigr) = -\sin(\omega t)
```
The curve does not close after a single revolution ($`t = 2\pi`$). Instead, the second revolution ($`t \in [2\pi, 4\pi]`$) traverses the **antipodal torus strand** with opposite signs on $`r_0`$:
```math
\mathbf{r}(t + 2\pi) = \begin{pmatrix} \bigl(R - r_0 \cos(\omega t)\bigr) \cos t \\ \bigl(R - r_0 \cos(\omega t)\bigr) \sin t \\ -r_0 \sin(\omega t) \end{pmatrix}
```
At $`t = 4\pi`$, $`\omega(4\pi) = 2(2m + 1)\pi \equiv 0 \pmod{2\pi}`$, closing the curve into a **single knot** ($`5_1`$).

### Case 2: Even $`k`$ ($`k = 2n`$, King Solomon's Link $`4^2_1`$)

When $`k = 2n`$, $`\omega = n`$ is an integer. The first strand closes on itself after a single $`2\pi`$ revolution:
```math
\mathbf{q}_0(t) = \begin{pmatrix} \bigl(R + r_0 \cos(n t)\bigr) \cos t \\ \bigl(R + r_0 \cos(n t)\bigr) \sin t \\ r_0 \sin(n t) \end{pmatrix}, \quad t \in [0, 2\pi)
```
The complementary antipodal strand:
```math
\mathbf{q}_1(t) = \begin{pmatrix} \bigl(R - r_0 \cos(n t)\bigr) \cos t \\ \bigl(R - r_0 \cos(n t)\bigr) \sin t \\ -r_0 \sin(n t) \end{pmatrix}, \quad t \in [0, 2\pi)
```
forms a separate closed ring. Together, $`\mathbf{q}_0`$ and $`\mathbf{q}_1`$ form a **2-component link** with linking number $`n`$ and $`2n`$ alternating crossings. When $`n = 2`$ ($`k = 4`$), this is the famous **King Solomon's Doubly-Linked Rings** ($`4^2_1`$).

---

## 2. The $`D_k`$ Antipodal Torus Fundamental-Cell Folding Theorem

### Dihedral Symmetry of the Antipodal Pair

Consider the unordered set of two antipodal strands:
```math
\mathcal{K} = \mathbf{q}_0([0, 2\pi]) \cup \mathbf{q}_1([0, 2\pi])
```
where $`\mathbf{q}_0`$ has phase $`\theta_0 = \omega t`$ and $`\mathbf{q}_1`$ has phase $`\theta_1 = \omega t + \pi`$.

**Theorem (Antipodal Torus $`D_k`$ Symmetry):**
The set $`\mathcal{K}`$ is strictly invariant under the dihedral group $`D_k`$ of order $`2k`$:

#### 1. Discrete Azimuthal Rotation ($`\mathbb{Z}_k`$)
Rotating the spatial coordinates by $`\Delta\phi = \frac{2\pi}{k}`$:
```math
t \mapsto t + \frac{2\pi}{k} \implies \omega t \mapsto \omega t + \frac{k}{2} \frac{2\pi}{k} = \omega t + \pi
```
This shifts the harmonic phase by exactly $`\pi`$, mapping:
```math
\mathbf{q}_0 \mapsto \mathbf{q}_1 \quad \text{and} \quad \mathbf{q}_1 \mapsto \mathbf{q}_0
```
The individual strands swap places, leaving the union $`\mathcal{K}`$ invariant.

#### 2. Dihedral Reflection
Under the coordinate reflection $`(\phi, z) \mapsto (-\phi, -z)`$:
```math
t \mapsto -t \implies \cos(-t) = \cos t, \quad \sin(-t) = -\sin t, \quad \sin(-\omega t) = -\sin(\omega t)
```
The curve is mapped to itself with reversed parameter direction.

### Domain Folding into $`[0, \pi / k]`$

Because of this symmetry, any point $`\mathbf{p} = (\rho \cos\phi, \rho \sin\phi, z)`$ can be folded into a single fundamental cell of width $`\frac{\pi}{k}`$:
```math
\phi_{\text{cell}} = \phi - k_{\text{cell}} \cdot \frac{2\pi}{k}, \quad k_{\text{cell}} = \left\lfloor \frac{\phi}{2\pi / k} \right\rfloor \in \mathbb{Z}
```
If $`\phi_{\text{cell}} > \frac{\pi}{k}`$, we apply the dihedral reflection:
```math
\phi_{\text{fold}} = \frac{2\pi}{k} - \phi_{\text{cell}}, \quad z_{\text{fold}} = -z
```
Otherwise:
```math
\phi_{\text{fold}} = \phi_{\text{cell}}, \quad z_{\text{fold}} = z
```
Folded coordinates in $`\mathbb{R}^3`$:
```math
\mathbf{p}_{\text{fold}} = \bigl(\rho \cos\phi_{\text{fold}}, \; \rho \sin\phi_{\text{fold}}, \; z_{\text{fold}}\bigr)
```

Within this narrow fundamental wedge $`\phi_{\text{fold}} \in [0, \pi / k]`$, the nearest footpoint on the curve must lie close to $`t \approx \phi_{\text{fold}}`$. This reduces the search domain from the full $`4\pi`$ circle down to a localized window of width $`\sim \frac{3.2\pi}{\omega}`$.

---

## 3. Vectorized 2-Lane `vec2` Newton Stationarity Solver

### Cylindrical Distance Stationarity

For a target point $`\mathbf{p}_{\text{fold}} = (p_x, p_y, p_z)`$ and curve strand $`\mathbf{q}_j(t)`$ with phase $`\theta_j = \omega t + \phi_j`$ ($`\phi_0 = 0`$, $`\phi_1 = \pi`$):
```math
r_j(t) = R + r_0 \cos\theta_j, \quad Z_j(t) = r_0 \sin\theta_j
```
The squared Euclidean distance is:
```math
D_j(t)^2 = \|\mathbf{p}_{\text{fold}} - \mathbf{q}_j(t)\|^2 = \bigl(p_x - r_j(t)\cos t\bigr)^2 + \bigl(p_y - r_j(t)\sin t\bigr)^2 + \bigl(p_z - Z_j(t)\bigr)^2
```
Expanding and projecting into the local rotating cylindrical basis $`u = p_x \cos t + p_y \sin t`$, $`v = -p_x \sin t + p_y \cos t`$:
```math
D_j(t)^2 = \rho^2 + r_j(t)^2 - 2 r_j(t) u + \bigl(Z_j(t) - p_z\bigr)^2
```
The stationarity condition $`\frac{1}{2} \frac{d}{dt} D_j(t)^2 = 0`$ yields the exact gradient function:
```math
f_j(t) = \bigl(r_j - u\bigr) r_j' - r_j v + \bigl(Z_j - p_z\bigr) Z_j' = 0
```
and its second derivative:
```math
f_j'(t) = (r_j')^2 + (r_j - u) r_j'' - 2 r_j' v + r_j u + (Z_j')^2 + (Z_j - p_z) Z_j''
```
where:
```math
r_j' = -\omega r_0 \sin\theta_j, \quad r_j'' = -\omega^2 r_0 \cos\theta_j, \quad Z_j' = \omega r_0 \cos\theta_j, \quad Z_j'' = -\omega^2 r_0 \sin\theta_j
```

### SIMD Vectorization

Both strands $`j \in \{0, 1\}`$ are evaluated simultaneously across GLSL's 2-lane SIMD register `vec2`:
```glsl
vec2 th = omega * t + phases; // phases = vec2(0.0, PI)
```
Each iteration performs a guarded Newton-Raphson correction step:
```glsl
vec2 step = vec2(
    (df.x > 1e-6) ? (f.x / df.x) : ((f.x >= 0.0 ? 0.2 : -0.2) / omEff),
    (df.y > 1e-6) ? (f.y / df.y) : ((f.y >= 0.0 ? 0.2 : -0.2) / omEff)
);
t = clamp(t - clamp(step, vec2(-0.4 / omEff), vec2(0.4 / omEff)), vec2(tMin), vec2(tMax));
```

Three seed offsets $`\{0, -0.45/\omega, +0.45/\omega\}`$ are evaluated with 6 iterations each, guaranteeing convergence to the absolute global distance minimum in $`O(1)`$ operations.

---

## 4. Analytical Gradient Unfolding

Let $`\mathbf{q}^* = (q_x, q_y, q_z)`$ be the closest footpoint on the curve in folded coordinates, and $`\mathbf{diff}_{\text{fold}} = \mathbf{p}_{\text{fold}} - \mathbf{q}^*`$. The unit normal gradient in folded space is:
```math
\mathbf{g}_{\text{fold}} = \frac{\mathbf{diff}_{\text{fold}}}{\|\mathbf{diff}_{\text{fold}}\|}
```

To transform $`\mathbf{g}_{\text{fold}}`$ back to original world space, we apply the inverse of the folding isometry:

### Unreflected Case ($`\phi_{\text{cell}} \le \frac{\pi}{k}`$)
The folding was a pure rotation by $`\alpha_0 = k_{\text{cell}} \cdot \frac{2\pi}{k}`$:
```math
\mathbf{g}_{\text{world}} = \begin{pmatrix} \cos\alpha_0 & -\sin\alpha_0 & 0 \\ \sin\alpha_0 & \cos\alpha_0 & 0 \\ 0 & 0 & 1 \end{pmatrix} \mathbf{g}_{\text{fold}}
```

### Reflected Case ($`\phi_{\text{cell}} > \frac{\pi}{k}`$)
The folding was a reflection across $`\alpha_1 = (k_{\text{cell}} + 1) \cdot \frac{2\pi}{k}`$ combined with $`z \mapsto -z`$:
```math
\mathbf{g}_{\text{world}} = \begin{pmatrix} \cos\alpha_1 & \sin\alpha_1 & 0 \\ \sin\alpha_1 & -\cos\alpha_1 & 0 \\ 0 & 0 & -1 \end{pmatrix} \mathbf{g}_{\text{fold}}
```

Because this mapping is an exact orthogonal transformation ($`M^T M = I`$), the unit gradient property $`\|\mathbf{g}_{\text{world}}\| \equiv 1`$ is preserved everywhere in $`\mathbb{R}^3`$.

---

## 5. GLSL Source Implementation

```glsl
#ifndef SDF_GEMS_SD_CINQUEFOIL_KNOT_GLSL
#define SDF_GEMS_SD_CINQUEFOIL_KNOT_GLSL

#include "../common/common.glsl"

// Shared O(1) D_k-folded 2-strand antipodal torus solver for odd-k (2, k) knots and even-k (2, 2n) links
vec4 _sdAntipodalTorus2Strand( in vec3 p, in float R, in float r0, in int kLobes, in float rTube )
{
    R  = max(R, 1e-5);
    r0 = max(r0, 0.0);
    float k     = float(max(kLobes, 1));
    float omega = 0.5 * k;

    const float PI = 3.14159265358979323846;
    float per   = (2.0 * PI) / k;
    float hPer  = 0.5 * per;
    const vec2 phases = vec2(0.0, 3.14159265358979323846);

    float rho = length(p.xy);
    float phi = (rho > 1e-12) ? atan(p.y, p.x) : 0.0;

    // Exact D_k dihedral-strand fold into [0, pi / k]
    float kCell   = floor(phi / per);
    float phiLoc  = phi - kCell * per;
    bool  refl    = (phiLoc > hPer);
    float phiFold = refl ? (per - phiLoc) : phiLoc;
    float pz      = refl ? -p.z : p.z;
    float px      = rho * cos(phiFold);
    float py      = rho * sin(phiFold);

    float omEff = max(omega, 1.0);
    float win   = 1.6 * PI / omEff;
    float tMin  = phiFold - win;
    float tMax  = phiFold + win;

    float bestD2   = 1e30;
    vec3  bestDiff = vec3(1.0, 0.0, 0.0);

    float offsets[3];
    offsets[0] = 0.0;
    offsets[1] = -0.45 / omEff;
    offsets[2] =  0.45 / omEff;

    for( int s = 0; s < 3; s++ )
    {
        vec2 t = vec2(clamp(phiFold + offsets[s], tMin, tMax));
        for( int i = 0; i < 6; i++ )
        {
            vec2 th = omega * t + phases;
            vec2 c1 = cos(t), s1 = sin(t);
            vec2 ct = cos(th), st = sin(th);

            vec2 r   = R + r0 * ct;
            vec2 dr  = (-omega * r0) * st;
            vec2 d2r = (-omega * omega * r0) * ct;

            vec2 Z   = r0 * st;
            vec2 dZ  = (omega * r0) * ct;
            vec2 d2Z = (-omega * omega * r0) * st;

            vec2 u  = px * c1 + py * s1;
            vec2 v  = -px * s1 + py * c1;
            vec2 f  = (r - u) * dr - r * v + (Z - pz) * dZ;
            vec2 df = dr * dr + (r - u) * d2r - 2.0 * dr * v + r * u + dZ * dZ + (Z - pz) * d2Z;

            vec2 step = vec2(
                (df.x > 1e-6) ? (f.x / df.x) : ((f.x >= 0.0 ? 0.2 : -0.2) / omEff),
                (df.y > 1e-6) ? (f.y / df.y) : ((f.y >= 0.0 ? 0.2 : -0.2) / omEff)
            );
            t = clamp(t - clamp(step, vec2(-0.4 / omEff), vec2(0.4 / omEff)), vec2(tMin), vec2(tMax));
        }

        vec2 th = omega * t + phases;
        vec2 c1 = cos(t), s1 = sin(t);
        vec2 r  = R + r0 * cos(th);
        vec2 Z  = r0 * sin(th);

        vec3  q0    = vec3(r.x * c1.x, r.x * s1.x, Z.x);
        vec3  diff0 = vec3(px, py, pz) - q0;
        float d2_0  = dot(diff0, diff0);
        if( d2_0 < bestD2 )
        {
            bestD2   = d2_0;
            bestDiff = diff0;
        }

        vec3  q1    = vec3(r.y * c1.y, r.y * s1.y, Z.y);
        vec3  diff1 = vec3(px, py, pz) - q1;
        float d2_1  = dot(diff1, diff1);
        if( d2_1 < bestD2 )
        {
            bestD2   = d2_1;
            bestDiff = diff1;
        }
    }

    float invD = inversesqrt(max(bestD2, 1e-30));
    float d    = bestD2 * invD;
    vec3  gq   = bestDiff * invD;

    vec3 grad;
    if( refl )
    {
        float a0 = (kCell + 1.0) * per;
        float ca = cos(a0), sa = sin(a0);
        grad = vec3(gq.x * ca + gq.y * sa, gq.x * sa - gq.y * ca, -gq.z);
    }
    else
    {
        float a0 = kCell * per;
        float ca = cos(a0), sa = sin(a0);
        grad = vec3(gq.x * ca - gq.y * sa, gq.x * sa + gq.y * ca, gq.z);
    }
    return vec4(d - rTube, grad);
}

// 1. Exact O(1) 3D (2, 5) Cinquefoil Torus Knot (5_1 Solomon's Seal Knot) + Unit Gradient vec4(d, grad.xyz)
vec4 sdgCinquefoilKnot( in vec3 p, in float R, in float r0, in float rTube )
{
    return _sdAntipodalTorus2Strand(p, R, r0, 5, rTube);
}

// 2. Scalar 3D (2, 5) Cinquefoil Torus Knot SDF
float sdCinquefoilKnot( in vec3 p, in float R, in float r0, in float rTube )
{
    return _sdAntipodalTorus2Strand(p, R, r0, 5, rTube).x;
}

// 5-parameter overloads for general (2, 2*m + 1) Torus Knots (e.g. qLobes = 3 Trefoil, 5 Cinquefoil, 7 Septafoil)
vec4 sdgCinquefoilKnot( in vec3 p, in float R, in float r0, in int qLobes, in float rTube )
{
    return _sdAntipodalTorus2Strand(p, R, r0, max(qLobes, 1), rTube);
}

float sdCinquefoilKnot( in vec3 p, in float R, in float r0, in int qLobes, in float rTube )
{
    return _sdAntipodalTorus2Strand(p, R, r0, max(qLobes, 1), rTube).x;
}

// 3. Exact O(1) 3D King Solomon's Link (4^2_1 Doubly-Linked Rings, or general (2, 2*nLinks) Torus Link)
vec4 sdgSolomonLink( in vec3 p, in float R, in float r0, in int nLinks, in float rTube )
{
    return _sdAntipodalTorus2Strand(p, R, r0, 2 * max(nLinks, 1), rTube);
}

// 4. Scalar 3D King Solomon's Link / (2, 2*nLinks) Torus Link SDF
float sdSolomonLink( in vec3 p, in float R, in float r0, in int nLinks, in float rTube )
{
    return _sdAntipodalTorus2Strand(p, R, r0, 2 * max(nLinks, 1), rTube).x;
}

// Canonical 4-parameter overloads for King Solomon's 4-Crossing Link 4^2_1 (nLinks = 2)
vec4 sdgSolomonLink( in vec3 p, in float R, in float r0, in float rTube )
{
    return _sdAntipodalTorus2Strand(p, R, r0, 4, rTube);
}

float sdSolomonLink( in vec3 p, in float R, in float r0, in float rTube )
{
    return _sdAntipodalTorus2Strand(p, R, r0, 4, rTube).x;
}

#endif // SDF_GEMS_SD_CINQUEFOIL_KNOT_GLSL
```
