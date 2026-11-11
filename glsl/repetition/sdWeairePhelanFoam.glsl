// ============================================================================
// sdWeairePhelanFoam, sdgWeairePhelanFoam, sdWeairePhelanWalls, opRepA15
// 0-Loop Exact 3D Weaire-Phelan A15 (Pm-3n) Minimal-Area Foam & Clathrate SDFs
//
// Implements Denis Weaire & Robert Phelan's (1994) space-filling minimal-area
// cellular partition (0.3% lower surface area than Lord Kelvin's 1887 BCC foam;
// realized architecturally in the Beijing National Aquatics Center / Water Cube),
// based on the crystallographic A15 (beta-W, space group Pm-3n, No. 223)
// Frank-Kasper clathrate with 8 cells per cubic period:
//   - 2 pyritohedral dodecahedra (12 pentagons) at Wyckoff 2a sites:
//       (0, 0, 0) and (1/2, 1/2, 1/2)
//   - 6 Goldberg tetrakaidecahedra (2 hexagons + 12 pentagons) at Wyckoff 6d:
//       (1/4, 1/2, 0), (3/4, 1/2, 0) and cyclic (x, y, z) permutations
//
// @author Discovered and derived using Gemini 4 Carbon (sdf-gems, 2026)
// @exact  Yes (0-loop, 0-sqrt exact A15 Voronoi/Laguerre repetition & wall
//         skeleton distance; 0-loop, 1-sqrt exact 3D Euclidean distance and
//         unit analytical gradient to all 92 Plateau-border struts per cell)
// ============================================================================

#ifndef SDF_GEMS_HAS_WEAIRE_PHELAN_FOAM
#define SDF_GEMS_HAS_WEAIRE_PHELAN_FOAM

#include "../common/common.glsl"

// Algebraic Laguerre power-weight shift (4 * 2^(1/3) - 5) / 16 equalizing the
// Wyckoff 2a dodecahedron and 6d tetrakaidecahedron volumes to V_2a = V_6d = s^3 / 8
// (unweighted A15 Voronoi delta = 0 gives V_2a = 125/1024 * s^3, V_6d = 129/1024 * s^3).
#ifndef SDF_WEAIRE_PHELAN_EQUAL_VOL_DELTA
#define SDF_WEAIRE_PHELAN_EQUAL_VOL_DELTA 0.0024802625
#endif

// ----------------------------------------------------------------------------
// Internal Helper: 48-Fold Pm-3n Fundamental Chamber Folding & Gradient Unfold
// Folds u in [-0.5, 0.5]^3 into the 1/48 convex chamber C_48:
//   0 <= w.y, w.z <= w.x <= 0.5,   w.x + w.y + w.z <= 0.75
// using 8 octant mirror reflections, the non-symmorphic (1/2, 1/2, 1/2) n-glide
// swap w -> 0.5 - w.yxz, and 3-fold C_3 cyclic axis permutation.
// ----------------------------------------------------------------------------
#ifndef SDF_GEMS_HAS_PM3N_FOLD
#define SDF_GEMS_HAS_PM3N_FOLD
vec3 sdfGemsFoldPm3n( in vec3 u, out vec3 sg, out bool glide, out int cyc )
{
    sg = vec3(greaterThanEqual(u, vec3(0.0))) * 2.0 - 1.0;
    vec3 w = abs(u);

    glide = (w.x + w.y + w.z > 0.75);
    if (glide) w = 0.5 - w.yxz;

    if (w.y > w.x && w.y >= w.z) {
        cyc = 1;
        return w.yzx;
    }
    if (w.z > w.x && w.z > w.y) {
        cyc = 2;
        return w.zxy;
    }
    cyc = 0;
    return w;
}

vec3 sdfGemsUnfoldPm3nGrad( in vec3 gw, in vec3 sg, in bool glide, in int cyc )
{
    vec3 g = (cyc == 1) ? gw.zxy : ((cyc == 2) ? gw.yzx : gw);
    if (glide) g = -g.yxz;
    return sg * g;
}
#endif // SDF_GEMS_HAS_PM3N_FOLD

// ----------------------------------------------------------------------------
// 1. opRepA15 — O(1) 0-Loop, 0-Sqrt A15 (Pm-3n) Voronoi & Laguerre Repetition
//
// Exploits 48-fold Pm-3n space-group symmetry to reduce all 216 neighbor sites
// in a 3x3x3 periodic supercell to 3 linear half-space comparisons in C_48:
//   - cellType = 0 : Wyckoff 2a pyritohedral dodecahedron
//   - cellType = 1 : Wyckoff 6d Goldberg tetrakaidecahedron
// Writes the exact fractional world site center to cellID and returns the local
// world-oriented offset q = p - cellID * scale.
// ----------------------------------------------------------------------------
vec3 opRepA15( in vec3 p, in float scale, in float laguerreDelta, out int cellType, out vec3 cellID )
{
    float s  = max(scale, 1e-12);
    vec3  u0 = p / s;
    vec3  n  = round(u0);
    vec3  u  = u0 - n;

    vec3 sg;
    bool glide;
    int  cyc;
    vec3 w = sdfGemsFoldPm3n(u, sg, glide, cyc);

    float d  = clamp(laguerreDelta, -0.05, 0.05);
    float la = -d;
    float lb = 0.3125 - w.x - 0.5 * w.z;
    float lc = 0.3125 - 0.5 * w.x - w.y;

    vec3 cw;
    if (la <= lb && la <= lc) {
        cellType = 0;
        cw = vec3(0.0);
    } else if (lb <= lc) {
        cellType = 1;
        cw = vec3(0.5, 0.0, 0.25);
    } else {
        cellType = 1;
        cw = vec3(0.25, 0.5, 0.0);
    }

    vec3 dw = w - cw;
    vec3 cg = (cyc == 1) ? cw.zxy : ((cyc == 2) ? cw.yzx : cw);
    vec3 dg = (cyc == 1) ? dw.zxy : ((cyc == 2) ? dw.yzx : dw);

    vec3 co = glide ? (0.5 - cg.yxz) : cg;
    vec3 do_ = glide ? (-dg.yxz) : dg;

    cellID = n + sg * co;
    return (sg * do_) * s;
}

vec3 opRepA15( in vec3 p, in float scale, out int cellType )
{
    vec3 cellID;
    return opRepA15(p, scale, 0.0, cellType, cellID);
}

vec3 opRepA15( in vec3 p, in float scale )
{
    int  cellType;
    vec3 cellID;
    return opRepA15(p, scale, 0.0, cellType, cellID);
}

// ----------------------------------------------------------------------------
// 2. sdWeairePhelanWalls & sdgWeairePhelanWalls
//
// 0-Loop, 0-Sqrt Polyhedral Cell-Wall Sheet Network of the A15 Voronoi /
// Equal-Volume Laguerre Diagram (thickened by half-thickness wallThickness).
// Pass laguerreDelta = 0.0 for unweighted A15 Voronoi walls, or
// SDF_WEAIRE_PHELAN_EQUAL_VOL_DELTA (0.0024802625) for equal-volume V_2a = V_6d
// Weaire-Phelan polyhedral walls (A / V^(2/3) = 5.29695 vs Kelvin 5.31474).
// Exact Euclidean distance to the 54-face wall skeleton at wallThickness = 0,
// and exact exterior / conservative 1-Lipschitz interior bound for > 0.
// ----------------------------------------------------------------------------
vec4 sdgWeairePhelanWalls( in vec3 p, in float scale, in float wallThickness, in float laguerreDelta )
{
    float s  = max(scale, 1e-12);
    vec3  u0 = p / s;
    vec3  u  = u0 - round(u0);

    vec3 sg;
    bool glide;
    int  cyc;
    vec3 w = sdfGemsFoldPm3n(u, sg, glide, cyc);

    float d  = clamp(laguerreDelta, -0.05, 0.05);
    float z  = 0.3125 + d;
    float la = -d;
    float lb = 0.3125 - w.x - 0.5 * w.z;
    float lc = 0.3125 - 0.5 * w.x - w.y;

    const float k5  = 0.89442719; // 2 / sqrt(5)
    const float h5  = 0.44721360; // 1 / sqrt(5)
    const float k6  = 0.40824829; // 1 / sqrt(6)
    const float k6b = 0.81649658; // 2 / sqrt(6)

    float dist_w;
    vec3  gw;

    if (la <= lb && la <= lc) {
        // Cell 2a at (0, 0, 0): 2 LP-certified active closest wall planes in C_48
        float dAB = (z - (w.x + 0.5 * w.z)) * k5;
        float dAC = (z - (0.5 * w.x + w.y)) * k5;
        if (dAB <= dAC) {
            dist_w = dAB;
            gw     = vec3(-k5, 0.0, -h5);
        } else {
            dist_w = dAC;
            gw     = vec3(-h5, -k5, 0.0);
        }
    } else if (lb <= lc) {
        // Cell 6d at (0.5, 0, 0.25): 5 LP-certified active closest wall planes in C_48
        dist_w = ((w.x + 0.5 * w.z) - z) * k5;
        gw     = vec3(k5, 0.0, h5);

        float dB1 = w.z;
        if (dB1 < dist_w) { dist_w = dB1; gw = vec3(0.0, 0.0, 1.0); }

        float dB2 = 0.5 - w.z;
        if (dB2 < dist_w) { dist_w = dB2; gw = vec3(0.0, 0.0, -1.0); }

        float dBC = (w.x - 2.0 * w.y + w.z) * k6;
        if (dBC < dist_w) { dist_w = dBC; gw = vec3(k6, -k6b, k6); }

        float dBD = (2.0 * w.x - w.y - w.z) * k6;
        if (dBD < dist_w) { dist_w = dBD; gw = vec3(k6b, -k6, -k6); }
    } else {
        // Cell 6d at (0.25, 0.5, 0): 3 LP-certified active closest wall planes in C_48
        dist_w = ((0.5 * w.x + w.y) - z) * k5;
        gw     = vec3(h5, k5, 0.0);

        float dC1 = 0.5 - w.x;
        if (dC1 < dist_w) { dist_w = dC1; gw = vec3(-1.0, 0.0, 0.0); }

        float dCB = (2.0 * w.y - w.x - w.z) * k6;
        if (dCB < dist_w) { dist_w = dCB; gw = vec3(-k6, k6b, -k6); }
    }

    vec3 grad = sdfGemsUnfoldPm3nGrad(gw, sg, glide, cyc);
    return vec4(max(dist_w, 0.0) * s - wallThickness, grad);
}

float sdWeairePhelanWalls( in vec3 p, in float scale, in float wallThickness, in float laguerreDelta )
{
    return sdgWeairePhelanWalls(p, scale, wallThickness, laguerreDelta).x;
}

float sdWeairePhelanWalls( in vec3 p, in float scale, in float wallThickness )
{
    return sdgWeairePhelanWalls(p, scale, wallThickness, SDF_WEAIRE_PHELAN_EQUAL_VOL_DELTA).x;
}

// ----------------------------------------------------------------------------
// 3. sdWeairePhelanFoam & sdgWeairePhelanFoam
//
// 0-Loop, 1-Sqrt Exact 3D Euclidean SDF & Analytical Unit Normal to the 92
// Plateau-Border Struts per Unit Cell of the A15 / Weaire-Phelan Open-Cell Foam.
//
// By the Right-Angled Mirror-Box Theorem on [0, 0.5]^3 combined with 48-fold
// Pm-3n folding into C_48, all 3,497 supercell struts collapse without error
// onto 8 constant-rational-direction segments in C_48:
//   - 3 axis-aligned segments along (0,1,0) and (1,0,0)
//   - 2 coordinate-plane 6c apex struts along (2,1,0) and (1,0,2)
//   - 1 body-diagonal 16i-16i strut along (1,1,1)
//   - 2 general 16i-24k struts along (2,-1,-4) and (-1,-4,2)
//
// The 3-argument overloads evaluate the unweighted A15 Voronoi strut network
// (laguerreDelta = 0.0); pass SDF_WEAIRE_PHELAN_EQUAL_VOL_DELTA to the 4-arg
// overloads for the equal-volume Laguerre Weaire-Phelan strut network.
// ----------------------------------------------------------------------------
float sdWeairePhelanFoam( in vec3 p, in float scale, in float strutRadius, in float laguerreDelta )
{
    float s = max(scale, 1e-12);
    vec3  u = p / s;
    vec3  w = abs(u - round(u));
    if (w.x + w.y + w.z > 0.75) w = 0.5 - w.yxz;
    if (w.y > w.x && w.y >= w.z)      w = w.yzx;
    else if (w.z > w.x && w.z > w.y)  w = w.zxy;

    float z = 0.3125 + clamp(laguerreDelta, -0.05, 0.05);
    float y = 0.5 * z;
    float a = (2.0 / 3.0) * z;

    // S1: (z, 0, 0) -> (z, y, 0)
    vec3  d1 = vec3(w.x - z, max(w.y - y, 0.0), w.z);
    float m2 = dot(d1, d1);

    // S2: (0, 0, z) -> (y, 0, z)
    vec3  d2 = vec3(max(w.x - y, 0.0), w.y, w.z - z);
    m2 = min(m2, dot(d2, d2));

    // S3: (0.5 - y, 0.5 - z, 0.5) -> (0.5, 0.5 - z, 0.5)
    vec3  d3 = vec3(min(w.x - (0.5 - y), 0.0), w.y - (0.5 - z), w.z - 0.5);
    m2 = min(m2, dot(d3, d3));

    // S4: (z, y, 0) -> (0.5, 0.25, 0), direction (-2, -1, 0) from (0.5, 0.25, 0)
    float t4 = clamp(0.25 - 0.4 * w.x - 0.2 * w.y, 0.0, 0.25 - y);
    vec3  d4 = vec3(w.x - (0.5 - 2.0 * t4), w.y - (0.25 - t4), w.z);
    m2 = min(m2, dot(d4, d4));

    // S5: (y, 0, z) -> (0.25, 0, 0.5), direction (-1, 0, -2) from (0.25, 0, 0.5)
    float t5 = clamp(0.25 - 0.2 * w.x - 0.4 * w.z, 0.0, 0.25 - y);
    vec3  d5 = vec3(w.x - (0.25 - t5), w.y, w.z - (0.5 - 2.0 * t5));
    m2 = min(m2, dot(d5, d5));

    // S6: (a, a, a) -> (0.5 - a, 0.5 - a, 0.5 - a), direction (1, 1, 1)
    float t6 = clamp((w.x + w.y + w.z) * (1.0 / 3.0), a, 0.5 - a);
    vec3  d6 = w - vec3(t6);
    m2 = min(m2, dot(d6, d6));

    // S7: (a, a, a) -> (z, y, 0), direction (2, -1, -4)
    float t7 = clamp((2.0 * w.x - w.y - 4.0 * w.z + 3.0 * a) * (1.0 / 21.0), 0.0, 0.25 * a);
    vec3  d7 = w - vec3(a + 2.0 * t7, a - t7, a - 4.0 * t7);
    m2 = min(m2, dot(d7, d7));

    // S8: (a, a, a) -> (y, 0, z), direction (-1, -4, 2)
    float t8 = clamp((-w.x - 4.0 * w.y + 2.0 * w.z + 3.0 * a) * (1.0 / 21.0), 0.0, 0.25 * a);
    vec3  d8 = w - vec3(a - t8, a - 4.0 * t8, a + 2.0 * t8);
    m2 = min(m2, dot(d8, d8));

    return sqrt(m2) * s - strutRadius;
}

vec4 sdgWeairePhelanFoam( in vec3 p, in float scale, in float strutRadius, in float laguerreDelta )
{
    float s  = max(scale, 1e-12);
    vec3  u0 = p / s;
    vec3  u  = u0 - round(u0);

    vec3 sg;
    bool glide;
    int  cyc;
    vec3 w = sdfGemsFoldPm3n(u, sg, glide, cyc);

    float z = 0.3125 + clamp(laguerreDelta, -0.05, 0.05);
    float y = 0.5 * z;
    float a = (2.0 / 3.0) * z;

    vec3  d_win = vec3(w.x - z, max(w.y - y, 0.0), w.z);
    float m2    = dot(d_win, d_win);
    int   k_win = 1;

    vec3  d2 = vec3(max(w.x - y, 0.0), w.y, w.z - z);
    float s2 = dot(d2, d2);
    if (s2 < m2) { m2 = s2; d_win = d2; k_win = 2; }

    vec3  d3 = vec3(min(w.x - (0.5 - y), 0.0), w.y - (0.5 - z), w.z - 0.5);
    float s3 = dot(d3, d3);
    if (s3 < m2) { m2 = s3; d_win = d3; k_win = 3; }

    float t4 = clamp(0.25 - 0.4 * w.x - 0.2 * w.y, 0.0, 0.25 - y);
    vec3  d4 = vec3(w.x - (0.5 - 2.0 * t4), w.y - (0.25 - t4), w.z);
    float s4 = dot(d4, d4);
    if (s4 < m2) { m2 = s4; d_win = d4; k_win = 4; }

    float t5 = clamp(0.25 - 0.2 * w.x - 0.4 * w.z, 0.0, 0.25 - y);
    vec3  d5 = vec3(w.x - (0.25 - t5), w.y, w.z - (0.5 - 2.0 * t5));
    float s5 = dot(d5, d5);
    if (s5 < m2) { m2 = s5; d_win = d5; k_win = 5; }

    float t6 = clamp((w.x + w.y + w.z) * (1.0 / 3.0), a, 0.5 - a);
    vec3  d6 = w - vec3(t6);
    float s6 = dot(d6, d6);
    if (s6 < m2) { m2 = s6; d_win = d6; k_win = 6; }

    float t7 = clamp((2.0 * w.x - w.y - 4.0 * w.z + 3.0 * a) * (1.0 / 21.0), 0.0, 0.25 * a);
    vec3  d7 = w - vec3(a + 2.0 * t7, a - t7, a - 4.0 * t7);
    float s7 = dot(d7, d7);
    if (s7 < m2) { m2 = s7; d_win = d7; k_win = 7; }

    float t8 = clamp((-w.x - 4.0 * w.y + 2.0 * w.z + 3.0 * a) * (1.0 / 21.0), 0.0, 0.25 * a);
    vec3  d8 = w - vec3(a - t8, a - 4.0 * t8, a + 2.0 * t8);
    float s8 = dot(d8, d8);
    if (s8 < m2) { m2 = s8; d_win = d8; k_win = 8; }

    float dist;
    vec3  gw;
    if (m2 > 1e-14) {
        float inv = inversesqrt(m2);
        dist = m2 * inv;
        gw   = d_win * inv;
    } else {
        dist = 0.0;
        // Exact unit normals orthogonal to the winning strut medial axis
        gw = (k_win == 1)                ? vec3(1.0, 0.0, 0.0) :
             (k_win == 2 || k_win == 3)  ? vec3(0.0, 1.0, 0.0) :
             (k_win == 4)                ? vec3(0.0, 0.0, 1.0) :
             (k_win == 5)                ? vec3(0.0, 1.0, 0.0) :
             (k_win == 6)                ? vec3(0.70710678, -0.70710678, 0.0) :
             (k_win == 7)                ? vec3(0.44721360, 0.89442719, 0.0)
                                         : vec3(0.0, 0.44721360, 0.89442719);
    }

    vec3 grad = sdfGemsUnfoldPm3nGrad(gw, sg, glide, cyc);
    return vec4(dist * s - strutRadius, grad);
}

float sdWeairePhelanFoam( in vec3 p, in float scale, in float strutRadius )
{
    return sdWeairePhelanFoam(p, scale, strutRadius, 0.0);
}

vec4 sdgWeairePhelanFoam( in vec3 p, in float scale, in float strutRadius )
{
    return sdgWeairePhelanFoam(p, scale, strutRadius, 0.0);
}

#endif // SDF_GEMS_HAS_WEAIRE_PHELAN_FOAM
