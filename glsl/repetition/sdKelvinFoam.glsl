// ============================================================================
// sdKelvinFoam, sdRhombicFoam, sdOctetTruss, sdBCCTruss
// 0-Loop, 1-Sqrt Exact 3D Voronoi Foams & Architected Crystal Truss Metamaterials
//
// @author Discovered using Gemini 4 Argon; implemented by Gemini 3.8 Flash (sdf-gems, 2026)
// @exact  Yes (exact Euclidean distance to infinite 3D strut networks minus radius r)
// ============================================================================

#ifndef SDF_GEMS_HAS_KELVIN_FOAM
#define SDF_GEMS_HAS_KELVIN_FOAM

#include "../common/common.glsl"
#include "opRepBCC.glsl"
#include "opRepFCC.glsl"

// Forward declarations for single-translation-unit compile order independence
vec3 opRepBCC( in vec3 p, in float s );
vec3 opRepFCC( in vec3 p, in float s );

// ----------------------------------------------------------------------------
// Internal Helper: 3-Element Descending Sort & Permutation Tracking
// Sorts |q| into descending order w.x >= w.y >= w.z >= 0, storing the original
// source component indices in idx such that w[i] = |q|[idx[i]].
// ----------------------------------------------------------------------------
#ifndef SDF_GEMS_HAS_SORT3_DESC
#define SDF_GEMS_HAS_SORT3_DESC
void sdfGemsSort3Desc( inout vec3 w, out ivec3 idx )
{
    idx = ivec3(0, 1, 2);
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; int ti = idx.x; idx.x = idx.y; idx.y = ti; }
    if (w.y < w.z) { float t = w.y; w.y = w.z; w.z = t; int ti = idx.y; idx.y = idx.z; idx.z = ti; }
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; int ti = idx.x; idx.x = idx.y; idx.y = ti; }
}

vec3 sdfGemsUnpermuteGrad3( in vec3 gw, in ivec3 idx, in vec3 q )
{
    vec3 ga;
    ga[idx.x] = gw.x;
    ga[idx.y] = gw.y;
    ga[idx.z] = gw.z;
    vec3 sgn = vec3(greaterThanEqual(q, vec3(0.0))) * 2.0 - 1.0;
    return sgn * ga;
}
#endif // SDF_GEMS_HAS_SORT3_DESC

// ----------------------------------------------------------------------------
// sdKelvinFoam — 0-Loop, 1-Sqrt Exact 3D Kelvin Tetrakaidecahedral Foam SDF
// Evaluates the exact 3D Euclidean distance to all 36 Plateau-border struts of
// Lord Kelvin's (1887) BCC truncated octahedron open-cell foam.
// ----------------------------------------------------------------------------
float sdKelvinFoam( in vec3 p, in float scale, in float r )
{
    vec3 q = opRepBCC(p, scale);
    vec3 w = abs(q);
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; }
    if (w.y < w.z) { float t = w.y; w.y = w.z; w.z = t; }
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; }

    // Edge 1 (square-hexagon strut): (0.5*s, 0.25*s - t1, t1), t1 in [0, 0.125*s]
    float t1 = clamp(0.5 * (0.25 * scale - w.y + w.z), 0.0, 0.125 * scale);
    vec3 d1 = w - vec3(0.5 * scale, 0.25 * scale - t1, t1);

    // Edge 2 (hexagon-hexagon strut): (0.5*s - t2, 0.25*s + t2, 0), t2 in [0, 0.125*s]
    float t2 = clamp(0.5 * (0.25 * scale - w.x + w.y), 0.0, 0.125 * scale);
    vec3 d2 = w - vec3(0.5 * scale - t2, 0.25 * scale + t2, 0.0);

    return sqrt(min(dot(d1, d1), dot(d2, d2))) - r;
}

// ----------------------------------------------------------------------------
// sdgKelvinFoam — Exact 3D Distance & Analytical Unit Normal of Kelvin Foam
// Returns vec4(d, normal.xyz).
// ----------------------------------------------------------------------------
vec4 sdgKelvinFoam( in vec3 p, in float scale, in float r )
{
    vec3 q = opRepBCC(p, scale);
    vec3 w = abs(q);
    ivec3 idx;
    sdfGemsSort3Desc(w, idx);

    float t1 = clamp(0.5 * (0.25 * scale - w.y + w.z), 0.0, 0.125 * scale);
    vec3 d1 = w - vec3(0.5 * scale, 0.25 * scale - t1, t1);
    float len1_sq = dot(d1, d1);

    float t2 = clamp(0.5 * (0.25 * scale - w.x + w.y), 0.0, 0.125 * scale);
    vec3 d2 = w - vec3(0.5 * scale - t2, 0.25 * scale + t2, 0.0);
    float len2_sq = dot(d2, d2);

    bool use1 = (len1_sq <= len2_sq);
    vec3 d_win = use1 ? d1 : d2;
    float dist = sqrt(use1 ? len1_sq : len2_sq);

    vec3 gw = (dist > 1e-7 * scale) ? (d_win / dist) : (use1 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 0.0, 1.0));
    vec3 grad = sdfGemsUnpermuteGrad3(gw, idx, q);

    return vec4(dist - r, grad);
}

// ----------------------------------------------------------------------------
// sdRhombicFoam — 0-Loop, 1-Sqrt Exact 3D FCC Rhombic Dodecahedral Foam SDF
// Evaluates the exact 3D Euclidean distance to all 24 struts of the FCC
// rhombic dodecahedron open-cell foam.
// ----------------------------------------------------------------------------
float sdRhombicFoam( in vec3 p, in float scale, in float r )
{
    vec3 q = opRepFCC(p, scale);
    vec3 w = abs(q);
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; }
    if (w.y < w.z) { float t = w.y; w.y = w.z; w.z = t; }
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; }

    // Edge A: along w.x = w.y = 0.25*s, w.z in [0, 0.25*s]
    float tA = clamp(w.z, 0.0, 0.25 * scale);
    vec3 dA = w - vec3(0.25 * scale, 0.25 * scale, tA);

    // Edge B: along w.y = w.z = t, w.x = 0.5*s - t for t in [0, 0.25*s]
    float tB = clamp((w.y + w.z - w.x + 0.5 * scale) / 3.0, 0.0, 0.25 * scale);
    vec3 dB = w - vec3(0.5 * scale - tB, tB, tB);

    return sqrt(min(dot(dA, dA), dot(dB, dB))) - r;
}

// ----------------------------------------------------------------------------
// sdgRhombicFoam — Exact 3D Distance & Analytical Unit Normal of Rhombic Foam
// Returns vec4(d, normal.xyz).
// ----------------------------------------------------------------------------
vec4 sdgRhombicFoam( in vec3 p, in float scale, in float r )
{
    vec3 q = opRepFCC(p, scale);
    vec3 w = abs(q);
    ivec3 idx;
    sdfGemsSort3Desc(w, idx);

    float tA = clamp(w.z, 0.0, 0.25 * scale);
    vec3 dA = w - vec3(0.25 * scale, 0.25 * scale, tA);
    float lenA_sq = dot(dA, dA);

    float tB = clamp((w.y + w.z - w.x + 0.5 * scale) / 3.0, 0.0, 0.25 * scale);
    vec3 dB = w - vec3(0.5 * scale - tB, tB, tB);
    float lenB_sq = dot(dB, dB);

    bool useA = (lenA_sq <= lenB_sq);
    vec3 d_win = useA ? dA : dB;
    float dist = sqrt(useA ? lenA_sq : lenB_sq);

    vec3 gw = (dist > 1e-7 * scale)
        ? (d_win / dist)
        : (useA ? vec3(0.70710678, -0.70710678, 0.0) : vec3(0.81649658, -0.40824829, -0.40824829));
    vec3 grad = sdfGemsUnpermuteGrad3(gw, idx, q);

    return vec4(dist - r, grad);
}

// ----------------------------------------------------------------------------
// sdOctetTruss — 0-Loop, 0-Clamp, 1-Sqrt Exact 3D FCC Octet Truss Metamaterial
// Evaluates the exact 3D Euclidean distance to Buckminster Fuller's (1961)
// infinite 12-connected FCC Octet Truss lattice: sqrt(0.5*(w1 - w2)^2 + w3^2) - r.
// ----------------------------------------------------------------------------
float sdOctetTruss( in vec3 p, in float scale, in float r )
{
    vec3 q = opRepFCC(p, scale);
    vec3 w = abs(q);
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; }
    if (w.y < w.z) { float t = w.y; w.y = w.z; w.z = t; }
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; }

    float diff = w.x - w.y;
    return sqrt(0.5 * diff * diff + w.z * w.z) - r;
}

// ----------------------------------------------------------------------------
// sdgOctetTruss — Exact 3D Distance & Analytical Unit Normal of FCC Octet Truss
// Returns vec4(d, normal.xyz).
// ----------------------------------------------------------------------------
vec4 sdgOctetTruss( in vec3 p, in float scale, in float r )
{
    vec3 q = opRepFCC(p, scale);
    vec3 w = abs(q);
    ivec3 idx;
    sdfGemsSort3Desc(w, idx);

    float diff = w.x - w.y;
    float dist = sqrt(0.5 * diff * diff + w.z * w.z);

    vec3 gw = (dist > 1e-7 * scale)
        ? vec3(0.5 * diff, -0.5 * diff, w.z) / dist
        : vec3(0.70710678, -0.70710678, 0.0);
    vec3 grad = sdfGemsUnpermuteGrad3(gw, idx, q);

    return vec4(dist - r, grad);
}

// ----------------------------------------------------------------------------
// sdBCCTruss — 0-Loop, 0-Clamp, 1-Sqrt Exact 3D BCC Body-Diagonal Truss
// Evaluates the exact 3D Euclidean distance to the 8-connected BCC body-diagonal
// isotropic truss network.
// ----------------------------------------------------------------------------
float sdBCCTruss( in vec3 p, in float scale, in float r )
{
    vec3 q = opRepBCC(p, scale);
    vec3 w = abs(q);
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; }
    if (w.y < w.z) { float t = w.y; w.y = w.z; w.z = t; }
    if (w.x < w.y) { float t = w.x; w.x = w.y; w.y = t; }

    // Strut 1: through (0,0,0) along (1,1,1)
    float m1 = (w.x + w.y + w.z) / 3.0;
    vec3 d1 = w - m1;

    // Strut 2: through (scale, 0, 0) along (-1, 1, 1) -> v = (scale - w.x, w.y, w.z) along (1,1,1)
    vec3 v = vec3(scale - w.x, w.y, w.z);
    float m2 = (v.x + v.y + v.z) / 3.0;
    vec3 d2 = v - m2;

    return sqrt(min(dot(d1, d1), dot(d2, d2))) - r;
}

// ----------------------------------------------------------------------------
// sdgBCCTruss — Exact 3D Distance & Analytical Unit Normal of BCC Truss
// Returns vec4(d, normal.xyz).
// ----------------------------------------------------------------------------
vec4 sdgBCCTruss( in vec3 p, in float scale, in float r )
{
    vec3 q = opRepBCC(p, scale);
    vec3 w = abs(q);
    ivec3 idx;
    sdfGemsSort3Desc(w, idx);

    float m1 = (w.x + w.y + w.z) / 3.0;
    vec3 d1 = w - m1;
    float len1_sq = dot(d1, d1);

    vec3 v = vec3(scale - w.x, w.y, w.z);
    float m2 = (v.x + v.y + v.z) / 3.0;
    vec3 d2 = v - m2;
    float len2_sq = dot(d2, d2);

    bool use1 = (len1_sq <= len2_sq);
    float dist = sqrt(use1 ? len1_sq : len2_sq);

    vec3 gw;
    if (use1) {
        gw = (dist > 1e-7 * scale) ? (d1 / dist) : vec3(0.70710678, -0.70710678, 0.0);
    } else {
        vec3 gv = (dist > 1e-7 * scale) ? (d2 / dist) : vec3(0.70710678, -0.70710678, 0.0);
        gw = vec3(-gv.x, gv.y, gv.z);
    }
    vec3 grad = sdfGemsUnpermuteGrad3(gw, idx, q);

    return vec4(dist - r, grad);
}

#endif // SDF_GEMS_HAS_KELVIN_FOAM
