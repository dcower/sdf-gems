// ============================================================================
// sdMeissnerTetrahedronV, sdMeissnerTetrahedronF, sdgMeissnerTetrahedronV,
// sdgMeissnerTetrahedronF — Exact SDFs of BOTH Meissner Tetrahedra
// (the classic 3D bodies of constant width obtained from the Reuleaux
// tetrahedron by smoothing three of its six circular edges)
//
// @author Claude Opus 5.5 (sdf-gems, 2026)
// @exact  Yes (0 trig; interior: 1-2 sqrt; exterior: 3-6 SFU incl. one extra
//         farthest-point evaluation for the sharp-feature test)
//
// Constant-Width Farthest-Point Duality:
//   For ANY convex body M of constant width s,
//       SDF_M(p) = max_{o in M} |p - o| - s        (exact, inside and outside)
//   Outside: nearest point q, outward normal n, opposite point o = q - s n lies
//   in M, |p - o| = d + s, and diam(M) = s bounds every other point.
//   Inside: M = ∩_{x in M} B(x, s) gives B(p, s - max|p - o|) ⊂ M, and the
//   opposite point of the nearest boundary point attains the bound.
//
// Meissner bodies as a ball-intersection over their sharp set:
//   M = ∩_{y in Y} B(y, s), Y = 4 tetrahedron vertices (edge s) ∪ 3 SHARP
//   Reuleaux edge arcs. The arc between va, vb is the short 70.53° arc of the
//   circle S(vc, s) ∩ S(vd, s) (center (vc + vd) / 2, radius s * sqrt(3) / 2).
//     Type V (vertex): the 3 smoothed edges meet at v0 -> sharp arcs v1v2v3.
//     Type F (face):   the 3 sharp arcs meet at v0 (v0v1, v0v2, v0v3).
//
// Exact SDF with G(p) = max_{y in Y} |p - y|:
//   * Every farthest point of M from p is a sharp point unless the nearest
//     point of M to p is itself sharp. Hence d = G - s exactly in the whole
//     interior and wherever the nearest boundary point is smooth.
//   * Otherwise the nearest point IS the nearest sharp point y. Exact test:
//     the opposite point o = y - s (p - y) / |p - y| lies in M  (G(o) <= s)
//     <=> d = |p - y| (since then G-s <= d <= |p - y| = max|p - o| - s).
//
// C3v folding: v0 = (1,1,1)k, v1 = (1,-1,-1)k, v2 = (-1,1,-1)k,
// v3 = (-1,-1,1)k, k = s / (2 sqrt 2). Both bodies are invariant under all
// coordinate permutations, so sorting x >= y >= z folds R^3 into one chamber.
//   Farthest Y point lies in the opposite chamber z >= y >= x:
//     V: v0, v3 or the half-arc v2v3 (center (k,0,0), plane y + z = 0)
//     F: v0, v3 or the arc v0v3      (center (0,0,-k), plane x = y)
//   Nearest sharp point lies in the same chamber x >= y >= z:
//     V: v0 or the half-arc v1v2 (center (0,0,k), plane x + y = 0) -> v1
//     F: the arc v0v1             (center (-k,0,0), plane y = z) -> v0 / v1
//   All "antipode / projection inside the 70.53° arc" tests reduce to
//   polynomial comparisons (cos(35.26°) = sqrt(2/3)).
// ============================================================================

#include "../common/common.glsl"

#ifndef SDF_GEMS_SD_MEISSNER_TETRAHEDRON_GLSL
#define SDF_GEMS_SD_MEISSNER_TETRAHEDRON_GLSL

// Relative tolerance on G(o)^2 <= s^2 for the sharp-feature test. A false
// accept only costs a second-order error; a false reject would be first order.
#define SDF_GEMS_MEISSNER_TOL 1e-4

// Sort descending (x >= y >= z): the C3v fold shared by both Meissner bodies.
vec3 sdMeissnerSort( vec3 q )
{
    if( q.x < q.y ) q.xy = q.yx;
    if( q.y < q.z ) q.yz = q.zy;
    if( q.x < q.y ) q.xy = q.yx;
    return q;
}

// Farthest sharp point o* of Y from a SORTED q: returns vec4(o*, |q - o*|^2).
vec4 sdMeissnerFar( vec3 q, float k, float rho, bool faceType )
{
    vec3  v0 = vec3(k), v3 = vec3(-k, -k, k);
    float d0 = dot2(q - v0), d3 = dot2(q - v3);
    vec4  best = (d0 >= d3) ? vec4(v0, d0) : vec4(v3, d3);
    if( faceType )
    {
        // Arc v0v3: center (0,0,-k), plane x = y; in-plane coords (A/sqrt2, B).
        float A = q.x + q.y, B = q.z + k;
        if( -B >= abs(A) )
        {
            float r  = sqrt(0.5 * A * A + B * B);
            float rp = r + rho, h = q.x - q.y;
            float g2 = 0.5 * h * h + rp * rp;
            vec3  o  = vec3(0.0, 0.0, -k) - (rho / max(r, 1e-30)) * vec3(0.5 * A, 0.5 * A, B);
            best = vec4(o, g2); // the antipode beats both arc endpoints v0, v3
        }
    }
    else
    {
        // Half-arc v2v3 near v3: center (k,0,0), plane y + z = 0.
        float X = q.x - k, T = q.y - q.z;
        if( X >= T )
        {
            float r  = sqrt(X * X + 0.5 * T * T);
            float rp = r + rho, h = q.y + q.z;
            float g2 = 0.5 * h * h + rp * rp;
            vec3  o  = vec3(k, 0.0, 0.0) - (rho / max(r, 1e-30)) * vec3(X, 0.5 * T, -0.5 * T);
            if( g2 > best.w ) best = vec4(o, g2);
        }
    }
    return best;
}

// Nearest sharp point of Y (vertex or arc point) to a SORTED q.
vec3 sdMeissnerNear( vec3 q, float k, float rho, bool faceType )
{
    if( faceType )
    {
        // Arc v0v1: center (-k,0,0), plane y = z, endpoints v0 (B > 0), v1 (B < 0).
        float A = q.x + k, B = q.y + q.z;
        if( A >= abs(B) )
        {
            float r = sqrt(A * A + 0.5 * B * B);
            return vec3(-k, 0.0, 0.0) + (rho / max(r, 1e-30)) * vec3(A, 0.5 * B, 0.5 * B);
        }
        return (B >= 0.0) ? vec3(k) : vec3(k, -k, -k);
    }
    // Half-arc v1v2 near v1: center (0,0,k), plane x + y = 0, or apex v0.
    float A = q.x - q.y, B = q.z - k;
    vec3  y = vec3(k, -k, -k);
    if( -B >= A )
    {
        float r = sqrt(0.5 * A * A + B * B);
        y = vec3(0.0, 0.0, k) + ((r > 1e-30) ? (rho / r) * vec3(0.5 * A, -0.5 * A, B) : vec3(0.0, 0.0, -rho));
    }
    vec3 v0 = vec3(k);
    return (dot2(q - v0) < dot2(q - y)) ? v0 : y;
}

float sdMeissnerTetrahedronImpl( vec3 p, float s, bool faceType )
{
    float k = 0.35355339 * s, rho = 0.86602540 * s;
    vec3  q = sdMeissnerSort(p);
    float G = sqrt(sdMeissnerFar(q, k, rho, faceType).w);
    if( G <= s ) return G - s;                         // whole interior: exact
    vec3  y  = sdMeissnerNear(q, k, rho, faceType);
    vec3  v  = q - y;
    float dy = length(v);
    vec3  o  = y - (s / max(dy, 1e-30)) * v;           // opposite point of y
    float Go = sdMeissnerFar(sdMeissnerSort(o), k, rho, faceType).w;
    return (Go <= s * s * (1.0 + SDF_GEMS_MEISSNER_TOL)) ? dy : G - s;
}

vec4 sdgMeissnerTetrahedronImpl( vec3 p, float s, bool faceType )
{
    float k = 0.35355339 * s, rho = 0.86602540 * s;
    vec3  q = p;
    bool  s1 = (q.x < q.y); if( s1 ) q.xy = q.yx;
    bool  s2 = (q.y < q.z); if( s2 ) q.yz = q.zy;
    bool  s3 = (q.x < q.y); if( s3 ) q.xy = q.yx;

    vec4  f = sdMeissnerFar(q, k, rho, faceType);
    float G = sqrt(f.w);
    float d = G - s;
    vec3  g = (q - f.xyz) / max(G, 1e-30);
    if( G > s )
    {
        vec3  y  = sdMeissnerNear(q, k, rho, faceType);
        vec3  v  = q - y;
        float dy = length(v);
        float iv = 1.0 / max(dy, 1e-30);
        float Go = sdMeissnerFar(sdMeissnerSort(y - (s * iv) * v), k, rho, faceType).w;
        if( Go <= s * s * (1.0 + SDF_GEMS_MEISSNER_TOL) )
        {
            d = dy;
            g = (dy > 1e-30) ? v * iv : normalize(q);
        }
    }
    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4(d, g);
}

// ----------------------------------------------------------------------------
// Public API. s = constant width (= tetrahedron edge length). The C3v axis is
// (1,1,1); v0 = (1,1,1) * s / (2 sqrt 2) is the apex.
//   ...V: smoothed edges meet at the apex v0 (sharp base triangle v1v2v3).
//   ...F: sharp edges meet at the apex v0 (smoothed base triangle).
// sdg* return vec4(distance, unit gradient).
// ----------------------------------------------------------------------------
float sdMeissnerTetrahedronV( vec3 p, float s ) { return sdMeissnerTetrahedronImpl(p, s, false); }
float sdMeissnerTetrahedronF( vec3 p, float s ) { return sdMeissnerTetrahedronImpl(p, s, true); }
vec4  sdgMeissnerTetrahedronV( vec3 p, float s ) { return sdgMeissnerTetrahedronImpl(p, s, false); }
vec4  sdgMeissnerTetrahedronF( vec3 p, float s ) { return sdgMeissnerTetrahedronImpl(p, s, true); }

#endif // SDF_GEMS_SD_MEISSNER_TETRAHEDRON_GLSL
