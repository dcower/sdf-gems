// ============================================================================
// suSupport — Legendre / support-function duality for convex SDFs
//
// A convex body K is determined by its support function h(n) = max_{x in K} n.x.
// The su*() functions return the support POINT x(n) = argmax_{x in K} n.x
// (= grad h(n)), so h(n) = dot(n, x(n)). In the dual:
//   Minkowski sum      x_{A (+) B}(n) = x_A(n) + x_B(n)
//   rotation / shift   x_{R K + t}(n) = R * x_K(n * R) + t        (n * R = R^T n)
//   offset by r        x(n) + r n                                 (n unit)
//   convex hull        the summand with the larger dot(n, x)      (suHull)
//
// For EVERY unit n and every p:   dot(n, p - x(n)) <= sd_K(p)   (a supporting plane).
// So every iterate of any search over n is a conservative sphere-tracing step.
// The search used here is GJK (Gilbert-Johnson-Keerthi) restricted to simplices of at
// most three support points (enough for an exterior point): the closest point v of the
// simplex to p supplies the next n = normalize(p - v).
// It also yields an UPPER bound |p - v| >= sd_K(p) for exterior p, i.e. a bracket.
//
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @source Support mappings + GJK: Gilbert, Johnson & Keerthi (1988); van den Bergen (1999,
//         2004: GJK ray cast by supporting-plane clipping). Nothing here is a new algorithm;
//         see docs/duality/suSupport.md for what was measured and what is claimed.
// @exact  su*: exact support points (verified against the repo SDFs).
//         sd*/sdg* built on SuState: Conservative bound (measured: never above the true
//         signed distance by more than FP32 rounding, 2.2e-7, for any iteration count),
//         converging to the exact exterior distance. NOT exact to 1e-5 at the default
//         SU_ITERS = 16: 96.4% - 100% of exterior points are within 1e-5, 99.86% - 100%
//         within 1e-4, worst 3.4e-4. Within ~1e-4 of the surface FP32 cannot resolve the
//         direction p - v, and the bound can be loose by up to 5e-4 there.
//         Interior: sdMinkowskiBoxBox is exact (15 separating axes); the generic interior
//         value is a conservative bound (never above the true value, i.e. the depth is
//         over-estimated) from a local search: Approximate, within 1e-4 for 94.6% - 99.7%
//         of interior points, worst 2.7e-2 (measured in the doc).
// ============================================================================

#ifndef SDF_GEMS_HAS_SU_SUPPORT
#define SDF_GEMS_HAS_SU_SUPPORT

#ifndef SU_ITERS
#define SU_ITERS 16          // iteration budget of the cold-start sd*() examples
#endif
#ifndef SU_TOL
#define SU_TOL 2e-6          // early exit when (upper - lower) <= SU_TOL * (1 + upper)
#endif

// ---------------------------------------------------------------------------
// 1. Support points. n must be unit length for the curved bodies (sphere, round
//    box, capsule, ellipsoid is scale-free, torus hull); polytopes, disk, cylinder
//    and cones only use the direction of n. Ties (n normal to a face / edge) return
//    one valid support point.
// ---------------------------------------------------------------------------
vec3 suSgn( in vec3 n )
{
    return vec3((n.x < 0.0) ? -1.0 : 1.0, (n.y < 0.0) ? -1.0 : 1.0, (n.z < 0.0) ? -1.0 : 1.0);
}

vec3 suSphere( in vec3 n, in float r ) { return r * n; }

vec3 suBox( in vec3 n, in vec3 b ) { return suSgn(n) * b; }

// Box of half-extents b offset by r (outer half-extents b + r).
vec3 suRoundBox( in vec3 n, in vec3 b, in float r ) { return suSgn(n) * b + r * n; }

vec3 suSegment( in vec3 n, in vec3 a, in vec3 b ) { return (dot(n, b - a) > 0.0) ? b : a; }

vec3 suCapsule( in vec3 n, in vec3 a, in vec3 b, in float r ) { return suSegment(n, a, b) + r * n; }

// Disk of radius r in the plane y = 0.
vec3 suDisk( in vec3 n, in float r )
{
    float l2 = dot(n.xz, n.xz);
    vec2  u  = (l2 > 1e-24) ? n.xz * inversesqrt(l2) : vec2(1.0, 0.0);
    return vec3(r * u.x, 0.0, r * u.y);
}

// Same convention as sdCappedCylinder(p, h, r): axis y, half-height h, radius r.
vec3 suCappedCylinder( in vec3 n, in float h, in float r )
{
    vec3 x = suDisk(n, r);
    x.y = (n.y < 0.0) ? -h : h;
    return x;
}

// Same convention as sdCone(p, c, h): apex at the origin, c = (sin, cos) of the half
// angle, base disk of radius h * c.x / c.y at y = -h.
vec3 suCone( in vec3 n, in vec2 c, in float h )
{
    float R = h * c.x / c.y;
    vec3  x = suDisk(n, R);
    x.y = -h;
    return (R * length(n.xz) - h * n.y > 0.0) ? x : vec3(0.0);
}

// Same convention as sdCappedCone(p, h, r1, r2): radius r1 at y = -h, r2 at y = +h.
vec3 suCappedCone( in vec3 n, in float h, in float r1, in float r2 )
{
    bool lo = (r1 - r2) * length(n.xz) > 2.0 * h * n.y;
    vec3 x  = suDisk(n, lo ? r1 : r2);
    x.y = lo ? -h : h;
    return x;
}

// Ellipsoid with semi-axes r: x(n) = r^2 n / |r n|.
vec3 suEllipsoid( in vec3 n, in vec3 r )
{
    vec3 q = r * n;
    return r * q * inversesqrt(max(dot(q, q), 1e-30));
}

// Convex hull of the torus sdTorus(p, t): a disk of radius t.x offset by t.y.
vec3 suTorusHull( in vec3 n, in vec2 t ) { return suDisk(n, t.x) + t.y * n; }

// Same convention as sdTetrahedron(p, s): vertices s * (+-1, +-1, +-1) with an even
// number of minus signs.
vec3 suTetrahedron( in vec3 n, in float s )
{
    vec3 g = suSgn(n), a = abs(n);
    if( g.x * g.y * g.z < 0.0 )
    {
        if( a.x <= a.y && a.x <= a.z ) g.x = -g.x;
        else if( a.y <= a.z )          g.y = -g.y;
        else                           g.z = -g.z;
    }
    return s * g;
}

// Same convention as sdOctahedron(p, s): vertices at distance s on the axes.
vec3 suOctahedron( in vec3 n, in float s )
{
    vec3 a = abs(n), g = suSgn(n);
    if( a.x >= a.y && a.x >= a.z ) return vec3(s * g.x, 0.0, 0.0);
    if( a.y >= a.z )               return vec3(0.0, s * g.y, 0.0);
    return vec3(0.0, 0.0, s * g.z);
}

// Same convention as sdIcosahedron(p, r): inradius r.
vec3 suIcosahedron( in vec3 n, in float r )
{
    const float PHI = 1.61803399;
    vec3 a = abs(n), g = suSgn(n) * (0.66158454 * r);
    vec3 c = vec3(a.z + PHI * a.y, a.x + PHI * a.z, a.y + PHI * a.x);
    if( c.x >= c.y && c.x >= c.z ) return vec3(0.0, PHI * g.y, g.z);
    if( c.y >= c.z )               return vec3(g.x, 0.0, PHI * g.z);
    return vec3(PHI * g.x, g.y, 0.0);
}

// Same convention as sdDodecahedron(p, r): inradius r.
vec3 suDodecahedron( in vec3 n, in float r )
{
    const float PHI = 1.61803399, IPH = 0.61803399;
    vec3 a = abs(n), g = suSgn(n) * (0.72654253 * r);
    float c0 = a.x + a.y + a.z;
    vec3  c  = vec3(IPH * a.y + PHI * a.z, IPH * a.z + PHI * a.x, IPH * a.x + PHI * a.y);
    float m  = max(c.x, max(c.y, c.z));
    if( c0 >= m )  return g;
    if( c.x >= m ) return vec3(0.0, IPH * g.y, PHI * g.z);
    if( c.y >= m ) return vec3(PHI * g.x, 0.0, IPH * g.z);
    return vec3(IPH * g.x, PHI * g.y, 0.0);
}

// Convex hull of two bodies: keep the support point that is further along n.
vec3 suHull( in vec3 n, in vec3 xa, in vec3 xb ) { return (dot(n, xa) > dot(n, xb)) ? xa : xb; }

// ---------------------------------------------------------------------------
// 2. Dual solver: triangle-simplex GJK with a certified bracket lo <= sd(p) <= ub (p outside).
// ---------------------------------------------------------------------------
struct SuState
{
    vec3  a, b, c;   // simplex: m affinely independent support points of K
    int   m;         // simplex size: 1 (point), 2 (segment) or 3 (triangle)
    vec3  v;         // closest point of the simplex to the query point (a point of K)
    vec3  n;         // unit direction that produced lo (outward normal estimate)
    vec3  p;         // query point of the last update (for warm starts)
    float lo;        // best lower bound  max_k dot(n_k, p - x(n_k))  <= sd(p)
    float ub;        // |p - v| >= sd(p): upper bound
    vec3  alt;       // direction to try next when the search has stalled (stall == 1)
    int   stall;     // 0: progressing, 1: stalled, try alt, 2: stalled for good (suDone)
};

// Barycentric weights of the point of triangle abc closest to p (Ericson, Real-Time
// Collision Detection 5.1.5). Weights are exactly 0 for vertices that do not support
// the closest feature, and always form a convex combination (also for degenerate
// triangles, where FMA contraction makes the region tests unreliable).
vec3 suTriWeights( in vec3 p, in vec3 a, in vec3 b, in vec3 c )
{
    vec3  ab = b - a, ac = c - a, ap = p - a;
    float d1 = dot(ab, ap), d2 = dot(ac, ap);
    if( d1 <= 0.0 && d2 <= 0.0 ) return vec3(1.0, 0.0, 0.0);
    vec3  bp = p - b;
    float d3 = dot(ab, bp), d4 = dot(ac, bp);
    if( d3 >= 0.0 && d4 <= d3 ) return vec3(0.0, 1.0, 0.0);
    float vc = d1 * d4 - d3 * d2;
    if( vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 )
    {
        float t = d1 / max(d1 - d3, 1e-30);
        return vec3(1.0 - t, t, 0.0);
    }
    vec3  cp = p - c;
    float d5 = dot(ab, cp), d6 = dot(ac, cp);
    if( d6 >= 0.0 && d5 <= d6 ) return vec3(0.0, 0.0, 1.0);
    float vb = d5 * d2 - d1 * d6;
    if( vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 )
    {
        float t = d2 / max(d2 - d6, 1e-30);
        return vec3(1.0 - t, 0.0, t);
    }
    float va = d3 * d6 - d5 * d4;
    if( va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 )
    {
        float t = (d4 - d3) / max((d4 - d3) + (d5 - d6), 1e-30);
        return vec3(0.0, 1.0 - t, t);
    }
    float den = 1.0 / max(va + vb + vc, 1e-30);
    float v   = clamp(vb * den, 0.0, 1.0);
    float w   = clamp(vc * den, 0.0, 1.0 - v);
    return vec3(1.0 - v - w, v, w);
}

// Replace the simplex by the vertices of (x0, x1, x2) that carry weight. Weights below
// SU_WEIGHT_EPS are rounding noise: keeping such a vertex would make suDir() use the
// normal of a face the closest point does not really lie inside.
#ifndef SU_WEIGHT_EPS
#define SU_WEIGHT_EPS 1e-6
#endif
void suAdopt( inout SuState s, in vec3 x0, in vec3 x1, in vec3 x2, in vec3 u )
{
    u  = mix(u, vec3(0.0), lessThan(u, vec3(SU_WEIGHT_EPS)));
    u /= (u.x + u.y + u.z);
    s.v = u.x * x0 + u.y * x1 + u.z * x2;
    int m = 0;
    if( u.x > 0.0 ) { s.a = x0; m = 1; }
    if( u.y > 0.0 ) { if( m == 0 ) s.a = x1; else s.b = x1; m++; }
    if( u.z > 0.0 ) { if( m == 0 ) s.a = x2; else if( m == 1 ) s.b = x2; else s.c = x2; m++; }
    s.m = m;
}

vec3 suSegWeights( in vec3 p, in vec3 a, in vec3 b )
{
    vec3  e = b - a;
    float t = clamp(dot(p - a, e) / max(dot(e, e), 1e-30), 0.0, 1.0);
    return vec3(1.0 - t, t, 0.0);
}

// First search direction: away from the origin (assumed to be in or near K).
vec3 suDir0( in vec3 p )
{
    float l2 = dot(p, p);
    return (l2 > 1e-24) ? p * inversesqrt(l2) : vec3(0.0, 1.0, 0.0);
}

// Cold start from one support point x0 = x(n0), n0 unit.
SuState suInit( in vec3 p, in vec3 n0, in vec3 x0 )
{
    SuState s;
    s.a = x0; s.b = x0; s.c = x0; s.m = 1; s.v = x0;
    s.n = n0; s.p = p; s.alt = n0; s.stall = 0;
    s.lo = dot(n0, p - x0);
    s.ub = length(p - x0);
    return s;
}

// Warm start at a new query point: the simplex stays inside K, and the old supporting
// plane is still a valid bound, lo' = lo + dot(n, p' - p), with no support evaluation.
void suRestart( inout SuState s, in vec3 p )
{
    s.lo += dot(s.n, p - s.p);
    s.p   = p;
    s.stall = 0;
    if( s.m == 3 )      suAdopt(s, s.a, s.b, s.c, suTriWeights(p, s.a, s.b, s.c));
    else if( s.m == 2 ) suAdopt(s, s.a, s.b, s.b, suSegWeights(p, s.a, s.b));
    else                s.v = s.a;
    s.ub  = length(p - s.v);
}

// Direction to evaluate the support point in next: normalize(p - v), computed from the
// closest FEATURE (face normal / edge-orthogonal part) instead of the difference p - v,
// which loses all its digits in FP32 when p is within ~1e-4 of the surface.
vec3 suDir( in SuState s, in vec3 p )
{
    if( s.stall == 1 ) return s.alt;
    vec3 r = p - s.a;
    if( s.m == 3 )
    {
        vec3  g  = cross(s.b - s.a, s.c - s.a);
        float l2 = dot(g, g);
        if( l2 > 1e-30 )
        {
            float il = inversesqrt(l2);
            float sg = dot(g, r) * il;
            if( abs(sg) < 1e-6 ) sg = dot(g, s.n);
            return g * ((sg < 0.0) ? -il : il);
        }
        r = p - s.v;
    }
    else if( s.m == 2 )
    {
        vec3 e = s.b - s.a;
        r -= e * (dot(r, e) / max(dot(e, e), 1e-30));
    }
    float l2 = dot(r, r);
    return (l2 > 1e-20) ? r * inversesqrt(l2) : s.n;
}

float suDist2( in vec3 p, in vec3 x0, in vec3 x1, in vec3 x2, in vec3 u )
{
    vec3 e = p - (u.x * x0 + u.y * x1 + u.z * x2);
    return dot(e, e);
}

// Accept a candidate feature if it is closer, or if it is no further than FP32 can
// resolve AND p projects strictly inside it. The second case matters: the upper bound
// is second-order insensitive to a tangential error of v, but the direction p - v, and
// with it the lower bound, is first-order sensitive (error * lever / distance).
// tol = FP32 noise of a squared distance: 2 |p - v| * (rounding error of v).
bool suAccept( in float dNew, in float dOld, in vec3 u, in int m, in float tol )
{
    if( dNew < dOld ) return true;
    bool inside = (m == 1) ? (min(u.x, u.y) > SU_WEIGHT_EPS) : (min(u.x, min(u.y, u.z)) > SU_WEIGHT_EPS);
    return inside && dNew <= dOld + tol;
}

// One iteration: w must be the support point x(n) of the body for the unit direction n.
//
// Stall handling. In FP32 the search can reach a fixed point that is not the optimum: when p is
// within rounding of the boundary between two features of the simplex (say an edge and the face
// next to it), the candidate that contains w is rejected or loses w again, the simplex and hence
// the next direction do not change, and every further iteration repeats the same support point.
// The value is still a valid lower bound but can be loose by (lever arm) * (angle error).
// When a step neither raises lo nor changes the simplex, the next direction is therefore taken
// from the full candidate feature through w (its face normal / edge-orthogonal direction) once;
// if that does not help either, suDone() ends the loop instead of burning the iteration budget.
void suStep( inout SuState s, in vec3 p, in vec3 n, in vec3 w )
{
    float lb = dot(n, p - w);
    bool  better = lb > s.lo;
    if( better ) { s.lo = lb; s.n = n; }
    vec3  e0 = p - s.v;
    float d0 = dot(e0, e0);
    float tl = 4e-7 * sqrt(d0 * (dot(p, p) + dot(w, w)));
    float ub0 = s.ub;
    int   m0  = s.m;
    bool  adopted = false;
    vec3  g;         // normal direction of the full candidate feature that contains w
    if( s.m == 1 )
    {
        vec3 u = suSegWeights(p, w, s.a);
        vec3 e = s.a - w;
        g = (p - w) - e * (dot(p - w, e) / max(dot(e, e), 1e-30));
        if( suAccept(suDist2(p, w, s.a, s.a, u), d0, u, 1, tl) ) { suAdopt(s, w, s.a, s.a, u); adopted = true; }
    }
    else if( s.m == 2 )
    {
        vec3 u = suTriWeights(p, w, s.a, s.b);
        g = cross(s.a - w, s.b - w);
        if( suAccept(suDist2(p, w, s.a, s.b, u), d0, u, 2, tl) ) { suAdopt(s, w, s.a, s.b, u); adopted = true; }
    }
    else
    {
        // p is outside the tetrahedron (w, a, b, c) whenever p is outside K, so the closest
        // point lies on one of the three faces that contain the new vertex w.
        vec3  u1 = suTriWeights(p, w, s.a, s.b);
        vec3  u2 = suTriWeights(p, w, s.b, s.c);
        vec3  u3 = suTriWeights(p, w, s.a, s.c);
        float d1 = suDist2(p, w, s.a, s.b, u1);
        float d2 = suDist2(p, w, s.b, s.c, u2);
        float d3 = suDist2(p, w, s.a, s.c, u3);
        if( d1 <= d2 && d1 <= d3 ) { g = cross(s.a - w, s.b - w); if( suAccept(d1, d0, u1, 3, tl) ) { suAdopt(s, w, s.a, s.b, u1); adopted = true; } }
        else if( d2 <= d3 )        { g = cross(s.b - w, s.c - w); if( suAccept(d2, d0, u2, 3, tl) ) { suAdopt(s, w, s.b, s.c, u2); adopted = true; } }
        else                       { g = cross(s.a - w, s.c - w); if( suAccept(d3, d0, u3, 3, tl) ) { suAdopt(s, w, s.a, s.c, u3); adopted = true; } }
    }
    s.ub = length(p - s.v);
    if( better || (adopted && (s.m != m0 || s.ub < ub0)) ) { s.stall = 0; return; }
    float l2 = dot(g, g);
    if( s.stall == 0 && l2 > 1e-30 )
    {
        float il = inversesqrt(l2);
        float sg = dot(g, p - w) * il;
        if( abs(sg) < 1e-6 ) sg = dot(g, s.n);
        s.alt   = g * ((sg < 0.0) ? -il : il);
        s.stall = 1;
    }
    else s.stall = 2;
}

bool suDone( in SuState s ) { return s.stall == 2 || s.ub - s.lo <= SU_TOL * (1.0 + s.ub); }

// Conservative advancement along a ray with unit direction rd (van den Bergen-style
// clipping against the supporting plane): the ray cannot reach K before travelling
// lo / -dot(n, rd) >= lo. Returns a huge step if the plane separates the ray from K.
float suRayStep( in SuState s, in vec3 rd )
{
    float c = -dot(s.n, rd);
    if( s.lo <= 0.0 ) return s.lo;
    return (c > 1e-6 * s.lo) ? s.lo / c : 1e30;
}

// Compose-your-own pattern. SUP is an expression of the unit direction `n` that returns
// the support point of your body, e.g.  suBox(n, b) + suEllipsoid(n, r).
//   float sdMine( vec3 p ) { SuState s; SU_SOLVE(s, p, suBox(n, b) + suEllipsoid(n, r), 16) return s.lo; }
// SU_RESUME continues from a state kept from a nearby point (previous march step).
#define SU_SOLVE(S, P, SUP, K) { vec3 n = suDir0(P); S = suInit(P, n, SUP); for( int su_i = 0; su_i < (K); su_i++ ) { if( suDone(S) ) break; n = suDir(S, P); suStep(S, P, n, SUP); } }
#define SU_RESUME(S, P, SUP, K) { suRestart(S, P); for( int su_i = 0; su_i < (K); su_i++ ) { if( suDone(S) ) break; vec3 n = suDir(S, P); suStep(S, P, n, SUP); } }

// Interior / not-certified-exterior points (lo <= 0). The interior distance
// -min_n (h(n) - n.p) is a NON-convex problem on the sphere, so this is a local search:
//   1. keep the 3 best of 26 fixed directions;
//   2. from each, "exterior probing": step out to q = p + t N (outside K), solve the
//      convex exterior problem at q, and re-use the supporting plane found there at p;
//      N <- that plane's normal; repeat up to SU_PROBES times while it improves.
// Every candidate is a supporting-plane value, so D stays a conservative bound (D <= sd,
// i.e. too negative inside) whatever happens; it equals the interior distance when one
// of the searches ends in the right basin (see the doc for measured rates).
// D, N: in/out (best value and its unit normal). K: GJK iterations per probe.
#ifndef SU_PROBES
#define SU_PROBES 4
#endif
#define SU_INTERIOR(D, N, P, SUP, K) { float su_f0 = -1e30, su_f1 = -1e30, su_f2 = -1e30; vec3 su_n0 = vec3(0.0, 1.0, 0.0), su_n1 = su_n0, su_n2 = su_n0; for( int su_j = 0; su_j < 27; su_j++ ) { if( su_j == 13 ) continue; vec3 n = normalize(vec3(float(su_j % 3) - 1.0, float((su_j / 3) % 3) - 1.0, float(su_j / 9) - 1.0)); float su_v = dot(n, (P) - (SUP)); if( su_v > su_f0 ) { su_f2 = su_f1; su_n2 = su_n1; su_f1 = su_f0; su_n1 = su_n0; su_f0 = su_v; su_n0 = n; } else if( su_v > su_f1 ) { su_f2 = su_f1; su_n2 = su_n1; su_f1 = su_v; su_n1 = n; } else if( su_v > su_f2 ) { su_f2 = su_v; su_n2 = n; } } float su_m; { vec3 n = su_n0; su_m = 0.02 * length(SUP); } for( int su_c = 0; su_c < 3; su_c++ ) { float su_d = (su_c == 0) ? su_f0 : ((su_c == 1) ? su_f1 : su_f2); vec3 su_g = (su_c == 0) ? su_n0 : ((su_c == 1) ? su_n1 : su_n2); for( int su_r = 0; su_r < SU_PROBES; su_r++ ) { vec3 su_q = (P) + (1.25 * max(-su_d, 0.0) + su_m) * su_g; SuState su_s; SU_SOLVE(su_s, su_q, SUP, K) float su_v = su_s.lo + dot(su_s.n, (P) - su_q); if( su_v > su_d ) { su_d = su_v; su_g = su_s.n; } else break; } if( su_d > D ) { D = su_d; N = su_g; } } }

// ---------------------------------------------------------------------------
// 3. Composed examples (none of these has a closed-form SDF, except the box-box
//    interior). Each has: the support point, a cold-start bound, and a warm-start
//    overload taking the state and an iteration budget.
// ---------------------------------------------------------------------------

// 3a. Box b1 (+) rotated box R * b2 (a zonotope-like polytope with up to 30 faces).
vec3 suMinkowskiBoxBox( in vec3 n, in vec3 b1, in mat3 R, in vec3 b2 )
{
    return suBox(n, b1) + R * suBox(n * R, b2);
}

// Exact interior distance: the faces of the sum are normal to the 3 + 3 axes of the
// boxes or to one of the 9 cross products of their edges (the separating-axis set).
float suMinkowskiBoxBoxInterior( in vec3 p, in vec3 b1, in mat3 R, in vec3 b2 )
{
    float d = -1e30;
    for( int i = 0; i < 3; i++ )
    {
        vec3 e = vec3(i == 0, i == 1, i == 2);
        d = max(d, abs(dot(e, p)) - dot(abs(e), b1) - dot(abs(e * R), b2));
        vec3 f = R[i];
        d = max(d, abs(dot(f, p)) - dot(abs(f), b1) - dot(abs(f * R), b2));
        for( int j = 0; j < 3; j++ )
        {
            vec3  g  = cross(e, R[j]);
            float l2 = dot(g, g);
            if( l2 > 1e-8 )
            {
                g *= inversesqrt(l2);
                d = max(d, abs(dot(g, p)) - dot(abs(g), b1) - dot(abs(g * R), b2));
            }
        }
    }
    return d;
}

float sdMinkowskiBoxBox( in vec3 p, in vec3 b1, in mat3 R, in vec3 b2, inout SuState s, in int iters )
{
    SU_RESUME(s, p, suMinkowskiBoxBox(n, b1, R, b2), iters)
    return s.lo;
}

float sdMinkowskiBoxBox( in vec3 p, in vec3 b1, in mat3 R, in vec3 b2 )
{
    float di = suMinkowskiBoxBoxInterior(p, b1, R, b2);
    if( di <= 0.0 ) return di;
    SuState s;
    SU_SOLVE(s, p, suMinkowskiBoxBox(n, b1, R, b2), SU_ITERS)
    return max(s.lo, di);
}

// 3b. Box b (+) rotated capped cylinder (half-height h, radius r).
vec3 suMinkowskiBoxCylinder( in vec3 n, in vec3 b, in mat3 R, in float h, in float r )
{
    return suBox(n, b) + R * suCappedCylinder(n * R, h, r);
}

float sdMinkowskiBoxCylinder( in vec3 p, in vec3 b, in mat3 R, in float h, in float r, inout SuState s, in int iters )
{
    SU_RESUME(s, p, suMinkowskiBoxCylinder(n, b, R, h, r), iters)
    return s.lo;
}

vec4 sdgMinkowskiBoxCylinder( in vec3 p, in vec3 b, in mat3 R, in float h, in float r )
{
    SuState s;
    SU_SOLVE(s, p, suMinkowskiBoxCylinder(n, b, R, h, r), SU_ITERS)
    float d = s.lo;
    vec3  g = s.n;
    if( d <= 0.0 ) SU_INTERIOR(d, g, p, suMinkowskiBoxCylinder(n, b, R, h, r), SU_ITERS)
    return vec4(d, g);
}

float sdMinkowskiBoxCylinder( in vec3 p, in vec3 b, in mat3 R, in float h, in float r )
{
    return sdgMinkowskiBoxCylinder(p, b, R, h, r).x;
}

// 3c. Ellipsoid (semi-axes re) (+) box b.
vec3 suMinkowskiEllipsoidBox( in vec3 n, in vec3 re, in vec3 b )
{
    return suEllipsoid(n, re) + suBox(n, b);
}

float sdMinkowskiEllipsoidBox( in vec3 p, in vec3 re, in vec3 b, inout SuState s, in int iters )
{
    SU_RESUME(s, p, suMinkowskiEllipsoidBox(n, re, b), iters)
    return s.lo;
}

float sdMinkowskiEllipsoidBox( in vec3 p, in vec3 re, in vec3 b )
{
    SuState s;
    SU_SOLVE(s, p, suMinkowskiEllipsoidBox(n, re, b), SU_ITERS)
    float d = s.lo;
    vec3  g = s.n;
    if( d <= 0.0 ) SU_INTERIOR(d, g, p, suMinkowskiEllipsoidBox(n, re, b), SU_ITERS)
    return d;
}

// 3d. Capsule (a, b, ra) (+) cone (c, h). The ball summand is peeled off:
//     K (+) ball(r) has sd = sd_K - r, so the solver only sees segment (+) cone.
vec3 suMinkowskiSegmentCone( in vec3 n, in vec3 a, in vec3 b, in vec2 c, in float h )
{
    return suSegment(n, a, b) + suCone(n, c, h);
}

float sdMinkowskiCapsuleCone( in vec3 p, in vec3 a, in vec3 b, in float ra, in vec2 c, in float h, inout SuState s, in int iters )
{
    SU_RESUME(s, p, suMinkowskiSegmentCone(n, a, b, c, h), iters)
    return s.lo - ra;
}

float sdMinkowskiCapsuleCone( in vec3 p, in vec3 a, in vec3 b, in float ra, in vec2 c, in float h )
{
    SuState s;
    SU_SOLVE(s, p, suMinkowskiSegmentCone(n, a, b, c, h), SU_ITERS)
    float d = s.lo;
    vec3  g = s.n;
    if( d <= 0.0 ) SU_INTERIOR(d, g, p, suMinkowskiSegmentCone(n, a, b, c, h), SU_ITERS)
    return d - ra;
}

// 3e. Convex hull of a sphere (centre cs, radius rs) and a box (centre cb, half-extents b).
vec3 suHullSphereBox( in vec3 n, in vec3 cs, in float rs, in vec3 cb, in vec3 b )
{
    return suHull(n, cs + rs * n, cb + suBox(n, b));
}

float sdHullSphereBox( in vec3 p, in vec3 cs, in float rs, in vec3 cb, in vec3 b, inout SuState s, in int iters )
{
    SU_RESUME(s, p, suHullSphereBox(n, cs, rs, cb, b), iters)
    return s.lo;
}

vec4 sdgHullSphereBox( in vec3 p, in vec3 cs, in float rs, in vec3 cb, in vec3 b )
{
    SuState s;
    vec3 q = p - cb;   // search around the box centre
    SU_SOLVE(s, q, suHullSphereBox(n, cs - cb, rs, vec3(0.0), b), SU_ITERS)
    float d = s.lo;
    vec3  g = s.n;
    if( d <= 0.0 ) SU_INTERIOR(d, g, q, suHullSphereBox(n, cs - cb, rs, vec3(0.0), b), SU_ITERS)
    return vec4(d, g);
}

float sdHullSphereBox( in vec3 p, in vec3 cs, in float rs, in vec3 cb, in vec3 b )
{
    return sdgHullSphereBox(p, cs, rs, cb, b).x;
}

#endif // SDF_GEMS_HAS_SU_SUPPORT
