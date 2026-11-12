// ============================================================================
// lcBounds — Certified two-sided bounds on an SDF from K cached samples
//   lcLower / lcUpper            McShane-Whitney bounds from values            (any 1-Lipschitz d)
//   lcClosestPoint               c = x - d * n                                 (EXACT SDF + unit gradient only)
//   lcUpperClosest               |d(x)| <= min_i |x - c_i|                     (c_i on the surface)
//   lcLowerConvex                d(x) >= max_i n_i . (x - c_i)                 (CONVEX body, exact SDF)
//   lcUpperConvexHull            d(x) <= dist(x, segments c_i c_j)             (CONVEX body, exact SDF)
//   lcInterval / lcIntervalConvex  [lo, hi] combining the above
//   lcRayStepConvex              Newton / ray-plane jump, never oversteps      (CONVEX body, exact SDF, d > 0)
//   lcRaySpanConvex              [t_in, t_out] of a ray against cached planes  (CONVEX body)
//   lcRaySkip                    march a ray through the union of cached empty balls (any 1-Lipschitz lower bound)
//   lcExample*                   4 convex repo primitives + lazy per-part marchers (preview / benchmark)
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  Conservative bounds (certificates), NOT distances. Each function states its validity condition.
//         Measured (FP32 GPU): violations <= 7.2e-7 under the stated conditions with samples taken OUTSIDE (d_i > 0).
//         The convex bounds are violated by up to 1.4 units on non-convex bodies; the closest-point bound fails for
//         non-exact (bound) SDFs. INSIDE samples are only safe away from the medial axis: there the unit normal is
//         ill-conditioned or a fallback, and c_i = x_i - d_i n_i leaves the surface (capsule: half-space bound wrong by
//         up to 1.1 for samples on the axis). Keep d_i > 0 samples for anything that must be certified.
//         See docs/optics/lcBounds.md.
// ============================================================================

#ifndef SDF_GEMS_LC_BOUNDS_GLSL
#define SDF_GEMS_LC_BOUNDS_GLSL

#include "../gradients3d/sdgBox.glsl"
#include "../gradients3d/sdgEllipsoid.glsl"
#include "../gradients3d/sdgCapsule.glsl"
#include "../gradients3d/sdgOctahedron.glsl"

#ifndef LC_K
#define LC_K 8            // cache capacity (samples); every function reads the first n <= LC_K entries
#endif

const float LC_INF = 1e30;

// ---------------------------------------------------------------------------
// 1. Value-only bounds (McShane 1934 / Whitney 1934). Samples s[i] = vec4(x_i, d_i).
//    VALID FOR: any 1-Lipschitz function d (exact SDFs and 1-Lipschitz bounds alike).
//    These are the tightest bounds obtainable from the values alone.
// ---------------------------------------------------------------------------
float lcLower( vec3 x, vec4 s[LC_K], int n )
{
    float lo = -LC_INF;
    for( int i = 0; i < n; i++ ) lo = max(lo, s[i].w - length(x - s[i].xyz));
    return lo;
}

float lcUpper( vec3 x, vec4 s[LC_K], int n )
{
    float hi = LC_INF;
    for( int i = 0; i < n; i++ ) hi = min(hi, s[i].w + length(x - s[i].xyz));
    return hi;
}

// ---------------------------------------------------------------------------
// 2. Closest surface point of a sample: c = x - d * grad d.
//    VALID FOR: exact SDFs with a unit gradient (dg = vec4(d, n) from an sdg* function), inside or outside.
//    FAILS FOR: bound SDFs (c is then not on the surface), and for inside samples on / near the medial axis
//    (normal undefined or FP32-ill-conditioned; error of c is about |d| * 1e-7 / distance to the medial axis).
// ---------------------------------------------------------------------------
vec3 lcClosestPoint( vec3 x, vec4 dg )
{
    return x - dg.x * dg.yzw;
}

// |d(x)| <= min_i |x - c_i|, i.e. -result <= d(x) <= result.
// VALID FOR: any closed set, any x, provided every c_i lies ON the surface.
// (Identical to the semiconcavity bound d^2(x) <= d_i^2 + 2 d_i n_i.(x - x_i) + |x - x_i|^2.)
float lcUpperClosest( vec3 x, vec3 c[LC_K], int n )
{
    float m = LC_INF;
    for( int i = 0; i < n; i++ ) { vec3 v = x - c[i]; m = min(m, dot(v, v)); }
    return sqrt(m);
}

// ---------------------------------------------------------------------------
// 3. Convex bodies: supporting half-spaces. d(x) >= max_i n_i . (x - c_i) for every x (inside too).
//    VALID FOR: exact SDF of a CONVEX body (the SDF is then a convex function; this is its gradient inequality).
//    FAILS FOR: non-convex bodies (torus, link, frame ...) and for unions of convex parts taken as one body.
//    For a union of convex parts keep one cache per part and take the min over parts of this bound.
// ---------------------------------------------------------------------------
float lcLowerConvex( vec3 x, vec3 c[LC_K], vec3 nrm[LC_K], int n )
{
    float lo = -LC_INF;
    for( int i = 0; i < n; i++ ) lo = max(lo, dot(nrm[i], x - c[i]));
    return lo;
}

// Convex bodies: every segment c_i c_j lies inside the body, so d(x) <= distance to the nearest such segment.
// VALID FOR: exact SDF of a CONVEX body. Never larger than lcUpperClosest. O(n^2).
float lcUpperConvexHull( vec3 x, vec3 c[LC_K], int n )
{
    float m = LC_INF;
    for( int i = 0; i < n; i++ )
    {
        vec3 v = x - c[i];
        m = min(m, dot(v, v));
        for( int j = i + 1; j < n; j++ )
        {
            vec3  e  = c[j] - c[i];
            float ee = dot(e, e);
            float h  = (ee > 0.0) ? clamp(dot(v, e) / ee, 0.0, 1.0) : 0.0;
            vec3  w  = v - h * e;
            m = min(m, dot(w, w));
        }
    }
    return sqrt(m);
}

// [lo, hi] for a general exact SDF: McShane-Whitney + closest points.
// (Whether this interval is the tightest possible from the data was not tested.)
vec2 lcInterval( vec3 x, vec4 s[LC_K], vec3 c[LC_K], int n )
{
    float u = lcUpperClosest(x, c, n);
    return vec2(max(lcLower(x, s, n), -u), min(lcUpper(x, s, n), u));
}

// [lo, hi] for the exact SDF of a convex body.
vec2 lcIntervalConvex( vec3 x, vec4 s[LC_K], vec3 c[LC_K], vec3 nrm[LC_K], int n )
{
    // Only the closest points themselves give the two-sided |d| <= |x - c_i|. The chords give an UPPER bound only
    // (a chord interior lies inside the body): -lcUpperConvexHull is not a lower bound (measured violation 0.34).
    float lo = max(max(lcLower(x, s, n), lcLowerConvex(x, c, nrm, n)), -lcUpperClosest(x, c, n));
    return vec2(lo, min(lcUpper(x, s, n), lcUpperConvexHull(x, c, n)));
}

// ---------------------------------------------------------------------------
// 4. Rays.
// ---------------------------------------------------------------------------
// Certified step for a ray x(t) = ro + t * rd (|rd| = 1) from a sample (d > 0, unit normal nrm) of a CONVEX body:
// the ray-plane jump d / (-nrm . rd), which is also the Newton step of g(t) = d(x(t)). Always >= d.
// Returns LC_INF when nrm . rd >= 0: the ray never enters the supporting half-space, so it misses the body.
// VALID FOR: exact SDF of a convex body, d > 0. On non-convex bodies it oversteps (measured: 1.8-3.2% of all pixels wrong,
// march passes up to 0.19 inside the surface).
float lcRayStepConvex( float d, vec3 nrm, vec3 rd )
{
    float g = dot(nrm, rd);
    return (g < 0.0) ? min(d / (-g), LC_INF) : LC_INF;
}

// Parameter span [t_in, t_out] of the ray inside all cached supporting half-spaces of ONE convex body.
// The body cannot be hit before t_in or after t_out; t_in > t_out certifies a miss with zero SDF evaluations.
// The planes are ray-independent, so they may come from other pixels, other rays or earlier frames.
// VALID FOR: exact SDF of a convex body (static, or planes transformed with it).
vec2 lcRaySpanConvex( vec3 ro, vec3 rd, vec3 c[LC_K], vec3 nrm[LC_K], int n )
{
    float tin = 0.0, tout = LC_INF;
    for( int i = 0; i < n; i++ )
    {
        float a = dot(nrm[i], ro - c[i]), b = dot(nrm[i], rd);
        if( b < 0.0 )      tin  = max(tin, min(a / (-b), LC_INF));
        else if( a > 0.0 ) return vec2(LC_INF, -LC_INF);          // outside and not approaching: miss
        else if( b > 0.0 ) tout = min(tout, -a / b);
    }
    return vec2(tin, tout);
}

// Advance t along the ray while x(t) lies in the union of the cached empty balls B(x_i, d_i), d_i > 0.
// Returns the largest t' >= t such that [t, t') is certified free of surface.
// VALID FOR: samples of any 1-Lipschitz lower bound of the distance (exact or conservative SDFs), any shape.
// Along the ray that produced the samples this never beats the last ball; it only helps for OTHER rays.
// CLAMP the result to tmax and do not cache far-field samples: a ball is only as empty as the FP32 absolute error
// of d_i allows (about 6e-8 * d_i). Unclamped pixel-to-pixel chaining reached d_i ~ 1e9 and skipped the whole scene.
float lcRaySkip( vec3 ro, vec3 rd, float t, vec4 s[LC_K], int n )
{
    for( int pass = 0; pass < n; pass++ )
    {
        float tn = t;
        for( int i = 0; i < n; i++ )
        {
            vec3  oc = s[i].xyz - ro;
            float b  = dot(oc, rd);
            float h  = b * b - dot(oc, oc) + s[i].w * abs(s[i].w);
            if( s[i].w > 0.0 && h > 0.0 )
            {
                h = sqrt(h);
                if( b - h <= t && b + h > tn ) tn = b + h;
            }
        }
        if( tn <= t ) break;
        t = tn;
    }
    return t;
}

// ---------------------------------------------------------------------------
// 5. Example: a union of 4 convex repo primitives marched with per-part certificates.
//    Convexity is lost under union, but each part keeps its own half-space bound.
// ---------------------------------------------------------------------------
vec4 sdgRoundBox( in vec3 p, in vec3 b, in float r );
vec4 sdgEllipsoidExact( vec3 p, vec3 r );
vec4 sdgCapsule( vec3 p, vec3 a, vec3 b, float r );
vec4 sdgOctahedron( vec3 p, float s );

#define LC_EX_PARTS 4

// Value + unit gradient of convex part j (all exact SDFs).
vec4 lcExamplePart( vec3 p, int j )
{
    if( j == 0 ) return sdgRoundBox(p - vec3(-0.62, -0.36, -0.15), vec3(0.34, 0.32, 0.30), 0.06);
    if( j == 1 ) return sdgEllipsoidExact(p - vec3(0.30, -0.42, 0.55), vec3(0.46, 0.26, 0.32));
    if( j == 2 ) return sdgCapsule(p, vec3(0.50, -0.50, -0.55), vec3(0.20, 0.42, -0.30), 0.17);
    return sdgOctahedron(p - vec3(-0.35, 0.33, -0.05), 0.36);
}

// Exact SDF + gradient of the union (4 part evaluations).
vec4 lcExampleScene( vec3 p )
{
    vec4 r = lcExamplePart(p, 0);
    for( int j = 1; j < LC_EX_PARTS; j++ ) { vec4 q = lcExamplePart(p, j); if( q.x < r.x ) r = q; }
    return r;
}

// Lazy per-part marcher. te[j] is a certified lower bound on the ray parameter at which part j can be hit
// (tmin, or lcRaySpanConvex(...).x from cached planes, or LC_INF for a part certified missed).
// Each iteration evaluates ONLY the part with the smallest bound, at that bound.
//   convex = true : te[j] advances by the half-space jump (requires convex parts)
//   convex = false: te[j] advances by d (plain Lipschitz ball; valid for any exact or conservative part SDF)
// Returns vec4(t, hit part index or -1, part evaluations, 0).
vec4 lcExampleMarchFrom( vec3 ro, vec3 rd, float te[LC_EX_PARTS], float tmax, float eps, bool convex )
{
    float evals = 0.0;
    for( int it = 0; it < 96 * LC_EX_PARTS; it++ )
    {
        int j = 0;
        for( int k = 1; k < LC_EX_PARTS; k++ ) if( te[k] < te[j] ) j = k;
        float t = te[j];
        if( t > tmax ) break;
        vec4 dg = lcExamplePart(ro + t * rd, j);
        evals += 1.0;
        if( dg.x < eps ) return vec4(t, float(j), evals, 0.0);
        te[j] = t + (convex ? lcRayStepConvex(dg.x, dg.yzw, rd) : dg.x);
    }
    return vec4(tmax, -1.0, evals, 0.0);
}

vec4 lcExampleMarchConvex( vec3 ro, vec3 rd, float tmin, float tmax, float eps )
{
    float te[LC_EX_PARTS];
    for( int j = 0; j < LC_EX_PARTS; j++ ) te[j] = tmin;
    return lcExampleMarchFrom(ro, rd, te, tmax, eps, true);
}

vec4 lcExampleMarchLazy( vec3 ro, vec3 rd, float tmin, float tmax, float eps )
{
    float te[LC_EX_PARTS];
    for( int j = 0; j < LC_EX_PARTS; j++ ) te[j] = tmin;
    return lcExampleMarchFrom(ro, rd, te, tmax, eps, false);
}

// Baseline: plain sphere tracing of the union (LC_EX_PARTS part evaluations per step).
vec4 lcExampleMarchPlain( vec3 ro, vec3 rd, float tmin, float tmax, float eps )
{
    float t = tmin, evals = 0.0;
    for( int it = 0; it < 256; it++ )
    {
        if( t > tmax ) break;
        float d = lcExampleScene(ro + t * rd).x;
        evals += float(LC_EX_PARTS);
        if( d < eps ) return vec4(t, 0.0, evals, 0.0);
        t += d;
    }
    return vec4(tmax, -1.0, evals, 0.0);
}

#endif // SDF_GEMS_LC_BOUNDS_GLSL
