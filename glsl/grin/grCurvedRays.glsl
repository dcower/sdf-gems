// ============================================================================
// grFishEyeAdvance, grFishEyeCurvature, grFishEyeChordStep, grFishEyeMarch,
// grCircleSphere, grCirclePlane, grFishEyeOpticalDist,
// grLuneburgAdvance, grLuneburgTauStep, grLuneburgExitTau, grLuneburgThrough,
// grLuneburgOpticalPath, grLuneburgMarch
// — Closed-form curved rays in gradient-index (GRIN) media, sphere-traced
//   along the curve with no ODE solver
// @author Claude Fable 5.2 (sdf-gems, 2026). The ray paths themselves are classical
//         (Maxwell 1854, Luneburg 1944); see docs/grin/grCurvedRays.md for prior art
//         and for what was measured.
// @exact  State updates (grFishEyeAdvance, grLuneburgAdvance, grLuneburgThrough,
//         grLuneburgExitTau), analytic hits (grCircleSphere, grCirclePlane) and optical
//         lengths (grFishEyeOpticalDist, grLuneburgOpticalPath) are closed forms,
//         verified against DOP853 integration of the ray ODE (see doc for numbers).
//         grFishEyeChordStep: exact up to FP32 rounding (lands on the unbounding sphere
//         within a relative 3e-5). grLuneburgExitTau is ill-conditioned for grazing entry.
//         grLuneburgTauStep: CONSERVATIVE bound (arc length <= d).
//         grFishEyeMarch / grLuneburgMarch: sphere tracing, converges from the outside
//         when GR_SCENE is a 1-Lipschitz SDF.
// ============================================================================
//
// Ray equation in a medium of index n(x):  d/ds (n dx/ds) = grad n.
//
//   Maxwell fish-eye   n(r) = n0 / (1 + r^2/a^2)       (all of space)
//       every ray is a circle of constant curvature  kappa = |grad_perp ln n|.
//   Luneburg lens      n(r) = sqrt(2 - r^2/R^2), r <= R (n = 1 outside, index matched)
//       in fictitious time tau (dx/dtau = p = n dx/ds, dp/dtau = -x/R^2) the ray is a
//       harmonic oscillator:  x(tau) = x0 cos(tau/R) + R p0 sin(tau/R).
//
// Sphere tracing: if d = SDF(x), the ball |y - x| < d is empty. Any path of arc length
// <= d that starts at x stays in that ball (chord <= arc), so a step never crosses a surface.

#ifndef SDF_GEMS_GR_CURVED_RAYS_GLSL
#define SDF_GEMS_GR_CURVED_RAYS_GLSL

#ifndef GR_MAX_STEPS
#define GR_MAX_STEPS 256
#endif
#ifndef GR_HIT_EPS
#define GR_HIT_EPS 1e-4
#endif

// ---------------------------------------------------------------------------
// Helpers that stay finite as the argument -> 0 (curvature -> 0, straight rays).
float grSinc( in float t )   // sin(t)/t
{
    return abs(t) < 1e-3 ? 1.0 - t * t / 6.0 : sin(t) / t;
}
float grCosc( in float t )   // (1 - cos t)/t^2, written without cancellation
{
    float h = grSinc(0.5 * t);
    return 0.5 * h * h;
}
float grAsinc( in float u )  // asin(u)/u, u in [0, 1]
{
    return u < 1e-3 ? 1.0 + u * u / 6.0 : asin(min(u, 1.0)) / u;
}

// ---------------------------------------------------------------------------
// Small demo scene (exact SDFs) so that the marchers below are complete programs.
// Override with  #define GR_SCENE myScene  before including this file.
float grDemoBox( in vec3 p, in vec3 b )
{
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}
float grDemoScene( in vec3 p )
{
    float d = p.y + 1.2;                                                  // ground plane
    d = min(d, grDemoBox(p - vec3(-0.9, -0.75, -2.4), vec3(0.45)));       // box behind
    d = min(d, length(p - vec3(1.1, -0.7, -2.0)) - 0.5);                  // sphere behind
    d = min(d, length(vec2(length(p.xz - vec2(0.0, 0.1)) - 0.45, p.y + 0.1)) - 0.09); // torus inside the lens
    return d;
}
#ifndef GR_SCENE
#define GR_SCENE grDemoScene
#endif

// ===========================================================================
// Maxwell fish-eye, n(r) = n0 / (1 + r^2/a^2). Rays do not depend on n0.
// ===========================================================================

// Curvature vector dv/ds = grad_perp ln n at (x, v), |v| = 1. Its length is constant
// along a ray (the ray is a circle of radius 1/|kappa|, centre x + kappa/|kappa|^2).
vec3 grFishEyeCurvature( in vec3 x, in vec3 v, in float a )
{
    vec3 g = x * (-2.0 / (a * a + dot(x, x)));
    return g - dot(g, v) * v;
}

// Advance the ray state by ARC LENGTH s (any sign). Exact. Degrades to x += s v as a -> inf.
void grFishEyeAdvance( inout vec3 x, inout vec3 v, in float s, in float a )
{
    vec3 kv = grFishEyeCurvature(x, v, a);
    float t = length(kv) * s;
    float sc = s * grSinc(t);
    x += v * sc + kv * (s * s * grCosc(t));
    v = normalize(v * cos(t) + kv * sc);
}

// Largest safe arc length for clearance d on a circle of curvature k: the arc whose
// CHORD is d (the step lands on the unbounding sphere). Returns -1 when k d >= 2:
// then the whole circular ray lies inside the empty ball and can never hit anything.
float grFishEyeChordStep( in float d, in float k )
{
    float u = 0.5 * k * d;
    return u >= 1.0 ? -1.0 : d * grAsinc(u);
}

// Smallest arc length s > sMin with  P + a sin(k s)/k + beta (1 - cos(k s))/k^2 = 0, or -1.
// (Substituting w = tan(k s/2)/k gives the quadratic (k^2 P + 2 beta) w^2 + 2 a w + P = 0,
//  which becomes the straight-ray equation 2 a w + P = 0, s = 2 w, at k = 0.)
float grCircleRoot( in float P, in float a, in float beta, in float k, in float sMin )
{
    float A = k * k * P + 2.0 * beta;
    float disc = a * a - A * P;
    if( disc < 0.0 ) return -1.0;
    float q = -(a + (a < 0.0 ? -1.0 : 1.0) * sqrt(disc));
    if( q == 0.0 ) return -1.0;
    float kk = max(k, 1e-18);
    const float PI = 3.14159265358979;
    // roots w1 = q/A, w2 = P/q; half angles h = atan(k w) taken in [0, pi). The denominators are
    // made positive first so that a small positive root keeps its relative precision in FP32.
    float h1 = atan(kk * (A < 0.0 ? -q : q), abs(A));
    float h2 = atan(kk * (q < 0.0 ? -P : P), abs(q));
    if( h1 < 0.0 ) h1 += PI;
    if( h2 < 0.0 ) h2 += PI;
    float s1 = 2.0 * h1 / kk, s2 = 2.0 * h2 / kk;
    float s = 1e30;
    if( s1 > sMin ) s = s1;
    if( s2 > sMin && s2 < s ) s = s2;
    return s < 1e29 ? s : -1.0;
}

// First hit (arc length, or -1) of the fish-eye ray from (x, v) with the sphere |y - c| = rad.
// No marching: a circle meets a sphere where it meets the radical plane (a quadratic).
float grCircleSphere( in vec3 x, in vec3 v, in float a, in vec3 c, in float rad, in float sMin )
{
    vec3 kv = grFishEyeCurvature(x, v, a);
    vec3 e = x - c;
    float le = length(e);
    return grCircleRoot(0.5 * (le - rad) * (le + rad), dot(e, v), 1.0 + dot(e, kv), length(kv), sMin);
}

// First hit (arc length, or -1) of the fish-eye ray with the plane dot(nrm, y) = h, |nrm| = 1.
float grCirclePlane( in vec3 x, in vec3 v, in float a, in vec3 nrm, in float h, in float sMin )
{
    vec3 kv = grFishEyeCurvature(x, v, a);
    return grCircleRoot(dot(nrm, x) - h, dot(nrm, v), dot(nrm, kv), length(kv), sMin);
}

// Optical distance (eikonal, |grad T| = n) between two points: the fish-eye optical metric
// is a round 3-sphere of radius n0 a / 2, so T = n0 a asin(chord). T - T0 is the optical
// SDF of the geodesic ball of optical radius T0 about y (a Euclidean ball, off centre).
float grFishEyeOpticalDist( in vec3 x, in vec3 y, in float a, in float n0 )
{
    float a2 = a * a;
    float q = a * length(x - y) * inversesqrt((a2 + dot(x, x)) * (a2 + dot(y, y)));
    return n0 * a * asin(min(q, 1.0));
}

// Sphere tracing along the circular ray. Returns the ARC LENGTH to the hit (x, v are left
// at the hit point, v = view tangent there), or -1 on a miss.
float grFishEyeMarch( inout vec3 x, inout vec3 v, in float a, in float sMax, out int steps )
{
    float k = length(grFishEyeCurvature(x, v, a));
    sMax = min(sMax, 6.2831853 / max(k, 1e-18));   // the ray is a closed loop
    float sTot = 0.0;
    steps = 0;
    for( int i = 0; i < GR_MAX_STEPS; i++ )
    {
        float d = GR_SCENE(x);
        if( d < GR_HIT_EPS ) return sTot;
        float s = grFishEyeChordStep(d, k);
        if( s < 0.0 ) return -1.0;                  // whole circle is inside the empty ball
        grFishEyeAdvance(x, v, s, a);
        sTot += s;
        steps = i + 1;
        if( sTot > sMax ) return -1.0;
    }
    return -1.0;
}

// ===========================================================================
// Luneburg lens of radius R at the origin, n(r) = sqrt(2 - r^2/R^2) for r <= R.
// n = 1 at r = R: the surface is index matched (no refraction, no Fresnel reflection).
// ===========================================================================

// Advance by FICTITIOUS TIME tau (d tau = ds / n) inside the lens. Exact.
void grLuneburgAdvance( inout vec3 x, inout vec3 v, in float tau, in float R )
{
    vec3 p = v * sqrt(max(2.0 - dot(x, x) / (R * R), 0.0));   // optical momentum n v
    float th = tau / R, c = cos(th);
    vec3 x1 = x * c + p * (tau * grSinc(th));
    v = normalize(p * c - x * (sin(th) / R));
    x = x1;
}

// Conservative fictitious-time step for clearance d: within the ball of radius d the index
// is at most nMax = n(max(r - d, 0)), and ds = n dtau, so tau = d / nMax travels arc <= d.
float grLuneburgTauStep( in vec3 x, in float d, in float R )
{
    float rm = max(length(x) - d, 0.0) / R;
    return d * inversesqrt(max(2.0 - rm * rm, 1.0));
}

// Fictitious time until the ray leaves the lens (first tau >= 0 with |x| = R). Exact:
// |x(tau)|^2 = R^2 + (r0^2 - R^2) cos(2 tau/R) + R (x0.p0) sin(2 tau/R).
float grLuneburgExitTau( in vec3 x, in vec3 v, in float R )
{
    float r = length(x);
    float n = sqrt(max(2.0 - r * r / (R * R), 0.0));
    return 0.5 * R * atan(max((R - r) * (R + r), 0.0), R * n * dot(x, v) + 1e-30);
}

// Whole traversal of an EMPTY lens, for a ray entering at x (|x| = R) with direction v:
// a quarter period of the oscillator. Exit point R v, exit direction -x/R.
void grLuneburgThrough( inout vec3 x, inout vec3 v, in float R )
{
    vec3 x1 = R * v;
    v = -x / R;
    x = x1;
}

// Optical path length  int n ds = int n^2 dtau  over fictitious time tau from (x, v). Exact.
float grLuneburgOpticalPath( in vec3 x, in vec3 v, in float tau, in float R )
{
    float r = length(x);
    float n = sqrt(max(2.0 - r * r / (R * R), 0.0));
    float th = tau / R, sn = sin(th);
    return tau + (R - r) * (R + r) / R * sn * cos(th) - n * dot(x, v) * sn * sn;
}

// Sphere tracing of GR_SCENE through and inside the lens. Returns the OPTICAL path length
// to the hit (x, v are left at the hit point, v = view tangent there), or -1 on a miss.
float grLuneburgMarch( inout vec3 x, inout vec3 v, in float R, in float lMax, out int steps )
{
    bool inside = dot(x, x) < R * R;
    float opl = 0.0;
    steps = 0;
    for( int i = 0; i < GR_MAX_STEPS; i++ )
    {
        float d = GR_SCENE(x);
        if( d < GR_HIT_EPS ) return opl;
        if( inside )
        {
            float tau = grLuneburgTauStep(x, d, R);
            float tauE = grLuneburgExitTau(x, v, R);
            if( tau >= tauE ) { tau = tauE; inside = false; }
            opl += grLuneburgOpticalPath(x, v, tau, R);
            grLuneburgAdvance(x, v, tau, R);
        }
        else
        {
            float b = dot(x, v), r = length(x);
            float h = b * b - (r - R) * (r + R);
            float tL = (b < 0.0 && h > 0.0) ? max(-b - sqrt(h), 0.0) : 1e30;   // lens entry
            if( d >= tL ) { d = tL; inside = true; }
            x += v * d;
            opl += d;
        }
        steps = i + 1;
        if( opl > lMax ) return -1.0;
    }
    return -1.0;
}

#endif // SDF_GEMS_GR_CURVED_RAYS_GLSL
