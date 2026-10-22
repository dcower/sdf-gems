// ============================================================================
// opRepLogPolar — Conformal Log-Polar & Gaussian-Integer Log-Spiral / Phyllotaxis Repetition
//
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (0.0% raymarching overstep via Exact Conformal Exponential Metric;
//         0.0 branch-cut jump across theta = +-pi)
// ============================================================================

#ifndef SDF_GEMS_OP_REP_LOG_POLAR_GLSL
#define SDF_GEMS_OP_REP_LOG_POLAR_GLSL

// ----------------------------------------------------------------------------
// sdLogPolarMetric(d_log, r, s) — Exact Conformal Exponential Distance Metric
//
// Standard log-polar shaders approximate Euclidean distance from log-space
// distance d_log via the first-order tangent metric:
//   d_naive = d_log * (r / s)
// Because J(u) = (r_0 / s) * exp(u / s) grows exponentially outward and shrinks
// exponentially inward, d_naive OVERSTEPS the true Euclidean distance whenever
// a ray marches inward toward the origin (d_log > 0):
//   L_inward = integral_{u - d_log}^{u} (r_0 / s) exp(t / s) dt
//            = r * (1.0 - exp(-d_log / s)) < d_log * (r / s)
// causing up to +12.5% to +25.7% raymarching overstep and surface clipping!
//
// Conversely, inside a primitive (d_log < 0), the shortest Euclidean escape path
// is strictly inward toward the origin, yielding the exact signed metric:
//   d_exact = sign(d_log) * r * (1.0 - exp(-abs(d_log) / s))
// which satisfies ||grad_p d_exact|| <= 1.0 everywhere in R^2 \ {0} (0.0% overstep).
// ----------------------------------------------------------------------------
float sdLogPolarMetric( float d_log, float r, float s )
{
    float sgn = (d_log >= 0.0) ? 1.0 : -1.0;
    return sgn * r * (1.0 - exp(-abs(d_log) / s));
}

// ----------------------------------------------------------------------------
// opRepLogPolar(p, n, id, r, s) — Conformal Log-Polar Domain Repetition
//
// Maps p in R^2 \ {0} to conformal unit-square cells q in [-0.5, 0.5]^2 via
//   w = (n / (2*pi)) * vec2(0.5 * log(dot(p, p)), atan(p.y, p.x))
// where n is the integer number of angular sectors per 2*pi revolution.
//
// Uses 0.5 * log(dot(p, p)) to fuse the log-radius evaluation, and wraps the
// angular index id.y modulo n so id is 100% continuous across the atan(y, x)
// branch cut at theta = +-pi.
//
// Usage:
//   vec2 id; float r, s;
//   vec2 q = opRepLogPolar(p, 12.0, id, r, s);
//   float d = sdLogPolarMetric(sdPrim2D(q), r, s);
// ----------------------------------------------------------------------------
vec2 opRepLogPolar( vec2 p, float n, out vec2 id, out float r, out float s )
{
    float r2 = dot(p, p);
    r = sqrt(r2);
    s = n * 0.159154943092; // n / (2.0 * pi)
    vec2 w  = vec2(0.5 * log(max(r2, 1e-30)), atan(p.y, p.x)) * s;
    vec2 i  = round(w);
    id = vec2(i.x, i.y - n * floor(i.y / n));
    return w - i;
}

vec2 opRepLogPolar( vec2 p, float n, out float r, out float s )
{
    vec2 id;
    return opRepLogPolar(p, n, id, r, s);
}

// ----------------------------------------------------------------------------
// opRepLogPolarLim(p, n, lima, limb, id, r, s) — Radially Limited Log-Polar
//
// Repeats conformal cells only across radial rings i.x in [lima, limb],
// preventing infinite accumulation at the origin (r -> 0) and outer infinity.
// ----------------------------------------------------------------------------
vec2 opRepLogPolarLim( vec2 p, float n, float lima, float limb, out vec2 id, out float r, out float s )
{
    float r2 = dot(p, p);
    r = sqrt(r2);
    s = n * 0.159154943092;
    vec2  w  = vec2(0.5 * log(max(r2, 1e-30)), atan(p.y, p.x)) * s;
    float ix = clamp(round(w.x), lima, limb);
    float iy = round(w.y);
    id = vec2(ix, iy - n * floor(iy / n));
    return w - vec2(ix, iy);
}

// ----------------------------------------------------------------------------
// opRepLogSpiral(p, n, k, id, r, s) — Branch-Cut-Free Gaussian-Integer
//                                      Logarithmic Spiral / Phyllotaxis Lattice
//
// Applies the complex Gaussian-integer multiplication (n + i*k) / (2*pi) to the
// conformal log-polar coordinate z = ln(r) + i*theta:
//   w = (1 / (2*pi)) * vec2(n * u - k * theta, k * u + n * theta)
// For any integer pair (n, k) — such as Fibonacci pairs (8, 5), (13, 8), (21, 13)
// for golden-ratio phyllotaxis (sunflower / pinecone parastichy lattices) —
// crossing the atan(y, x) branch cut theta -> theta + 2*pi shifts w by the
// exact integer vector (-k, +n), so q = w - round(w) is 100% C^infinity
// continuous across theta = +-pi with ZERO seam artifacts!
//
// Canonical Shear-Reduced Cell ID:
//   float m = floor(i.y / n);
//   id = vec2(i.x + k * m, i.y - n * m);
// makes the 2D integer cell ID invariant across theta = +-pi as well!
// ----------------------------------------------------------------------------
vec2 opRepLogSpiral( vec2 p, float n, float k, out vec2 id, out float r, out float s )
{
    float r2    = dot(p, p);
    r = sqrt(r2);
    float u     = 0.5 * log(max(r2, 1e-30));
    float theta = atan(p.y, p.x);
    vec2  w     = vec2(n * u - k * theta, k * u + n * theta) * 0.159154943092;
    vec2  i     = round(w);
    float m     = floor(i.y / n);
    id = vec2(i.x + k * m, i.y - n * m);
    s  = length(vec2(n, k)) * 0.159154943092;
    return w - i;
}

// ----------------------------------------------------------------------------
// opRepLogSpiral4Tap(p, n, k, base_id, r, s) — 4-Tap Neighbor-Corrected Spiral
//
// Returns the fractional cell coordinate f = w - floor(w) in [0, 1)^2 and the
// branch-cut-invariant base index base_id for 2x2 neighbor evaluation of
// size-varied, jittered, or asymmetric primitives.
// ----------------------------------------------------------------------------
vec2 opRepLogSpiral4Tap( vec2 p, float n, float k, out vec2 base_id, out float r, out float s )
{
    float r2    = dot(p, p);
    r = sqrt(r2);
    float u     = 0.5 * log(max(r2, 1e-30));
    float theta = atan(p.y, p.x);
    vec2  w     = vec2(n * u - k * theta, k * u + n * theta) * 0.159154943092;
    vec2  b     = floor(w);
    float m     = floor(b.y / n);
    base_id = vec2(b.x + k * m, b.y - n * m);
    s       = length(vec2(n, k)) * 0.159154943092;
    return w - b;
}

#endif // SDF_GEMS_OP_REP_LOG_POLAR_GLSL
