// ============================================================================
// pfGaussOccupancy — Gaussian-prefiltered occupancy O_s(p) = (G_s * 1_K)(p):
//   the indicator function of a shape K blurred by an isotropic Gaussian of
//   standard deviation s (per axis). O_s is the quantity that Gaussian-pixel
//   anti-aliasing, drop shadows, blur, and LOD / cone-footprint queries need.
//
//   pfErf, pfErfc, pfPhi            erf / erfc / normal CDF (GLSL has none)
//   pfHalfSpace, pfSlab             half-space and slab                     (closed form)
//   pfBox2D, pfBox3D                axis-aligned box, separable erf product (closed form)
//   pfBall3D                        3D ball: erf + exp                      (closed form)
//   pfDisk2D                        2D disk = 1 - MarcumQ1: no elementary form,
//                                   12-node Gauss-Legendre of a windowed Rice density
//   pfCappedCylinder3D              disk x slab                             (separable)
//   pfRightTri, pfOwenT             Gaussian mass of a right triangle at the origin;
//                                   Owen's T function (6-node Gauss-Legendre + reflection)
//   pfEdge2D, pfTriangle2D,
//   pfPolygon2D                     any simple polygon (convex or not): sum over edges
//   pfSdfOccupancy                  the usual SDF approximation Phi(-d/s) (exact for a half-space only)
//   pfCurvatureCorrected            Phi(-d/s) - s H phi(d/s): first curvature correction (approximate)
//   pfProbit, pfDistance            D = -s * Phi^-1(O): a 1-Lipschitz pseudo-distance to the
//                                   level set O_s = 1/2 of ANY blurred shape (conservative bound)
//
// Sign conventions: SDF d > 0 outside; O_s -> 1 deep inside, -> 0 far outside.
// H = 0.5 * Laplacian(d) at the query point: the MEAN curvature (k1 + k2) / 2 of the SDF
// level set through p; positive for convex bodies (ball: H = 1/|p|; 2D disk: H = 1/(2|p|)).
//
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  Not SDFs. Measured FP32 (GPU) max abs error vs independent FP64 ground truth, for
//         s from 0.01 to 2 x the shape size, evaluated at the FP32-rounded inputs:
//         pfErf 1.3e-7; half-space / slab / box 1.5e-7; pfBall3D 2.3e-7 (closed forms);
//         pfDisk2D 2.4e-7, pfRightTri 3.9e-8, pfOwenT 5.7e-8, pfPolygon2D 7.0e-7,
//         pfTriangle2D 9.9e-7 (fixed-node quadrature). Add ~0.4 ulp(|p|) / s of unavoidable
//         FP32 position rounding (pfCappedCylinder3D: 2.5e-6 at s = 0.01).
//         pfCurvatureCorrected: APPROXIMATE (error O((s k)^2), see docs).
//         pfDistance: CONSERVATIVE 1-Lipschitz bound, exact for a half-space.
//         s is clamped to >= 1e-15. See docs/volume/pfGaussOccupancy.md.
// ============================================================================

#ifndef SDF_GEMS_HAS_PF_GAUSS_OCCUPANCY
#define SDF_GEMS_HAS_PF_GAUSS_OCCUPANCY

// Maximum vertex count of pfPolygon2D (override before including this file).
#ifndef SDF_GEMS_PF_POLY_MAX
#define SDF_GEMS_PF_POLY_MAX 16
#endif

// ---------------------------------------------------------------------------
// 1. erf / erfc / Phi. Two Chebyshev-fitted branches (FP64 fit error 1.6e-9):
//      |x| <  1 : erf(x)  = x P6(x^2)
//      |x| >= 1 : erfc(x) = exp(-x^2) t P7(t),  t = 1 / (1 + |x| / 2)
//    pfErfc keeps RELATIVE accuracy in the upper tail (no 1 - erf cancellation).
// ---------------------------------------------------------------------------
#ifndef SDF_GEMS_HAS_PF_ERF
#define SDF_GEMS_HAS_PF_ERF
float pfErfSmall( in float x )   // |x| <= 1
{
    float u = x * x;
    return x * (1.1283791659 + u * (-0.37612626904 + u * (0.11283597463 + u * (-0.026854328796
             + u * (0.0051893119423 + u * (-0.00080188549892 + u * 0.000078824970970))))));
}

float pfErfcLarge( in float a )  // a >= 1
{
    float t = 1.0 / (1.0 + 0.5 * a);
    float P = 0.28209475535 + t * (0.28210247329 + t * (0.24656491612 + t * (0.17990245839
            + t * (0.059930209400 + t * (0.081285845356 + t * (-0.22310510079 + t * 0.091713959875))))));
    return exp(-a * a) * t * P;
}

float pfErf( in float x )
{
    float a = abs(x);
    if( a < 1.0 ) return pfErfSmall(x);
    float c = pfErfcLarge(a);
    return (x < 0.0) ? (c - 1.0) : (1.0 - c);
}

float pfErfc( in float x )
{
    float a = abs(x);
    if( a < 1.0 ) return 1.0 - pfErfSmall(x);
    float c = pfErfcLarge(a);
    return (x < 0.0) ? (2.0 - c) : c;
}

// Standard normal CDF Phi(x) = P(Z <= x).
float pfPhi( in float x )
{
    return 0.5 * pfErfc(-0.70710678119 * x);
}
#endif // SDF_GEMS_HAS_PF_ERF

// ---------------------------------------------------------------------------
// 2. Half-space, slab, boxes (closed form).
// ---------------------------------------------------------------------------

// Half-space with signed distance d (d > 0 outside): O = Phi(-d / s).
float pfHalfSpace( in float d, in float s )
{
    return pfPhi(-d / max(s, 1e-15));
}

// The usual SDF-based approximation of the blurred occupancy of ANY shape.
// Exact only where the shape looks like a half-space at scale s.
float pfSdfOccupancy( in float d, in float s )
{
    return pfPhi(-d / max(s, 1e-15));
}

// Slab |x| <= b in 1D: O = Phi((b - x)/s) - Phi((-b - x)/s). Written with erfc of
// |x| so that the far tail keeps relative accuracy.
float pfSlab( in float x, in float b, in float s )
{
    float k = 0.70710678119 / max(s, 1e-15);
    float a = abs(x);
    return 0.5 * (pfErfc((a - b) * k) - pfErfc((a + b) * k));
}

// Axis-aligned box with half-extents b, centred at the origin: separable product.
float pfBox2D( in vec2 p, in vec2 b, in float s )
{
    return pfSlab(p.x, b.x, s) * pfSlab(p.y, b.y, s);
}

float pfBox3D( in vec3 p, in vec3 b, in float s )
{
    return pfSlab(p.x, b.x, s) * pfSlab(p.y, b.y, s) * pfSlab(p.z, b.z, s);
}

// ---------------------------------------------------------------------------
// 3. 3D ball of radius R centred at the origin (closed form).
//    With a = |p| / s, b = R / s (the CDF of a non-central chi with 3 dof):
//      O = 1/2 [erf((b - a)/sqrt2) + erf((b + a)/sqrt2)]
//          - [exp(-(b - a)^2 / 2) - exp(-(b + a)^2 / 2)] / (a sqrt(2 pi)).
//    The second term is evaluated as sqrt(2/pi) b exp(-(b-a)^2/2) g(2ab),
//    g(x) = (1 - exp(-x)) / x, which is finite at the centre a = 0.
// ---------------------------------------------------------------------------
float pfBall3D( in float rho, in float R, in float s )
{
    float is = 1.0 / max(s, 1e-15);
    float a = rho * is, b = R * is;
    float dl = (R - rho) * is;              // b - a, subtracted BEFORE scaling (FP32 conditioning)
    float x = 2.0 * a * b;
    float g = (x < 0.25) ? (1.0 + x * (-0.5 + x * (0.16666667 + x * (-0.041666667 + x * (0.0083333333 - x * 0.0013888889)))))
                         : (1.0 - exp(-x)) / x;
    float E = exp(-0.5 * dl * dl);
    float e = 0.5 * (pfErfc(-0.70710678119 * dl) - pfErfc(0.70710678119 * (a + b)));
    return e - 0.79788456080 * b * E * g;
}

float pfBall3D( in vec3 p, in float R, in float s )
{
    return pfBall3D(length(p), R, s);
}

// ---------------------------------------------------------------------------
// 4. 2D disk of radius R (also the infinite 3D cylinder): O = 1 - Q_1(a, b),
//    a = |p| / s, b = R / s, Q_1 the Marcum Q-function. No elementary closed
//    form exists. O is the CDF of the Rice density
//      f(r) = r exp(-(r - a)^2 / 2) I0e(a r),   I0e(x) = exp(-x) I_0(x),
//    which is a smooth bump of unit width around r = a. We integrate f with 12
//    Gauss-Legendre nodes over the part of [a - 6, a + 6] between r = b and the
//    nearer tail (complement when b >= a), so the window is at most 6 wide for
//    every R / s. FP64 algorithm error 2.1e-8.
// ---------------------------------------------------------------------------
float pfBesselI0e( in float x )  // exp(-x) I_0(x), x >= 0; fit rel. error 6e-8
{
    if( x < 4.0 )
    {
        float u = x * x;
        float P = 1.0 + u * (0.25 + u * (1.5624999995e-2 + u * (4.3402778147e-4 + u * (6.7816825320e-6
                + u * (6.7817187845e-8 + u * (4.7090120082e-10 + u * (2.4070834473e-12
                + u * (9.1622413985e-15 + u * 3.5349075117e-17))))))));
        return exp(-x) * P;
    }
    float v = 4.0 / x;
    float P = 0.39894230049 + v * (0.012463766893 + v * (0.0018345820992 + v * (-0.00032521099787
            + v * (0.0037836616459 + v * (-0.0086719348180 + v * (0.011224397168
            + v * (-0.0067080204716 + v * 0.0014602869160)))))));
    return P * inversesqrt(x);
}

const float PF_GL12_X[12] = float[12](
    0.0092196829, 0.0479413718, 0.1150486629, 0.2063410229, 0.3160842505, 0.4373832957,
    0.5626167043, 0.6839157495, 0.7936589771, 0.8849513371, 0.9520586282, 0.9907803171);
const float PF_GL12_W[12] = float[12](
    0.0235876682, 0.0534696630, 0.0800391643, 0.1015837134, 0.1167462683, 0.1245735229,
    0.1245735229, 0.1167462683, 0.1015837134, 0.0800391643, 0.0534696630, 0.0235876682);

float pfDisk2D( in float rho, in float R, in float s )
{
    float is = 1.0 / max(s, 1e-15);
    float a  = rho * is;
    float dl = (R - rho) * is;              // b - a, subtracted BEFORE scaling (FP32 conditioning)
    bool  inside = (dl >= 0.0);
    // integration window in t = r - a
    float lo = inside ? dl : max(-6.0, -a);
    float hi = inside ? max(dl, 6.0) : dl;
    float w  = max(hi - lo, 0.0);
    float I  = 0.0;
    for( int i = 0; i < 12; i++ )
    {
        float t = lo + w * PF_GL12_X[i];
        float r = a + t;
        I += PF_GL12_W[i] * r * exp(-0.5 * t * t) * pfBesselI0e(a * r);
    }
    I *= w;
    return clamp(inside ? (1.0 - I) : I, 0.0, 1.0);
}

float pfDisk2D( in vec2 p, in float R, in float s )
{
    return pfDisk2D(length(p), R, s);
}

// Capped cylinder (axis y, radius R, half-height hh): disk x slab, separable.
float pfCappedCylinder3D( in vec3 p, in float R, in float hh, in float s )
{
    return pfDisk2D(length(p.xz), R, s) * pfSlab(p.y, hh, s);
}

// ---------------------------------------------------------------------------
// 5. Right triangle, Owen's T, polygons. All in units of s (unit Gaussian at 0).
//    pfRightTri(h, t), h, t >= 0: standard bivariate normal mass of the right
//    triangle (0,0), (h,0), (h,t). For t <= h, with a = t / h,
//      A(h, t) = 1/(2 pi) int_0^a (1 - exp(-h^2 (1 + x^2) / 2)) / (1 + x^2) dx
//              = atan(a) / (2 pi) - T(h, a)            (T = Owen's T function),
//    a smooth integrand on [0, a] with a <= 1 (6 Gauss-Legendre nodes, FP64 error
//    1.7e-9). For t > h the two triangles tile the rectangle [0,h] x [0,t]:
//      A(h, t) = erf(h / sqrt2) erf(t / sqrt2) / 4 - A(t, h).
// ---------------------------------------------------------------------------
const float PF_GL6_X[6] = float[6](0.0337652429, 0.1693953068, 0.3806904070, 0.6193095930, 0.8306046932, 0.9662347571);
const float PF_GL6_W[6] = float[6](0.0856622462, 0.1803807865, 0.2339569673, 0.2339569673, 0.1803807865, 0.0856622462);

float pfRightTriCore( in float h, in float a )  // 0 <= a <= 1
{
    float hh = 0.5 * h * h;
    float I = 0.0;
    for( int i = 0; i < 6; i++ )
    {
        float x = a * PF_GL6_X[i];
        float q = 1.0 + x * x;
        I += PF_GL6_W[i] * (1.0 - exp(-hh * q)) / q;
    }
    return 0.15915494309 * a * I;
}

float pfRightTri( in float h, in float t )
{
    h = min(h, 1e15); t = min(t, 1e15);
    if( t <= h ) return pfRightTriCore(h, (h > 0.0) ? t / h : 0.0);
    return 0.25 * pfErf(0.70710678119 * h) * pfErf(0.70710678119 * t) - pfRightTriCore(t, h / t);
}

// Owen's T function T(h, a) = 1/(2 pi) int_0^a exp(-h^2 (1 + x^2) / 2) / (1 + x^2) dx.
// Absolute accuracy ~1e-7 (not relative: T is formed as atan(a)/(2 pi) - A).
float pfOwenT( in float h, in float a )
{
    float ah = abs(h), aa = abs(a);
    float T = 0.15915494309 * atan(aa) - pfRightTri(ah, aa * ah);
    return (a < 0.0) ? -T : T;
}

// Signed unit-Gaussian mass of the triangle (0, v1, v2); positive when counter-clockwise.
// (v1, v2 are already relative to the query point and divided by s.)
float pfEdge2D( in vec2 v1, in vec2 v2 )
{
    vec2  d  = v2 - v1;
    float l2 = dot(d, d);
    if( l2 <= 0.0 ) return 0.0;
    vec2  e  = d * inversesqrt(l2);
    vec2  c  = (dot(v1, v1) < dot(v2, v2)) ? v1 : v2;   // nearer endpoint: shorter lever arm for rounding
    float h  = c.x * e.y - c.y * e.x;     // signed distance from the origin to the edge line
    float t1 = dot(v1, e), t2 = dot(v2, e);
    float ah = abs(h);
    float A1 = pfRightTri(ah, abs(t1)); if( t1 < 0.0 ) A1 = -A1;
    float A2 = pfRightTri(ah, abs(t2)); if( t2 < 0.0 ) A2 = -A2;
    return (h < 0.0) ? (A1 - A2) : (A2 - A1);
}

// Triangle (a, b, c), any orientation.
float pfTriangle2D( in vec2 p, in vec2 a, in vec2 b, in vec2 c, in float s )
{
    float is = 1.0 / max(s, 1e-15);
    a = (a - p) * is; b = (b - p) * is; c = (c - p) * is;
    return clamp(abs(pfEdge2D(a, b) + pfEdge2D(b, c) + pfEdge2D(c, a)), 0.0, 1.0);
}

// Simple polygon v[0..n-1] (convex or concave, either orientation).
float pfPolygon2D( in vec2 p, in vec2 v[SDF_GEMS_PF_POLY_MAX], in int n, in float s )
{
    float is = 1.0 / max(s, 1e-15);
    float O = 0.0;
    vec2  a = (v[n - 1] - p) * is;
    for( int i = 0; i < SDF_GEMS_PF_POLY_MAX; i++ )
    {
        if( i >= n ) break;
        vec2 b = (v[i] - p) * is;
        O += pfEdge2D(a, b);
        a = b;
    }
    return clamp(abs(O), 0.0, 1.0);
}

// ---------------------------------------------------------------------------
// 6. Curvature-corrected SDF occupancy (APPROXIMATE). O_s solves the heat equation
//    dO/dt = Laplacian(O) / 2 with t = s^2. Inserting u = Phi(-d/s) leaves the residual
//    phi(d/s) Laplacian(d) / (2 s); the next term of the expansion cancels it:
//        O ~= Phi(-d/s) - s H phi(d/s),      H = Laplacian(d) / 2  AT THE QUERY POINT
//    (3D: mean curvature (k1 + k2)/2 of the SDF level set through p, e.g.
//    sdhMeanCurvature of glsl/gradients3d/sdhJet3D.glsl; 2D: half the level-set
//    curvature). H > 0 for convex bodies. For a 3D ball this formula is exact up to
//    exp(-(R + |p|)^2 / (2 s^2)). In general the error drops from ~0.4 s H to
//    O((s k)^2); useless at sharp corners. s H is capped at +-1 because H diverges
//    at medial axes (see docs for measured errors and the failure regime s k > 0.3).
// ---------------------------------------------------------------------------
float pfCurvatureCorrected( in float d, in float H, in float s )
{
    s = max(s, 1e-15);
    float x = d / s;
    float c = clamp(s * H, -1.0, 1.0);
    return clamp(pfPhi(-x) - 0.39894228040 * c * exp(-0.5 * x * x), 0.0, 1.0);
}

// ---------------------------------------------------------------------------
// 7. Probit pseudo-distance. For ANY measurable set K, D(p) = -s Phi^-1(O_s(p))
//    is 1-Lipschitz (half-spaces are extremal: |grad O| <= phi(Phi^-1(O)) / s),
//    equals the exact signed distance for a half-space, and its zero set is the
//    level set O_s = 1/2. D - s * Phi^-1(1 - level) bounds the distance to any
//    other level set. FP32: O is clamped to [1e-4, 1 - 1e-4], so |D| <= 3.72 s
//    (the ~3e-7 absolute error of O would otherwise be amplified in the tails).
// ---------------------------------------------------------------------------
float pfProbit( in float u )   // Phi^-1(u); published single-precision erfinv polynomials (see docs)
{
    float m = clamp(min(u, 1.0 - u), 1e-7, 0.5);
    float w = -log(4.0 * m * (1.0 - m));
    float x = 1.0 - 2.0 * m;
    float q;
    if( w < 5.0 )
    {
        w -= 2.5;
        q = 2.81022636e-08;
        q = 3.43273939e-07 + q * w;
        q = -3.5233877e-06 + q * w;
        q = -4.39150654e-06 + q * w;
        q = 0.00021858087 + q * w;
        q = -0.00125372503 + q * w;
        q = -0.00417768164 + q * w;
        q = 0.246640727 + q * w;
        q = 1.50140941 + q * w;
    }
    else
    {
        w = sqrt(w) - 3.0;
        q = -0.000200214257;
        q = 0.000100950558 + q * w;
        q = 0.00134934322 + q * w;
        q = -0.00367342844 + q * w;
        q = 0.00573950773 + q * w;
        q = -0.0076224613 + q * w;
        q = 0.00943887047 + q * w;
        q = 1.00167406 + q * w;
        q = 2.83297682 + q * w;
    }
    float z = 1.41421356237 * q * x;
    return (u < 0.5) ? -z : z;
}

float pfDistance( in float O, in float s )
{
    return -s * pfProbit(clamp(O, 1e-4, 1.0 - 1e-4));
}

#endif // SDF_GEMS_HAS_PF_GAUSS_OCCUPANCY
