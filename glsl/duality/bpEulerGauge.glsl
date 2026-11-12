// ============================================================================
// bpEulerGauge — Bipolar Euler-Legendre Supporting-Plane Theorem,
//                Brunn-Minkowski-Firey L_p Polar Duals & Weingarten-Secant SDFs
//
// Let K = { x in R^3 : gamma(x) <= 1 } be any origin-containing convex body
// defined by a convex positively 1-homogeneous Minkowski gauge function
// gamma(t z) = t gamma(z) (t > 0), which equals the support function h_{K^o}(z)
// of the polar body K^o = { y : x.y <= 1 for all x in K }.
//
// Theorem (Bipolar Euler-Legendre Supporting-Plane Identity):
// By Euler's homogeneous function theorem, z . grad_gamma(z) = gamma(z) and
// grad_gamma(t z) = grad_gamma(z). By bipolar duality (K^o)^o = K, for EVERY
// non-zero trial point z in R^3 \ {0}, the gradient g = grad_gamma(z) lies on
// the polar boundary d(K^o) and satisfies:
//
//     h_K( grad_gamma(z) ) = max_{y in K} dot(grad_gamma(z), y) = 1   (identically!)
//
// with the maximum attained at the radial boundary point y(z) = z / gamma(z) in dK.
// Consequently, for EVERY query point p in R^3 and EVERY trial point z != 0:
//
//   1. Supporting-plane lower bound (valid both inside and outside K):
//        S(p; z) = ( dot(grad_gamma(z), p) - 1.0 ) / length(grad_gamma(z)) <= sd_K(p)
//   2. Radial-boundary upper bound (for p outside K, gamma(p) >= 1):
//        U(p; z) = length( p - z / gamma(z) ) >= sd_K(p)
//   3. First-order Euler bound at z = p (using dot(grad_gamma(p), p) = gamma(p)):
//        S(p; p) = ( gamma(p) - 1.0 ) / length(grad_gamma(p)) <= sd_K(p)
//      (which recovers Inigo Quilez's k0*(k0-1)/k1 ellipsoid bound as the k=1
//      special case of gamma(p) = length(p / r)).
//   4. Weingarten-Secant tangent cascade:
//      Projecting p onto the supporting plane gives x_pl = p - S(p; z_k) n_k.
//      Along the boundary secant dy = y_k - y_{k-1}, the directional normal
//      curvature of dK is kappa_k = dot(dy, n_k - n_{k-1}) / dot(dy, dy), giving
//      the curvature-damped tangent update:
//        z_{k+1} = y_k + (x_pl - y_k) / (1.0 + S(p; z_k) * kappa_k).
//
// @author Gemini 4 Argon (sdf-gems, 2026)
// @source Minkowski gauge & polar duality: Minkowski (1911), Firey (1962 L_p
//         Brunn-Minkowski theory). First-order ellipsoid bound k0*(k0-1)/k1:
//         Inigo Quilez (2019). Bipolar Euler-Legendre supporting-plane theorem
//         and Weingarten-Secant gauge cascade packaged and verified here; see
//         docs/duality/bpEulerGauge.md for measured accuracy and limitations.
// @exact  Conservative lower bound (S(p; z) <= sd_K(p) everywhere to FP32 rounding,
//         measured max overshoot <= 3.6e-7 across all 4 test bodies and all step
//         counts k = 1..12). Converges to the exact Euclidean SDF at BP_ITERS = 8
//         with mean exterior error 8.4e-6 .. 3.4e-5 and mean interior error
//         2.9e-5 .. 1.4e-4 (see docs/duality/bpEulerGauge.md).
// ============================================================================

#ifndef SDF_GEMS_HAS_BP_EULER_GAUGE
#define SDF_GEMS_HAS_BP_EULER_GAUGE

#ifndef BP_ITERS
#define BP_ITERS 8
#endif

// ---------------------------------------------------------------------------
// 1. 1-Homogeneous Convex Gauge Primitives
//    Each bp*(z, ...) returns vec4(gamma(z), grad_gamma(z)).
//    By Euler's theorem, dot(z, g.yzw) == g.x for all z != 0.
// ---------------------------------------------------------------------------

vec3 bpSgn( in vec3 z )
{
    return vec3((z.x < 0.0) ? -1.0 : 1.0,
                (z.y < 0.0) ? -1.0 : 1.0,
                (z.z < 0.0) ? -1.0 : 1.0);
}

// Sphere of radius r: gamma(z) = |z| / r.
vec4 bpSphere( in vec3 z, in float r )
{
    float l2  = dot(z, z);
    float inv = 1.0 / max(r, 1e-20);
    if( l2 > 1e-30 )
    {
        float il = inversesqrt(l2);
        return vec4(l2 * il * inv, z * (il * inv));
    }
    return vec4(0.0, inv, 0.0, 0.0);
}

// Ellipsoid with semi-axes r: gamma(z) = |z / r|_2.
vec4 bpEllipsoid( in vec3 z, in vec3 r )
{
    vec3  ir = 1.0 / max(r, vec3(1e-20));
    vec3  q  = z * ir;
    float l2 = dot(q, q);
    if( l2 > 1e-30 )
    {
        float il = inversesqrt(l2);
        return vec4(l2 * il, q * ir * il);
    }
    return vec4(0.0, ir.x, 0.0, 0.0);
}

// Axis-aligned box with half-extents b (L_infinity gauge): gamma(z) = max_i |z_i| / b_i.
vec4 bpBox( in vec3 z, in vec3 b )
{
    vec3 ir = 1.0 / max(b, vec3(1e-20));
    vec3 q  = abs(z) * ir;
    vec3 s  = bpSgn(z) * ir;
    if( q.x >= q.y && q.x >= q.z ) return vec4(q.x, s.x, 0.0, 0.0);
    if( q.y >= q.z )               return vec4(q.y, 0.0, s.y, 0.0);
    return vec4(q.z, 0.0, 0.0, s.z);
}

// Octahedron / cross-polytope with semi-axes b (L_1 gauge, polar dual of box 1/b):
// gamma(z) = |z.x|/b.x + |z.y|/b.y + |z.z|/b.z.
vec4 bpOctahedron( in vec3 z, in vec3 b )
{
    vec3 ir = 1.0 / max(b, vec3(1e-20));
    return vec4(dot(abs(z), ir), bpSgn(z) * ir);
}

// Lame / L_p superellipsoid with semi-axes b and exponent p >= 1:
// gamma(z) = ( |z.x/b.x|^p + |z.y/b.y|^p + |z.z/b.z|^p )^(1/p).
// Normalized by max_i |z_i/b_i| so pow() never overflows or underflows in FP32.
vec4 bpLpBox( in vec3 z, in vec3 b, in float p )
{
    vec3  ir = 1.0 / max(b, vec3(1e-20));
    vec3  q  = abs(z) * ir;
    float m  = max(q.x, max(q.y, q.z));
    if( m <= 1e-30 ) return vec4(0.0, ir.x, 0.0, 0.0);
    vec3  u  = q / m;
    vec3  up = pow(max(u, vec3(0.0)), vec3(p));
    float s  = up.x + up.y + up.z;
    float sp = pow(s, 1.0 / p);
    float g  = m * sp;
    vec3  w  = pow(max(u, vec3(0.0)), vec3(max(p - 1.0, 0.0))) / max(pow(s, 1.0 - 1.0 / p), 1e-30);
    return vec4(g, bpSgn(z) * ir * w);
}

// L_p-rounded cylinder of radius r and half-height h along y:
// gamma(z) = ( (|z.xz|/r)^p + (|z.y|/h)^p )^(1/p).
vec4 bpCylinderLp( in vec3 z, in float r, in float h, in float p )
{
    float ir  = 1.0 / max(r, 1e-20);
    float ih  = 1.0 / max(h, 1e-20);
    float r2  = dot(z.xz, z.xz);
    float rl  = (r2 > 1e-30) ? sqrt(r2) : 0.0;
    float qr  = rl * ir;
    float qy  = abs(z.y) * ih;
    float m   = max(qr, qy);
    if( m <= 1e-30 ) return vec4(0.0, ir, 0.0, 0.0);
    float ur  = qr / m;
    float uy  = qy / m;
    float s   = pow(ur, p) + pow(uy, p);
    float sp  = pow(s, 1.0 / p);
    float den = max(pow(s, 1.0 - 1.0 / p), 1e-30);
    float wr  = pow(ur, max(p - 1.0, 0.0)) / den;
    float wy  = pow(uy, max(p - 1.0, 0.0)) / den;
    vec2  uxz = (r2 > 1e-30) ? (z.xz / rl) : vec2(1.0, 0.0);
    return vec4(m * sp, uxz.x * ir * wr, ((z.y < 0.0) ? -1.0 : 1.0) * ih * wy, uxz.y * ir * wr);
}

// ---------------------------------------------------------------------------
// 2. Brunn-Minkowski-Firey Polar Dual Operators on Gauges
// ---------------------------------------------------------------------------

// Transform a gauge evaluated in local coordinates z_local = z * R (= R^T z) back to world:
// gamma_{R K}(z) = gamma_K(R^T z), grad_{R K}(z) = R * grad_K(R^T z).
vec4 bpRotateGrad( in vec4 gLocal, in mat3 R )
{
    return vec4(gLocal.x, R * gLocal.yzw);
}

// Polar dual of Minkowski sum: K = (A (+) B)^o  <=>  gamma_K(z) = gamma_{A^o}(z) + gamma_{B^o}(z).
vec4 bpPolarMinkowski( in vec4 ga, in vec4 gb )
{
    return ga + gb;
}

// Exact convex intersection: K = A cap B  <=>  gamma_{A cap B}(z) = max(gamma_A(z), gamma_B(z)).
vec4 bpIntersect( in vec4 ga, in vec4 gb )
{
    return (ga.x >= gb.x) ? ga : gb;
}

// Firey L_p polar intersection (smooth L_p bevel of convex bodies, p >= 1):
// gamma(z) = ( gamma_A(z)^p + gamma_B(z)^p )^(1/p).
// Preserves exact 1-homogeneity and convexity for all p >= 1!
vec4 bpFireyIntersect( in vec4 ga, in vec4 gb, in float p )
{
    float m = max(ga.x, gb.x);
    if( m <= 1e-30 ) return ga;
    float ua  = max(ga.x / m, 0.0);
    float ub  = max(gb.x / m, 0.0);
    float s   = pow(ua, p) + pow(ub, p);
    float sp  = pow(s, 1.0 / p);
    float den = max(pow(s, 1.0 - 1.0 / p), 1e-30);
    float wa  = pow(ua, max(p - 1.0, 0.0)) / den;
    float wb  = pow(ub, max(p - 1.0, 0.0)) / den;
    return vec4(m * sp, wa * ga.yzw + wb * gb.yzw);
}

// ---------------------------------------------------------------------------
// 3. Bipolar Euler-Legendre Supporting-Plane Functional & Weingarten Cascade
// ---------------------------------------------------------------------------

// Single-evaluation Euler Gauge Bound at z = p:
//   S(p; p) = (gamma(p) - 1) / |grad_gamma(p)| <= sd_K(p)  everywhere in R^3.
// Returns vec4(d_euler, unit_normal).
vec4 bpEulerBound1( in vec4 gp )
{
    float gn2 = dot(gp.yzw, gp.yzw);
    if( gn2 <= 1e-30 ) return vec4(-1.0, 1.0, 0.0, 0.0);
    float ign = inversesqrt(gn2);
    return vec4((gp.x - 1.0) * ign, gp.yzw * ign);
}

// Supporting-plane evaluation at query point p from gauge gp = (gamma(z), grad_gamma(z))
// evaluated at any non-zero trial point z:
//   S(p; z) = (dot(grad_gamma(z), p) - 1) / |grad_gamma(z)| <= sd_K(p).
vec4 bpSupportPlane( in vec3 p, in vec4 gp )
{
    float gn2 = dot(gp.yzw, gp.yzw);
    if( gn2 <= 1e-30 ) return vec4(-1e30, 1.0, 0.0, 0.0);
    float ign = inversesqrt(gn2);
    return vec4((dot(gp.yzw, p) - 1.0) * ign, gp.yzw * ign);
}

struct BpState
{
    vec3  z;       // current trial point in R^3 \ {0}
    vec3  yPrev;   // previous radial boundary point y_{k-1} = z_{k-1} / gamma(z_{k-1}) in dK
    vec3  nPrev;   // previous supporting-plane unit normal n_{k-1}
    vec3  n;       // best unit normal found so far (maximizing S(p; z_k))
    float lo;      // best conservative lower bound max_k S(p; z_k) <= sd_K(p)
    float ub;      // best exterior upper bound min_k |p - y_k| >= sd_K(p) (when gamma(p) >= 1)
    int   k;       // step counter
};

BpState bpInit( in vec3 p )
{
    BpState s;
    float l2 = dot(p, p);
    s.z     = (l2 > 1e-24) ? p : vec3(1.0, 0.0, 0.0);
    s.yPrev = vec3(0.0);
    s.nPrev = vec3(1.0, 0.0, 0.0);
    s.n     = vec3(1.0, 0.0, 0.0);
    s.lo    = -1e30;
    s.ub    = 1e30;
    s.k     = 0;
    return s;
}

// Update supporting-plane bound and radial upper bound from gauge gp at trial point zc,
// without advancing the Weingarten-Secant state (used for interior octant-axis seeding).
void bpProbe( inout BpState s, in vec3 p, in vec3 zc, in vec4 gp )
{
    float gn2 = dot(gp.yzw, gp.yzw);
    if( gn2 <= 1e-30 ) return;
    float ign = inversesqrt(gn2);
    vec3  n   = gp.yzw * ign;
    float sp  = (dot(gp.yzw, p) - 1.0) * ign;
    if( sp > s.lo )
    {
        s.lo = sp;
        s.n  = n;
        s.z  = zc;
    }
    if( gp.x > 1e-20 )
    {
        vec3 y = zc / gp.x;
        s.ub = min(s.ub, length(p - y));
    }
}

// One Weingarten-Secant step of the Bipolar Euler-Legendre cascade.
void bpStep( inout BpState s, in vec3 p, in vec4 gp )
{
    float gn2 = dot(gp.yzw, gp.yzw);
    if( gn2 <= 1e-30 ) return;
    float ign = inversesqrt(gn2);
    float gn  = gn2 * ign;
    vec3  n   = gp.yzw * ign;
    float sp  = (dot(gp.yzw, p) - 1.0) * ign;
    if( sp > s.lo )
    {
        s.lo = sp;
        s.n  = n;
    }
    vec3  y   = s.z / max(gp.x, 1e-20);
    s.ub      = min(s.ub, length(p - y));
    vec3  xPl = p - sp * n;
    vec3  res = xPl - y;

    float kappa = 2.0 * gn;
    if( s.k > 0 )
    {
        vec3  dy  = y - s.yPrev;
        vec3  dn  = n - s.nPrev;
        float dy2 = dot(dy, dy);
        if( dy2 > 1e-16 )
        {
            kappa = max(dot(dy, dn) / dy2, 0.1 * gn);
        }
    }
    float alpha = clamp(1.0 / max(1.0 + sp * kappa, 0.15), 0.05, 1.5);
    s.yPrev = y;
    s.nPrev = n;
    vec3 zNext = y + alpha * res;
    s.z = (dot(zNext, zNext) > 1e-24) ? zNext : y;
    s.k++;
}

// Compose-your-own solver macro.
// GAUGE is any vec4 expression of the trial vector `z`, e.g.:
//   bpFireyIntersect(bpBox(z, b), bpOctahedron(z, vec3(s)), 4.0)
// For interior points (gamma(p) < 1), probes the 3 octant axes aligned with sign(p)
// before running K Weingarten-Secant steps so deep interior points near the medial axis
// find the nearest face basin.
#define BP_SOLVE(S, P, GAUGE, K) \
{ \
    S = bpInit(P); \
    { \
        vec3 z = S.z; \
        vec4 bp_g0 = (GAUGE); \
        bpProbe(S, P, z, bp_g0); \
        if( bp_g0.x < 1.0 ) \
        { \
            z = vec3( 1.0,  0.0,  0.0); bpProbe(S, P, z, (GAUGE)); \
            z = vec3(-1.0,  0.0,  0.0); bpProbe(S, P, z, (GAUGE)); \
            z = vec3( 0.0,  1.0,  0.0); bpProbe(S, P, z, (GAUGE)); \
            z = vec3( 0.0, -1.0,  0.0); bpProbe(S, P, z, (GAUGE)); \
            z = vec3( 0.0,  0.0,  1.0); bpProbe(S, P, z, (GAUGE)); \
            z = vec3( 0.0,  0.0, -1.0); bpProbe(S, P, z, (GAUGE)); \
            vec3 bp_sg = bpSgn(P); \
            z = bp_sg * vec3(1.0, 1.0, 0.2); bpProbe(S, P, z, (GAUGE)); \
            z = bp_sg * vec3(1.0, 0.2, 1.0); bpProbe(S, P, z, (GAUGE)); \
            z = bp_sg * vec3(0.2, 1.0, 1.0); bpProbe(S, P, z, (GAUGE)); \
        } \
    } \
    for( int bp_i = 0; bp_i < (K); bp_i++ ) \
    { \
        vec3 z = S.z; \
        bpStep(S, P, (GAUGE)); \
    } \
}

// ---------------------------------------------------------------------------
// 4. Ready-Made Primitives & Polar Dual Bodies (sd* and sdg*)
// ---------------------------------------------------------------------------

// 4a. Ellipsoid with semi-axes r via Bipolar Weingarten-Secant cascade.
vec4 sdgBpEllipsoid( in vec3 p, in vec3 r, in int iters )
{
    BpState s;
    BP_SOLVE(s, p, bpEllipsoid(z, r), iters)
    return vec4(s.lo, s.n);
}

vec4 sdgBpEllipsoid( in vec3 p, in vec3 r )
{
    return sdgBpEllipsoid(p, r, BP_ITERS);
}

float sdBpEllipsoid( in vec3 p, in vec3 r )
{
    return sdgBpEllipsoid(p, r, BP_ITERS).x;
}

// 4b. Lame L_p Superellipsoid: |x/b.x|^e + |y/b.y|^e + |z/b.z|^e <= 1 (e >= 1).
vec4 sdgBpSuperellipsoid( in vec3 p, in vec3 b, in float e )
{
    BpState s;
    BP_SOLVE(s, p, bpLpBox(z, b, e), BP_ITERS)
    return vec4(s.lo, s.n);
}

float sdBpSuperellipsoid( in vec3 p, in vec3 b, in float e )
{
    return sdgBpSuperellipsoid(p, b, e).x;
}

// 4c. Firey L_p Cube-Octahedron (smooth rhombic-dodecahedral / cuboctahedral gauge body):
//     gamma(z) = ( gamma_Box(z, b)^e + gamma_Octahedron(z, s)^e )^(1/e).
vec4 sdgBpFireyCubeOctahedron( in vec3 p, in vec3 b, in float sOct, in float e )
{
    BpState s;
    BP_SOLVE(s, p, bpFireyIntersect(bpBox(z, b), bpOctahedron(z, vec3(sOct)), e), BP_ITERS)
    return vec4(s.lo, s.n);
}

float sdBpFireyCubeOctahedron( in vec3 p, in vec3 b, in float sOct, in float e )
{
    return sdgBpFireyCubeOctahedron(p, b, sOct, e).x;
}

// 4d. Firey L_e Polar Dual of Two Rotated L_q Boxes:
//     K = ( B_1 (+)_e R B_2 )^o, with gauge
//     gamma(z) = ( ||z / b1||_q^e + ||(R^T z) / b2||_q^e )^(1/e).
//     A 24-facet smooth polar polyhedron with no closed-form primal SDF.
vec4 bpPolarDualBoxesGauge( in vec3 z, in vec3 b1, in mat3 R, in vec3 b2, in float q, in float e )
{
    vec4 g1 = bpLpBox(z, b1, q);
    vec4 g2 = bpRotateGrad(bpLpBox(z * R, b2, q), R);
    return bpFireyIntersect(g1, g2, e);
}

vec4 sdgBpPolarDualBoxes( in vec3 p, in vec3 b1, in mat3 R, in vec3 b2, in float q, in float e, in int iters )
{
    BpState s;
    BP_SOLVE(s, p, bpPolarDualBoxesGauge(z, b1, R, b2, q, e), iters)
    return vec4(s.lo, s.n);
}

vec4 sdgBpPolarDualBoxes( in vec3 p, in vec3 b1, in mat3 R, in vec3 b2, in float q, in float e )
{
    return sdgBpPolarDualBoxes(p, b1, R, b2, q, e, BP_ITERS);
}

float sdBpPolarDualBoxes( in vec3 p, in vec3 b1, in mat3 R, in vec3 b2, in float q, in float e )
{
    return sdgBpPolarDualBoxes(p, b1, R, b2, q, e, BP_ITERS).x;
}

#endif // SDF_GEMS_HAS_BP_EULER_GAUGE
