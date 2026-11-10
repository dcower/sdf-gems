// ============================================================================
// tpGaussBonnet — integrands for topological invariants of SDF scenes:
//   Gaussian / mean curvature of level sets (Gauss-Bonnet -> Euler characteristic, genus),
//   the polygonal Gauss linking number (closed form per segment pair), and solid angles.
//
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  The curvature formulas and the segment-pair linking term are exact closed forms.
//         The INTEGRALS (Euler characteristic, total mean curvature) are NOT computed here:
//         they are Monte Carlo / quadrature sums of these integrands, done by the caller
//         (see docs/topology/tpGaussBonnet.md for the estimators and their measured error bars).
//         tpGenusChain is an implicit surface (log-sum-exp union), 1-Lipschitz, NOT an exact SDF.
// @source Goldman, "Curvature formulas for implicit curves and surfaces", CAGD 22 (2005);
//         van Oosterom & Strackee, IEEE Trans. Biomed. Eng. 30 (1983) (triangle solid angle);
//         Klenin & Langowski, Biopolymers 54 (2000) (polygonal linking number).
// ============================================================================

#ifndef SDF_GEMS_TP_GAUSS_BONNET_GLSL
#define SDF_GEMS_TP_GAUSS_BONNET_GLSL

// ---------------------------------------------------------------------------
// 1. Curvature of the level set through p from the gradient g and Hessian H
// ---------------------------------------------------------------------------

// Gaussian curvature of the level set of ANY C^2 implicit function f:
//   K = g^T adj(H) g / |g|^4          (Goldman 2005, eq. 4.1)
// H must be symmetric. Independent of the sign and scale of f.
float tpGaussCurvatureImplicit( vec3 g, mat3 H )
{
    vec3 c = vec3(dot(g, cross(H[1], H[2])), dot(g, cross(H[2], H[0])), dot(g, cross(H[0], H[1])));
    float g2 = dot(g, g);
    return dot(g, c) / max(g2 * g2, 1e-30);
}

// Mean curvature (k1 + k2) / 2 of the level set of any C^2 implicit function f, positive where the
// level set is convex as seen from the side the gradient points to (sphere of radius r: +1 / r):
//   Hm = (|g|^2 tr H - g^T H g) / (2 |g|^3)
float tpMeanCurvatureImplicit( vec3 g, mat3 H )
{
    float g2 = dot(g, g);
    float tr = H[0][0] + H[1][1] + H[2][2];
    return (g2 * tr - dot(g, H * g)) / max(2.0 * g2 * sqrt(g2), 1e-30);
}

// Exact SDF only (|grad d| = 1  =>  H g = 0, H has eigenvalues {0, k1, k2}):
// K = k1 k2 = sum of the three principal 2x2 minors of H. No gradient needed.
float tpGaussCurvatureSDF( mat3 H )
{
    return H[0][0] * H[1][1] - H[0][1] * H[0][1]
         + H[1][1] * H[2][2] - H[1][2] * H[1][2]
         + H[2][2] * H[0][0] - H[2][0] * H[2][0];
}

// Exact SDF only: (k1 + k2) / 2 = tr(H) / 2 = laplacian(d) / 2.
float tpMeanCurvatureSDF( mat3 H )
{
    return 0.5 * (H[0][0] + H[1][1] + H[2][2]);
}

// Gradient and Hessian from a 19-point central stencil of step h (for SDFs without an analytic jet).
//   c            = f(p)
//   ap, am       = f(p + h e_i), f(p - h e_i)                       for i = x, y, z
//   dpp, dpm,
//   dmp, dmm     = f(p + s1 h e_i + s2 h e_j) for (i, j) = (x, y), (y, z), (z, x) and signs (s1, s2)
// Measured in FP32 for |p| ~ 1 (see doc): rounding noise ~ ulp(f) / h^2 biases the integral for h <= 3e-4, truncation
// biases it for h >~ (smallest curvature radius) / 20. Usable window: 1e-3 <= h <= 3e-3 for features >= 0.05.
// Feed the result to tpGaussCurvatureImplicit, NOT tpGaussCurvatureSDF (the latter needs H g = 0, which a
// finite-difference Hessian only satisfies to O(h^2): measured bias 10x - 100x larger).
void tpJetFromStencil( float h, float c, vec3 ap, vec3 am, vec3 dpp, vec3 dpm, vec3 dmp, vec3 dmm,
                       out vec3 g, out mat3 H )
{
    g = (ap - am) / (2.0 * h);
    vec3 dd = (ap - 2.0 * c + am) / (h * h);
    vec3 od = (dpp - dpm - dmp + dmm) / (4.0 * h * h);   // (xy, yz, zx)
    H = mat3(dd.x, od.x, od.z,
             od.x, dd.y, od.y,
             od.z, od.y, dd.z);
}

// ---------------------------------------------------------------------------
// 2. Concrete example: analytic 2-jet of a torus and of a genus-n chain of tori
// ---------------------------------------------------------------------------

// Exact SDF of the ring torus (axis y, t = (R, r), R > r) with analytic gradient and Hessian.
// H = (1 / u) m m^T + ((rho - R) / (u rho)) e_phi e_phi^T,  m = meridian tangent, u = distance to the centre circle.
float tpTorusJet( vec3 p, vec2 t, out vec3 g, out mat3 H )
{
    float rho = max(length(p.xz), 1e-12);
    vec3 er = vec3(p.x, 0.0, p.z) / rho;
    vec3 ep = vec3(-p.z, 0.0, p.x) / rho;
    vec2 q = vec2(rho - t.x, p.y);
    float u = max(length(q), 1e-12);
    g = (q.x * er + vec3(0.0, q.y, 0.0)) / u;
    vec3 m = (-q.y * er + vec3(0.0, q.x, 0.0)) / u;
    H = outerProduct(m, m) / u + outerProduct(ep, ep) * (q.x / (u * rho));
    return u - t.y;
}

// Log-sum-exp (C-infinity) union of n <= 8 tori t = (R, r) lying in the xz-plane with centres spaced
// 2 c apart along x. For r < c < R + r and small k the surface is closed with genus n
// (Euler characteristic 2 - 2 n). Returns f (1-Lipschitz, f <= min d_i, NOT an exact SDF) and
// its exact analytic gradient g and Hessian H:
//   f = -k log sum_i exp(-d_i / k),  g = sum w_i g_i,
//   H = sum w_i H_i - (1 / k) (sum w_i g_i g_i^T - g g^T),   w = softmax(-d / k).
float tpGenusChain( vec3 p, int n, vec2 t, float c, float k, out vec3 g, out mat3 H )
{
    float d[8]; vec3 gi[8]; mat3 Hi[8];
    float m = 1e30;
    for (int i = 0; i < 8; ++i)
    {
        if (i >= n) break;
        float xc = c * float(2 * i - (n - 1));
        d[i] = tpTorusJet(p - vec3(xc, 0.0, 0.0), t, gi[i], Hi[i]);
        m = min(m, d[i]);
    }
    float sw = 0.0;
    g = vec3(0.0);
    H = mat3(0.0);
    mat3 G2 = mat3(0.0);
    for (int i = 0; i < 8; ++i)
    {
        if (i >= n) break;
        float w = exp(-(d[i] - m) / k);
        sw += w;
        g += w * gi[i];
        H += w * Hi[i];
        G2 += w * outerProduct(gi[i], gi[i]);
    }
    g /= sw;
    H = H / sw - (G2 / sw - outerProduct(g, g)) / k;
    return m - k * log(sw);
}

// Distance-only convenience wrapper (for ray marching the example surface).
float tpGenusChain( vec3 p, int n, vec2 t, float c, float k )
{
    vec3 g; mat3 H;
    return tpGenusChain(p, n, t, c, k, g, H);
}

// Gauss-Bonnet integrand of the example: vec4(f, K |grad f|, Hm |grad f|, |grad f|).
// Coarea formula: for uniform samples p of a box of volume V,
//   V / (b - a) * mean( [a < f(p) < b] * .yzw )  ->  (2 pi chi, total mean curvature, area),
// each averaged over the level sets { f = t }, a < t < b. chi is an integer only if no critical value of f lies in (a, b).
vec4 tpGenusChainIntegrand( vec3 p, int n, vec2 t, float c, float k )
{
    vec3 g; mat3 H;
    float f = tpGenusChain(p, n, t, c, k, g, H);
    float gl = length(g);
    return vec4(f, tpGaussCurvatureImplicit(g, H) * gl, tpMeanCurvatureImplicit(g, H) * gl, gl);
}

// ---------------------------------------------------------------------------
// 3. Solid angles and the polygonal Gauss linking number
// ---------------------------------------------------------------------------

// Signed solid angle of the triangle (a, b, c) seen from the origin (van Oosterom & Strackee 1983).
// Sign = sign(det(a, b, c)): positive when a, b, c appear counter-clockwise to a viewer on the far side of the
// triangle looking back at the origin.
// Range [-2 pi, 2 pi]; 1 atan, 3 sqrt. If the origin lies IN the triangle the value (+-2 pi) is ambiguous.
float tpSolidAngleTri( vec3 a, vec3 b, vec3 c )
{
    float la = length(a), lb = length(b), lc = length(c);
    float num = dot(a, cross(b, c));
    float den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
    // atan(0, 0) is undefined in GLSL (NaN on some GPUs): degenerate triangles with a vertex at the origin return 0.
    return (abs(num) + abs(den) > 0.0) ? 2.0 * atan(num, den) : 0.0;
}

// Exact Gauss linking double integral of two straight segments a0->a1 and b0->b1:
//   (1 / 4 pi) int int (a - b) . (da x db) / |a - b|^3
// = (signed solid angle of the parallelogram { b(t) - a(s) } seen from the origin) / (4 pi).
// Summed over all segment pairs of two disjoint closed polygons this is EXACTLY their integer linking
// number (not an approximation in the number of segments). Each term lies in (-1/2, 1/2).
// Undefined if the segments intersect. For the writhe of one polygon, sum over i != j and SKIP
// adjacent segments (they share a vertex; their exact contribution is 0).
float tpLinkingSegSeg( vec3 a0, vec3 a1, vec3 b0, vec3 b1 )
{
    vec3 v00 = b0 - a0, v10 = b0 - a1, v11 = b1 - a1, v01 = b1 - a0;
    return (tpSolidAngleTri(v00, v10, v11) + tpSolidAngleTri(v00, v11, v01)) * 0.07957747154594767;
}

// Example curves: the two antipodal torus strands used by sdCinquefoilKnot / sdSolomonLink,
//   q_s(t) = ((R + s r0 cos(k t / 2)) cos t, (R + s r0 cos(k t / 2)) sin t, s r0 sin(k t / 2)),  s = +-1.
// For even k = 2 n they are two closed curves forming the torus link T(2, 2 n).
vec3 tpTorusStrand( float tt, float R, float r0, int k, float s )
{
    float th = 0.5 * float(k) * tt;
    float w = R + s * r0 * cos(th);
    return vec3(w * cos(tt), w * sin(tt), s * r0 * sin(th));
}

#ifndef TP_LINK_SEGMENTS
#define TP_LINK_SEGMENTS 24
#endif

// Fixed-N polygonal linking number of the two strands of the torus link T(2, 2 n) (k = 2 n lobes):
// TP_LINK_SEGMENTS^2 calls of tpLinkingSegSeg. Returns Lk = -n (both strands traversed with increasing t) up to
// FP32 rounding (measured <= 1e-5), provided TP_LINK_SEGMENTS is large enough that the inscribed polygons are
// isotopic to the curves (measured for R = 0.5, r0 = 0.2: N >= 2 n + 1 in FP64 for n <= 4; confirmed on the
// GPU at N = 6, 6, 8, 12 for n = 1, 2, 3, 4). Too few segments give a WRONG INTEGER, not a fractional value.
float tpLinkingTorusLink( float R, float r0, int n )
{
    const float TAU = 6.283185307179586;
    float lk = 0.0;
    vec3 a0 = tpTorusStrand(0.0, R, r0, 2 * n, 1.0);
    for (int i = 0; i < TP_LINK_SEGMENTS; ++i)
    {
        vec3 a1 = tpTorusStrand(TAU * float(i + 1) / float(TP_LINK_SEGMENTS), R, r0, 2 * n, 1.0);
        vec3 b0 = tpTorusStrand(0.0, R, r0, 2 * n, -1.0);
        for (int j = 0; j < TP_LINK_SEGMENTS; ++j)
        {
            vec3 b1 = tpTorusStrand(TAU * float(j + 1) / float(TP_LINK_SEGMENTS), R, r0, 2 * n, -1.0);
            lk += tpLinkingSegSeg(a0, a1, b0, b1);
            b0 = b1;
        }
        a0 = a1;
    }
    return lk;
}

// Solid angle subtended at p by the regular TP_LINK_SEGMENTS-gon inscribed in the circle of radius R
// in the xy-plane (fan of triangles from the centre; counter-clockwise seen from +z).
// Omega / (4 pi) is the generalized winding number of the flat spanning disc: it jumps by exactly 1
// across the disc and is harmonic elsewhere. Omega = -2 pi (sign z)(1 - |z| / sqrt(z^2 + R^2)) on the axis
// in the limit of many segments.
float tpSolidAngleCircle( vec3 p, float R )
{
    const float TAU = 6.283185307179586;
    float om = 0.0;
    vec3 a = vec3(R, 0.0, 0.0) - p;
    for (int i = 0; i < TP_LINK_SEGMENTS; ++i)
    {
        float t1 = TAU * float(i + 1) / float(TP_LINK_SEGMENTS);
        vec3 b = vec3(R * cos(t1), R * sin(t1), 0.0) - p;
        om += tpSolidAngleTri(-p, a, b);
        a = b;
    }
    return om;
}

#endif // SDF_GEMS_TP_GAUSS_BONNET_GLSL
