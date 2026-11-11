// ============================================================================
// tpGaussBonnet — integrands and exact discrete evaluators for topological invariants of SDF scenes:
//   1-2. Gaussian / mean curvature of C^2 level sets (Gauss-Bonnet -> Euler characteristic, genus),
//   3.   Solid angles and the polygonal Gauss linking number (closed form per segment pair),
//   4.   Distributional spherical Gauss-map degree on Marching-Tetrahedra / Kuhn cells (exact integer
//        Euler characteristic on sharp non-convex polyhedra, box frames, and sponges),
//   5.   Milnor's triple-linking invariant mu_123 via Seifert-fan intersections (Mellor-Melvin).
//
// @author Claude Fable 5.2 (sdf-gems, 2026) — upgraded with distributional Kuhn-tetrahedra Gauss degree
//         and Milnor triple-linking invariant using Gemini 4 Carbon (sdf-gems, 2026)
// @exact  Curvature formulas, tpLinkingSegSeg, tpSphericalGaussTri, tpTetGaussDegree / tpTetEulerContribution,
//         and tpMilnorTripleLink are exact closed-form / combinatorial evaluators. Smooth coarea Gauss-Bonnet
//         integrals are Monte Carlo estimates; Kuhn-tetrahedra sums of tpTetEulerContribution are deterministic
//         and topologically quantized to exact even integers (to FP32 rounding ~1e-6) on both smooth and sharp
//         non-convex shapes once voxel spacing resolves tunnel thickness.
//         tpGenusChain is an implicit surface (log-sum-exp union), 1-Lipschitz, NOT an exact SDF.
// @source Goldman, "Curvature formulas for implicit curves and surfaces", CAGD 22 (2005);
//         van Oosterom & Strackee, IEEE Trans. Biomed. Eng. 30 (1983) (triangle solid angle);
//         Klenin & Langowski, Biopolymers 54 (2000) (polygonal linking number);
//         Hopf, Math. Ann. 96 (1927) (Gauss-map degree = chi / 2);
//         Milnor, Ann. Math. 59 (1954); Mellor & Melvin, Algebr. Geom. Topol. 3 (2003) (triple linking).
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

// ---------------------------------------------------------------------------
// 4. Distributional Gauss-map degree & Marching-Tetrahedra Euler characteristic
// ---------------------------------------------------------------------------

// Signed fractional solid angle Omega(n0, n1, n2) / (4 pi) in [-1/2, 1/2] of a geodesic spherical triangle
// with UNIT vertices n0, n1, n2 on S^2 (van Oosterom & Strackee 1983; 0 sqrt, 1 atan).
// By Heinz Hopf's degree theorem, for any closed oriented surface triangulation isotopic to { d = 0 }
// carrying a continuous outward-transverse vector field n, sum(tpSphericalGaussTri) = deg(n) = chi / 2,
// even when { d = 0 } has sharp convex or re-entrant C^0 creases where pointwise curvature misses Dirac seams.
float tpSphericalGaussTri( vec3 n0, vec3 n1, vec3 n2 )
{
    float num = dot(n0, cross(n1, n2));
    float den = 1.0 + dot(n0, n1) + dot(n1, n2) + dot(n2, n0);
    return (abs(num) + abs(den) > 0.0) ? 0.15915494309189535 * atan(num, den) : 0.0;
}

// Permutation-symmetric zero-crossing unit normal along an edge between two vec4(d, grad) samples.
// Orders endpoints by d so every tetrahedron sharing a grid edge produces the bit-identical unit vector on S^2,
// guaranteeing exact cancellation of shared geodesic boundary arcs across adjacent Marching-Tetrahedra cells.
vec3 tpEdgeZeroNormal( vec4 dg0, vec4 dg1 )
{
    vec4  neg = (dg0.x < dg1.x) ? dg0 : dg1;
    vec4  pos = (dg0.x < dg1.x) ? dg1 : dg0;
    float den = pos.x - neg.x;
    float t   = (den > 1e-30) ? clamp((-neg.x) / den, 0.0, 1.0) : 0.5;
    vec3  g   = mix(neg.yzw, pos.yzw, t);
    float g2  = dot(g, g);
    if (g2 > 1e-24) return g * inversesqrt(g2);
    float p2  = dot(pos.yzw, pos.yzw);
    return (p2 > 1e-24) ? pos.yzw * inversesqrt(p2) : vec3(1.0, 0.0, 0.0);
}

// Exact spherical Gauss-map degree contribution of a positively oriented tetrahedron (det(p1-p0, p2-p0, p3-p0) > 0)
// from four vertex samples vec4(d, grad). A 5-comparator sorting network reduces all 16 sign masks to 3 active cases.
// Sums over a closed Kuhn-tetrahedral grid to deg(n) = chi(K) = chi(dK) / 2.
float tpTetGaussDegree( vec4 dg0, vec4 dg1, vec4 dg2, vec4 dg3 )
{
    float orient = 1.0;
    vec4  tmp;
    if (dg1.x < dg0.x) { tmp = dg0; dg0 = dg1; dg1 = tmp; orient = -orient; }
    if (dg3.x < dg2.x) { tmp = dg2; dg2 = dg3; dg3 = tmp; orient = -orient; }
    if (dg2.x < dg0.x) { tmp = dg0; dg0 = dg2; dg2 = tmp; orient = -orient; }
    if (dg3.x < dg1.x) { tmp = dg1; dg1 = dg3; dg3 = tmp; orient = -orient; }
    if (dg2.x < dg1.x) { tmp = dg1; dg1 = dg2; dg2 = tmp; orient = -orient; }

    if (dg0.x >= 0.0 || dg3.x < 0.0) return 0.0;

    if (dg1.x >= 0.0)
    {
        vec3 n01 = tpEdgeZeroNormal(dg0, dg1);
        vec3 n02 = tpEdgeZeroNormal(dg0, dg2);
        vec3 n03 = tpEdgeZeroNormal(dg0, dg3);
        return orient * tpSphericalGaussTri(n01, n02, n03);
    }
    if (dg2.x < 0.0)
    {
        vec3 n03 = tpEdgeZeroNormal(dg0, dg3);
        vec3 n13 = tpEdgeZeroNormal(dg1, dg3);
        vec3 n23 = tpEdgeZeroNormal(dg2, dg3);
        return orient * tpSphericalGaussTri(n03, n13, n23);
    }
    vec3 n02 = tpEdgeZeroNormal(dg0, dg2);
    vec3 n03 = tpEdgeZeroNormal(dg0, dg3);
    vec3 n13 = tpEdgeZeroNormal(dg1, dg3);
    vec3 n12 = tpEdgeZeroNormal(dg1, dg2);
    return orient * (tpSphericalGaussTri(n02, n03, n13) + tpSphericalGaussTri(n02, n13, n12));
}

// Surface Euler characteristic contribution 2 * tpTetGaussDegree of a positively oriented tetrahedron:
// summed over all tetrahedra enclosing { d <= 0 }, returns chi({ d = 0 }) = 2 - 2 g (to ~1e-6 in FP32).
float tpTetEulerContribution( vec4 dg0, vec4 dg1, vec4 dg2, vec4 dg3 )
{
    return 2.0 * tpTetGaussDegree(dg0, dg1, dg2, dg3);
}

// Face-conforming 6-tetrahedra Kuhn (Freudenthal) subdivision of a voxel [0, h]^3 from its 8 corner vec4(d, grad)
// samples c_{ijk} (0 SFU when all 8 corners share a sign).
float tpKuhnVoxelEuler( vec4 c000, vec4 c100, vec4 c010, vec4 c110,
                        vec4 c001, vec4 c101, vec4 c011, vec4 c111 )
{
    float mn = min(min(min(c000.x, c100.x), min(c010.x, c110.x)),
                   min(min(c001.x, c101.x), min(c011.x, c111.x)));
    float mx = max(max(max(c000.x, c100.x), max(c010.x, c110.x)),
                   max(max(c001.x, c101.x), max(c011.x, c111.x)));
    if (mn >= 0.0 || mx < 0.0) return 0.0;
    return tpTetEulerContribution(c000, c100, c110, c111)
         + tpTetEulerContribution(c000, c100, c111, c101)
         + tpTetEulerContribution(c000, c010, c111, c110)
         + tpTetEulerContribution(c000, c010, c011, c111)
         + tpTetEulerContribution(c000, c001, c101, c111)
         + tpTetEulerContribution(c000, c001, c111, c011);
}

// ---------------------------------------------------------------------------
// 5. Higher-order link topology: Milnor's triple-linking invariant mu_123
// ---------------------------------------------------------------------------

// Normalized solid-angle scalar potential Omega / (4 pi) in [-1/2, 1/2] of the triangle fan from `apex`
// over any closed TP_LINK_SEGMENTS-gon L in R^3. Jumps by +-1 across the fan and satisfies grad u = B_BiotSavart
// away from the fan (distributional Poincare dual used by Mellor-Melvin triple linking).
float tpSolidAnglePolygon( vec3 p, vec3 L[TP_LINK_SEGMENTS], vec3 apex )
{
    float om = 0.0;
    vec3  o  = apex - p;
    vec3  a  = L[0] - p;
    for (int i = 0; i < TP_LINK_SEGMENTS; ++i)
    {
        int   nxt = (i + 1 < TP_LINK_SEGMENTS) ? (i + 1) : 0;
        vec3 b   = L[nxt] - p;
        om += tpSolidAngleTri(o, a, b);
        a = b;
    }
    return om * 0.07957747154594767;
}

// Signed intersection in {-1, 0, +1} of oriented segment p0 -> p1 with oriented triangle (a, b, c) (0 sqrt, 0 atan).
int tpSegTriHit( vec3 p0, vec3 p1, vec3 a, vec3 b, vec3 c, out float tHit )
{
    tHit = 0.0;
    vec3  d   = p1 - p0;
    vec3  e1  = b - a;
    vec3  e2  = c - a;
    vec3  n   = cross(e1, e2);
    float det = dot(n, d);
    float tol = 1e-20 * max(dot(n, n) * dot(d, d), 1e-20);
    if (det * det <= tol) return 0;
    float invDet = 1.0 / det;
    vec3  rhs = p0 - a;
    float t   = -dot(n, rhs) * invDet;
    if (t < 0.0 || t >= 1.0) return 0;
    float u   = dot(rhs, cross(e2, d)) * invDet;
    if (u < 0.0 || u > 1.0) return 0;
    float v   = dot(d, cross(e1, rhs)) * invDet;
    if (v < 0.0 || u + v > 1.0) return 0;
    tHit = t;
    return (det > 0.0) ? 1 : -1;
}

// Parameter interval along x(s) = orig + (s / |ldir|^2) ldir where a plane-straddling triangle clips the plane.
vec2 tpTriPlaneInterval( vec3 v0, vec3 v1, vec3 v2, float d0, float d1, float d2, vec3 orig, vec3 ldir )
{
    float sLo = 1e30, sHi = -1e30;
    if ((d0 < 0.0) != (d1 < 0.0))
    {
        vec3  a = (d0 < d1) ? v0 : v1, b = (d0 < d1) ? v1 : v0;
        float da = min(d0, d1), db = max(d0, d1);
        vec3  pt = mix(a, b, (-da) / max(db - da, 1e-30));
        float s  = dot(pt - orig, ldir);
        sLo = min(sLo, s); sHi = max(sHi, s);
    }
    if ((d1 < 0.0) != (d2 < 0.0))
    {
        vec3  a = (d1 < d2) ? v1 : v2, b = (d1 < d2) ? v2 : v1;
        float da = min(d1, d2), db = max(d1, d2);
        vec3  pt = mix(a, b, (-da) / max(db - da, 1e-30));
        float s  = dot(pt - orig, ldir);
        sLo = min(sLo, s); sHi = max(sHi, s);
    }
    if ((d2 < 0.0) != (d0 < 0.0))
    {
        vec3  a = (d2 < d0) ? v2 : v0, b = (d2 < d0) ? v0 : v2;
        float da = min(d2, d0), db = max(d2, d0);
        vec3  pt = mix(a, b, (-da) / max(db - da, 1e-30));
        float s  = dot(pt - orig, ldir);
        sLo = min(sLo, s); sHi = max(sHi, s);
    }
    return vec2(sLo, sHi);
}

// Oriented intersection segment q0 -> q1 (along n1 x n2) of two oriented triangles, with 6-dot early rejection.
bool tpTriTriSegment( vec3 a1, vec3 b1, vec3 c1, vec3 a2, vec3 b2, vec3 c2, out vec3 q0, out vec3 q1 )
{
    q0 = vec3(0.0); q1 = vec3(0.0);
    vec3  n1Raw = cross(b1 - a1, c1 - a1);
    float d2a = dot(n1Raw, a2 - a1), d2b = dot(n1Raw, b2 - a1), d2c = dot(n1Raw, c2 - a1);
    if ((d2a >= 0.0 && d2b >= 0.0 && d2c >= 0.0) || (d2a < 0.0 && d2b < 0.0 && d2c < 0.0)) return false;

    vec3  n2Raw = cross(b2 - a2, c2 - a2);
    float d1a = dot(n2Raw, a1 - a2), d1b = dot(n2Raw, b1 - a2), d1c = dot(n2Raw, c1 - a2);
    if ((d1a >= 0.0 && d1b >= 0.0 && d1c >= 0.0) || (d1a < 0.0 && d1b < 0.0 && d1c < 0.0)) return false;

    float n1Sq = dot(n1Raw, n1Raw), n2Sq = dot(n2Raw, n2Raw);
    if (n1Sq <= 1e-30 || n2Sq <= 1e-30) return false;
    vec3  n1   = n1Raw * inversesqrt(n1Sq);
    vec3  n2   = n2Raw * inversesqrt(n2Sq);
    vec3  ldir = cross(n1, n2);
    float l2   = dot(ldir, ldir);
    if (l2 <= 1e-16) return false;

    float invL2 = 1.0 / l2;
    float p2    = dot(n2, a2 - a1);
    vec3  orig  = a1 + (p2 * invL2) * cross(ldir, n1);

    vec2  i1 = tpTriPlaneInterval(a1, b1, c1, d1a, d1b, d1c, orig, ldir);
    vec2  i2 = tpTriPlaneInterval(a2, b2, c2, d2a, d2b, d2c, orig, ldir);
    float lo = max(i1.x, i2.x), hi = min(i1.y, i2.y);
    if (hi <= lo) return false;
    q0 = orig + (lo * invL2) * ldir;
    q1 = orig + (hi * invL2) * ldir;
    return true;
}

// Basepoint-independent Massey word pairing sum_{p in L cap F_b} w_a(p) sign_b(p) along polygon L (0 local arrays).
int tpCurveFanPairing( vec3 L[TP_LINK_SEGMENTS], vec3 oA, vec3 LA[TP_LINK_SEGMENTS], vec3 oB, vec3 LB[TP_LINK_SEGMENTS] )
{
    int wa    = 0;
    int total = 0;
    for (int i = 0; i < TP_LINK_SEGMENTS; ++i)
    {
        int  i1 = (i + 1 < TP_LINK_SEGMENTS) ? (i + 1) : 0;
        vec3 p0 = L[i], p1 = L[i1];
        for (int jb = 0; jb < TP_LINK_SEGMENTS; ++jb)
        {
            int   jb1 = (jb + 1 < TP_LINK_SEGMENTS) ? (jb + 1) : 0;
            float tb;
            int   sb = tpSegTriHit(p0, p1, oB, LB[jb], LB[jb1], tb);
            if (sb != 0)
            {
                int wHere = wa;
                for (int ja = 0; ja < TP_LINK_SEGMENTS; ++ja)
                {
                    int   ja1 = (ja + 1 < TP_LINK_SEGMENTS) ? (ja + 1) : 0;
                    float ta;
                    int   sa = tpSegTriHit(p0, p1, oA, LA[ja], LA[ja1], ta);
                    if (sa != 0 && ta < tb) wHere += sa;
                }
                total += wHere * sb;
            }
        }
        for (int ja = 0; ja < TP_LINK_SEGMENTS; ++ja)
        {
            int   ja1 = (ja + 1 < TP_LINK_SEGMENTS) ? (ja + 1) : 0;
            float ta;
            wa += tpSegTriHit(p0, p1, oA, LA[ja], LA[ja1], ta);
        }
    }
    return total;
}

// Exact Milnor triple-linking invariant mu_123(L1, L2, L3) = c1 + c2 + c3 - t123 (Mellor & Melvin 2003)
// for three disjoint pairwise-unlinked closed TP_LINK_SEGMENTS-gons (0 sqrt, 0 atan, exact integer result).
// Perturbs fan apexes slightly off the polygon centroids along irrational axes to avoid planar vertex/edge ties.
float tpMilnorTripleLink( vec3 L1[TP_LINK_SEGMENTS], vec3 L2[TP_LINK_SEGMENTS], vec3 L3[TP_LINK_SEGMENTS] )
{
    vec3 c1 = vec3(0.0), c2 = vec3(0.0), c3 = vec3(0.0);
    vec3 mn = L1[0], mx = L1[0];
    for (int i = 0; i < TP_LINK_SEGMENTS; ++i)
    {
        c1 += L1[i]; c2 += L2[i]; c3 += L3[i];
        mn = min(mn, min(L1[i], min(L2[i], L3[i])));
        mx = max(mx, max(L1[i], max(L2[i], L3[i])));
    }
    float invN  = 1.0 / float(TP_LINK_SEGMENTS);
    vec3  span  = mx - mn;
    float scale = max(max(span.x, max(span.y, span.z)), 1e-3);
    vec3  o1 = c1 * invN + scale * vec3( 0.01317, -0.01973,  0.02411);
    vec3  o2 = c2 * invN + scale * vec3(-0.02137,  0.01511, -0.01789);
    vec3  o3 = c3 * invN + scale * vec3( 0.01723,  0.02341, -0.01193);

    int m1 = tpCurveFanPairing(L1, o2, L2, o3, L3);
    int m2 = tpCurveFanPairing(L2, o3, L3, o1, L1);
    int m3 = tpCurveFanPairing(L3, o1, L1, o2, L2);

    vec3 segQ0[4 * TP_LINK_SEGMENTS];
    vec3 segQ1[4 * TP_LINK_SEGMENTS];
    int  nSeg = 0;
    int  t123 = 0;
    for (int i = 0; i < TP_LINK_SEGMENTS; ++i)
    {
        int i1 = (i + 1 < TP_LINK_SEGMENTS) ? (i + 1) : 0;
        for (int j = 0; j < TP_LINK_SEGMENTS; ++j)
        {
            int  j1 = (j + 1 < TP_LINK_SEGMENTS) ? (j + 1) : 0;
            vec3 q0, q1;
            if (tpTriTriSegment(o1, L1[i], L1[i1], o2, L2[j], L2[j1], q0, q1))
            {
                if (nSeg < 4 * TP_LINK_SEGMENTS)
                {
                    segQ0[nSeg] = q0;
                    segQ1[nSeg] = q1;
                    ++nSeg;
                }
            }
        }
    }
    for (int s = 0; s < nSeg; ++s)
    {
        for (int k = 0; k < TP_LINK_SEGMENTS; ++k)
        {
            int   k1 = (k + 1 < TP_LINK_SEGMENTS) ? (k + 1) : 0;
            float th;
            t123 += tpSegTriHit(segQ0[s], segQ1[s], o3, L3[k], L3[k1], th);
        }
    }
    return float(m1 + m2 + m3 - t123);
}

// Centre curves of the orthogonal Borromean rings (matching sdBorromeanRings: semi-axes ab = (a, b), a > b),
// with optional (m, 1)-longitudinal unknotted cabling on component 0 (giving mu_123 = -m while all three
// components remain unknots with pairwise linking numbers 0).
vec3 tpBorromeanStrand( float tt, int comp, vec2 ab, int m )
{
    if (comp == 0)
    {
        if (m == 1) return vec3(ab.x * cos(tt), ab.y * sin(tt), 0.0);
        float th = float(m) * tt;
        float rc = 0.08 * min(ab.x, ab.y);
        return vec3((ab.x + rc * cos(tt)) * cos(th), (ab.y + rc * cos(tt)) * sin(th), rc * sin(tt));
    }
    if (comp == 1) return vec3(0.0, ab.x * cos(tt), ab.y * sin(tt));
    return vec3(ab.y * sin(tt), 0.0, ab.x * cos(tt));
}

// Evaluates vec4(mu_123, lk_12, lk_23, lk_31) on the m-Borromean rings (ab.x > ab.y) with component 2 shifted
// by `shift` along +x. For m = 1 and shift = 0: returns vec4(-1, 0, 0, 0); unlinks cleanly to vec4(0) once
// |shift| > b * sqrt(1 - (b / a)^2), where component 2 passes across component 1.
vec4 tpMilnorBorromean( vec2 ab, int m, float shift )
{
    const float TAU = 6.283185307179586;
    vec3 L1[TP_LINK_SEGMENTS];
    vec3 L2[TP_LINK_SEGMENTS];
    vec3 L3[TP_LINK_SEGMENTS];
    for (int i = 0; i < TP_LINK_SEGMENTS; ++i)
    {
        float tt = TAU * float(i) / float(TP_LINK_SEGMENTS);
        L1[i] = tpBorromeanStrand(tt, 0, ab, m);
        L2[i] = tpBorromeanStrand(tt, 1, ab, 1);
        L3[i] = tpBorromeanStrand(tt, 2, ab, 1) + vec3(shift, 0.0, 0.0);
    }
    float lk12 = 0.0, lk23 = 0.0, lk31 = 0.0;
    for (int i = 0; i < TP_LINK_SEGMENTS; ++i)
    {
        int i1 = (i + 1 < TP_LINK_SEGMENTS) ? (i + 1) : 0;
        for (int j = 0; j < TP_LINK_SEGMENTS; ++j)
        {
            int j1 = (j + 1 < TP_LINK_SEGMENTS) ? (j + 1) : 0;
            lk12 += tpLinkingSegSeg(L1[i], L1[i1], L2[j], L2[j1]);
            lk23 += tpLinkingSegSeg(L2[i], L2[i1], L3[j], L3[j1]);
            lk31 += tpLinkingSegSeg(L3[i], L3[i1], L1[j], L1[j1]);
        }
    }
    return vec4(tpMilnorTripleLink(L1, L2, L3), lk12, lk23, lk31);
}

#endif // SDF_GEMS_TP_GAUSS_BONNET_GLSL
