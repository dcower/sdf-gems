// ============================================================================
// trAmoeba — Tropical Algebraic Geometry, Viro's Maslov Dequantization & Algebraic Amoebas
//   sdTropicalLine2D, sdgTropicalLine2D, sdTropicalLineTop2Bound2D,
//   sdTropicalPlane3D, sdgTropicalPlane3D, sdTropicalSkeleton3D,
//   sdTropicalEllipticCurve2D, sdTropicalQuartic3D,
//   sdAmoebaLine2D, sdgAmoebaLine2D, sdAmoebaPlane3D, sdgAmoebaPlane3D, sdAmoebaOctahedron3D
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  sdTropicalLine2D, sdgTropicalLine2D, sdTropicalPlane3D, sdgTropicalPlane3D,
//         sdTropicalSkeleton3D, sdTropicalEllipticCurve2D: Exact Euclidean distance to the tropical
//         hypersurface (corner locus / spine) by the Exact Tropical Spine Theorem (De Morgan convex-vacuum duality).
//         sdTropicalLineTop2Bound2D: Naive unnormalized top-2 bisector formula (P_(1)-P_(2))/||w_(1)-w_(2)||;
//         shown to be a discontinuous UPPER bound (overestimator by up to sqrt(2)) near tropical vertices.
//         sdTropicalQuartic3D: Exact 1-Lipschitz conservative bound for a thickened bounded tropical quartic K3 shell.
//         sdAmoebaLine2D, sdgAmoebaLine2D, sdAmoebaPlane3D, sdgAmoebaPlane3D, sdAmoebaOctahedron3D:
//         Exact zero-level set of the phase-free Algebraic Amoeba A_h = {P_(1) <= h ln sum_{j>=2} exp(P_(j)/h)}
//         and globally 1-Lipschitz conservative signed distance bound (normalized by Newton polytope diameter).
// ============================================================================

#ifndef SDF_GEMS_HAS_TR_AMOEBA
#define SDF_GEMS_HAS_TR_AMOEBA

// ---------------------------------------------------------------------------
// 1. 2D Tropical Line: P_trop(p) = max(a.x, a.y + p.x, a.z + p.y)
//    Weights: w_0 = (0,0), w_1 = (1,0), w_2 = (0,1).
//    Pairwise weight distances: ||w_1 - w_0|| = 1, ||w_2 - w_0|| = 1, ||w_1 - w_2|| = sqrt(2).
//
//    Exact Tropical Spine Theorem:
//    Every tropical vacuum cell C_i = {p : m_i(p) > m_j(p) for all j != i} is a convex polyhedron,
//    and R^n \ T(P) = U_i C_i. By De Morgan's law, R^n \ C_i = U_{j != i} H_{ij}^-, so for p in C_(1),
//        dist(p, T(P)) = dist(p, dC_(1)) = min_{j != (1)} (m_(1)(p) - m_j(p)) / ||w_(1) - w_j||
//    is the EXACT Euclidean distance everywhere in R^n without any vertex projections.
// ---------------------------------------------------------------------------

vec3 sdgTropicalLine2D( in vec2 p, in vec3 a )
{
    const float INV_SQRT2 = 0.70710678118654752440;
    float m0 = a.x;
    float m1 = a.y + p.x;
    float m2 = a.z + p.y;

    if (m0 >= m1 && m0 >= m2)
    {
        float d1 = m0 - m1;
        float d2 = m0 - m2;
        return (d1 <= d2) ? vec3(d1, -1.0, 0.0) : vec3(d2, 0.0, -1.0);
    }
    if (m1 >= m0 && m1 >= m2)
    {
        float d0 = m1 - m0;
        float d2 = (m1 - m2) * INV_SQRT2;
        return (d0 <= d2) ? vec3(d0, 1.0, 0.0) : vec3(d2, INV_SQRT2, -INV_SQRT2);
    }
    float d0 = m2 - m0;
    float d1 = (m2 - m1) * INV_SQRT2;
    return (d0 <= d1) ? vec3(d0, 0.0, 1.0) : vec3(d1, -INV_SQRT2, INV_SQRT2);
}

float sdTropicalLine2D( in vec2 p, in vec3 a )
{
    return sdgTropicalLine2D(p, a).x;
}

// Naive top-2 bisector formula (P_(1) - P_(2)) / ||w_(1) - w_(2)||: sorts by algebraic value BEFORE
// dividing by ||w_i - w_j||. Discontinuous and overestimates by up to sqrt(2) near the tropical vertex.
float sdTropicalLineTop2Bound2D( in vec2 p, in vec3 a )
{
    const float INV_SQRT2 = 0.70710678118654752440;
    float m0 = a.x, m1 = a.y + p.x, m2 = a.z + p.y;
    if (m0 >= m1 && m0 >= m2) return m0 - max(m1, m2);
    if (m1 >= m0 && m1 >= m2) return (m0 >= m2) ? (m1 - m0) : (m1 - m2) * INV_SQRT2;
    return (m0 >= m1) ? (m2 - m0) : (m2 - m1) * INV_SQRT2;
}

// ---------------------------------------------------------------------------
// 2. 3D Tropical Plane: P_trop(p) = max(a.x, a.y + p.x, a.z + p.y, a.w + p.z)
//    2-skeleton T_2(P) (6 planar sectors meeting along 4 rays): exact via min_{j != (1)}.
//    1-skeleton T_1(P) (4-ray tropical tripod/tetrapod where P_(1) = P_(2) = P_(3)): exact via ray projection.
// ---------------------------------------------------------------------------

vec4 sdgTropicalPlane3D( in vec3 p, in vec4 a )
{
    const float INV_SQRT2 = 0.70710678118654752440;
    float m0 = a.x, m1 = a.y + p.x, m2 = a.z + p.y, m3 = a.w + p.z;

    if (m0 >= m1 && m0 >= m2 && m0 >= m3)
    {
        float d1 = m0 - m1, d2 = m0 - m2, d3 = m0 - m3;
        if (d1 <= d2 && d1 <= d3) return vec4(d1, -1.0, 0.0, 0.0);
        if (d2 <= d3)             return vec4(d2, 0.0, -1.0, 0.0);
        return vec4(d3, 0.0, 0.0, -1.0);
    }
    if (m1 >= m0 && m1 >= m2 && m1 >= m3)
    {
        float d0 = m1 - m0, d2 = (m1 - m2) * INV_SQRT2, d3 = (m1 - m3) * INV_SQRT2;
        if (d0 <= d2 && d0 <= d3) return vec4(d0, 1.0, 0.0, 0.0);
        if (d2 <= d3)             return vec4(d2, INV_SQRT2, -INV_SQRT2, 0.0);
        return vec4(d3, INV_SQRT2, 0.0, -INV_SQRT2);
    }
    if (m2 >= m0 && m2 >= m1 && m2 >= m3)
    {
        float d0 = m2 - m0, d1 = (m2 - m1) * INV_SQRT2, d3 = (m2 - m3) * INV_SQRT2;
        if (d0 <= d1 && d0 <= d3) return vec4(d0, 0.0, 1.0, 0.0);
        if (d1 <= d3)             return vec4(d1, -INV_SQRT2, INV_SQRT2, 0.0);
        return vec4(d3, 0.0, INV_SQRT2, -INV_SQRT2);
    }
    float d0 = m3 - m0, d1 = (m3 - m1) * INV_SQRT2, d2 = (m3 - m2) * INV_SQRT2;
    if (d0 <= d1 && d0 <= d2) return vec4(d0, 0.0, 0.0, 1.0);
    if (d1 <= d2)             return vec4(d1, -INV_SQRT2, 0.0, INV_SQRT2);
    return vec4(d2, 0.0, -INV_SQRT2, INV_SQRT2);
}

float sdTropicalPlane3D( in vec3 p, in vec4 a )
{
    return sdgTropicalPlane3D(p, a).x;
}

// Exact Euclidean distance to the codimension-2 tropical 1-skeleton T_1(P) = {P_(1) = P_(2) = P_(3)}
// (4 rays emanating from the tropical vertex v = (a.x - a.y, a.x - a.z, a.x - a.w) in directions
// (-1,0,0), (0,-1,0), (0,0,-1), and (1,1,1)). Here 3-way edge / 4-way vertex projection IS required.
float sdTropicalSkeleton3D( in vec3 p, in vec4 a )
{
    vec3  q  = p - (a.xxx - a.yzw);
    float d1 = length(vec3(q.x - min(q.x, 0.0), q.y, q.z));
    float d2 = length(vec3(q.x, q.y - min(q.y, 0.0), q.z));
    float d3 = length(vec3(q.x, q.y, q.z - min(q.z, 0.0)));
    float t4 = max((q.x + q.y + q.z) * (1.0 / 3.0), 0.0);
    float d4 = length(q - vec3(t4));
    return min(min(d1, d2), min(d3, d4));
}

// ---------------------------------------------------------------------------
// 3. 2D Tropical Elliptic Curve (Smooth Genus-1 Cubic in R^2):
//    Newton polytope = degree-3 triangle {(i, j) in Z_>=0^2 : i + j <= 3} (10 monomials),
//    with concave quadratic lifting a_{i,j} = -c * (i^2 + j^2 + i*j - 3*(i + j)) + b * delta_{(1,1)}.
//    The interior lattice point (1, 1) generates a compact hexagonal cycle (genus 1) with 9 infinite tentacles.
//    By the Exact Tropical Spine Theorem, the minimum over the 9 other monomials is 100% exact everywhere!
// ---------------------------------------------------------------------------

float sdTropicalEllipticCurve2D( in vec2 p, in float c, in float b )
{
    vec2  W[10];
    float M[10];
    W[0] = vec2(0.0, 0.0); M[0] = 0.0;
    W[1] = vec2(1.0, 0.0); M[1] = 2.0 * c + p.x;
    W[2] = vec2(0.0, 1.0); M[2] = 2.0 * c + p.y;
    W[3] = vec2(2.0, 0.0); M[3] = 2.0 * c + 2.0 * p.x;
    W[4] = vec2(1.0, 1.0); M[4] = 3.0 * c + b + p.x + p.y;
    W[5] = vec2(0.0, 2.0); M[5] = 2.0 * c + 2.0 * p.y;
    W[6] = vec2(3.0, 0.0); M[6] = 3.0 * p.x;
    W[7] = vec2(2.0, 1.0); M[7] = 2.0 * c + 2.0 * p.x + p.y;
    W[8] = vec2(1.0, 2.0); M[8] = 2.0 * c + p.x + 2.0 * p.y;
    W[9] = vec2(0.0, 3.0); M[9] = 3.0 * p.y;

    int   i1 = 0;
    float p1 = M[0];
    for (int i = 1; i < 10; i++)
    {
        if (M[i] > p1) { p1 = M[i]; i1 = i; }
    }

    vec2  w1 = W[i1];
    float d  = 1e30;
    for (int j = 0; j < 10; j++)
    {
        if (j != i1)
        {
            vec2 dw = w1 - W[j];
            d = min(d, (p1 - M[j]) * inversesqrt(dot(dw, dw)));
        }
    }
    return d;
}

// ---------------------------------------------------------------------------
// 4. 3D Tropical Quartic K3 Hypersurface (Symmetric Octahedral-Tetrahedral Spine):
//    In octant coordinates q = abs(p), uses the interior monomial w_0 = (0,0,0) (m_0 = r),
//    the 3 axial weights (1,0,0), (0,1,0), (0,0,1), and the octant diagonal weight (1,1,1) (m_4 = (q.x+q.y+q.z) - 0.35*r),
//    giving a bounded rhombic-dodecahedral / octahedral K3 core with radiating tropical sheets, thickened by th
//    and bounded inside an octahedral envelope of radius 1.45 * r.
// ---------------------------------------------------------------------------

float sdTropicalQuartic3D( in vec3 p, in float r, in float th )
{
    vec3  q = abs(p);
    vec3  W[5];
    float M[5];
    W[0] = vec3(0.0, 0.0, 0.0); M[0] = 0.55 * r;
    W[1] = vec3(1.0, 0.0, 0.0); M[1] = q.x;
    W[2] = vec3(0.0, 1.0, 0.0); M[2] = q.y;
    W[3] = vec3(0.0, 0.0, 1.0); M[3] = q.z;
    W[4] = vec3(1.0, 1.0, 1.0); M[4] = (q.x + q.y + q.z) - 0.85 * r;

    int   i1 = 0;
    float p1 = M[0];
    for (int i = 1; i < 5; i++)
    {
        if (M[i] > p1) { p1 = M[i]; i1 = i; }
    }

    vec3  w1 = W[i1];
    float d_spine = 1e30;
    for (int j = 0; j < 5; j++)
    {
        if (j != i1)
        {
            vec3 dw = w1 - W[j];
            d_spine = min(d_spine, (p1 - M[j]) * inversesqrt(dot(dw, dw)));
        }
    }
    float d_shell = d_spine - th;
    float d_env   = (q.x + q.y + q.z - 1.35 * r) * 0.57735026918962576451;
    return max(d_shell, d_env);
}

// ---------------------------------------------------------------------------
// 5. Viro's Maslov Dequantization & Algebraic Amoeba of a 2D Line:
//    A_h = { p in R^2 : exists theta in T^3 s.t. sum_{i=0}^2 exp((a_i + w_i . p)/h + i theta_i) = 0 }.
//    By the polygon triangle inequality, p in A_h iff exp(P_(1)/h) <= exp(P_(2)/h) + exp(P_(3)/h), i.e.,
//        g_h(p) = P_(1)(p) - h * ln(exp(P_(2)(p)/h) + exp(P_(3)(p)/h)) <= 0.
//    Each vacuum component V_i(h) = {g_h > 0 in C_i} is strictly convex because g_h is concave on C_i,
//    and ||grad g_h(p)|| = ||w_(1) - sum_{j!=1} sigma_j w_j|| <= D(Delta) = sqrt(2) everywhere, so
//    d(p) = g_h(p) / sqrt(2) is globally C^0 and globally 1-Lipschitz on all of R^2.
// ---------------------------------------------------------------------------

vec3 sdgAmoebaLine2D( in vec2 p, in vec3 a, in float h )
{
    const float INV_SQRT2 = 0.70710678118654752440;
    float m0 = a.x, m1 = a.y + p.x, m2 = a.z + p.y;
    vec2  w0 = vec2(0.0, 0.0), w1 = vec2(1.0, 0.0), w2 = vec2(0.0, 1.0);

    float p1 = m0, p2 = m1, p3 = m2;
    vec2  v1 = w0, v2 = w1, v3 = w2;
    if (p2 > p1) { float t = p1; p1 = p2; p2 = t; vec2 tv = v1; v1 = v2; v2 = tv; }
    if (p3 > p2) { float t = p2; p2 = p3; p3 = t; vec2 tv = v2; v2 = v3; v3 = tv; }
    if (p2 > p1) { float t = p1; p1 = p2; p2 = t; vec2 tv = v1; v1 = v2; v2 = tv; }

    float hm   = max(h, 1e-6);
    float e32  = exp(max((p3 - p2) / hm, -80.0));
    float lse  = p2 + hm * log(1.0 + e32);
    float g    = p1 - lse;
    float s3   = e32 / (1.0 + e32);
    vec2  grad = v1 - ((1.0 - s3) * v2 + s3 * v3);
    float glen = length(grad);
    vec2  n    = (glen > 1e-12) ? (grad / glen) : vec2(1.0, 0.0);
    return vec3(g * INV_SQRT2, n);
}

float sdAmoebaLine2D( in vec2 p, in vec3 a, in float h )
{
    return sdgAmoebaLine2D(p, a, h).x;
}

// ---------------------------------------------------------------------------
// 6. 3D Algebraic Amoeba of a Plane:
//    A_h = { p in R^3 : P_(1)(p) - h * ln(exp(P_(2)/h) + exp(P_(3)/h) + exp(P_(4)/h)) <= 0 }.
//    Newton polytope diameter D(Delta) = sqrt(2), so d(p) = g_h(p) / sqrt(2) is globally 1-Lipschitz.
// ---------------------------------------------------------------------------

vec4 sdgAmoebaPlane3D( in vec3 p, in vec4 a, in float h )
{
    const float INV_SQRT2 = 0.70710678118654752440;
    float m[4];
    vec3  w[4];
    m[0] = a.x;       w[0] = vec3(0.0, 0.0, 0.0);
    m[1] = a.y + p.x; w[1] = vec3(1.0, 0.0, 0.0);
    m[2] = a.z + p.y; w[2] = vec3(0.0, 1.0, 0.0);
    m[3] = a.w + p.z; w[3] = vec3(0.0, 0.0, 1.0);

    // Sort descending so m[0] >= m[1] >= m[2] >= m[3]
    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            if (m[j] > m[i])
            {
                float tm = m[i]; m[i] = m[j]; m[j] = tm;
                vec3  tw = w[i]; w[i] = w[j]; w[j] = tw;
            }
        }
    }

    float hm   = max(h, 1e-6);
    float e2   = exp(max((m[2] - m[1]) / hm, -80.0));
    float e3   = exp(max((m[3] - m[1]) / hm, -80.0));
    float sumE = 1.0 + e2 + e3;
    float lse  = m[1] + hm * log(sumE);
    float g    = m[0] - lse;

    vec3  w_avg = (w[1] + e2 * w[2] + e3 * w[3]) / sumE;
    vec3  grad  = w[0] - w_avg;
    float glen  = length(grad);
    vec3  n     = (glen > 1e-12) ? (grad / glen) : vec3(1.0, 0.0, 0.0);
    return vec4(g * INV_SQRT2, n);
}

float sdAmoebaPlane3D( in vec3 p, in vec4 a, in float h )
{
    return sdgAmoebaPlane3D(p, a, h).x;
}

// ---------------------------------------------------------------------------
// 7. 3D Octahedral Algebraic Amoeba (`sdAmoebaOctahedron3D`):
//    Symmetric 3D Amoeba on octant coordinates q = abs(p) with central hub a.x = 0.22 * r,
//    8 diagonal tentacles along (+-1, +-1, +-1), and 12 hyperbolic-logarithmic vacuum arches,
//    smoothly bounded inside an octahedral envelope of radius r. Globally 1-Lipschitz.
// ---------------------------------------------------------------------------

float sdAmoebaOctahedron3D( in vec3 p, in float r, in float h )
{
    vec3  q    = abs(p);
    float d_am = sdAmoebaPlane3D(q, vec4(0.22 * r, 0.0, 0.0, 0.0), h);
    float d_oc = (q.x + q.y + q.z - 1.35 * r) * 0.57735026918962576451;
    return max(d_am, d_oc);
}

#endif // SDF_GEMS_HAS_TR_AMOEBA
