// ============================================================================
// hiProjectiveMetric — Hilbert Projective Metric, Funk Weak Metrics & Geodesic Balls in Convex Polyhedra
//   hiChordDistanceCube, hiFunkDistanceCube, hiChordDistanceTetra, hiFunkDistanceTetra,
//   sdHilbertBallCube, sdgHilbertBallCube, sdHilbertBallOcta, sdgHilbertBallOcta,
//   sdHilbertBallTetra, sdgHilbertBallTetra, sdFunkBallTetra, sdgFunkBallTetra,
//   sdSymFunkBallTetra, sdgSymFunkBallTetra, sdHilbertBallAt, sdgHilbertBallAt
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  sdHilbertBallCube / sdHilbertBallOcta / sdHilbertBallTetra / sdFunkBallTetra / sdSymFunkBallTetra:
//         Exact Euclidean SDFs (`||grad d|| = 1` everywhere in R^3, 0/1 sqrt via Weyl-chamber folding).
//         sdHilbertBallAt: Exact inside and on all 18 face prisms of the off-center Hilbert polyhedron in a box;
//         1-Lipschitz conservative Euclidean bound outside edges/vertices.
// ============================================================================

#ifndef SDF_GEMS_HAS_HI_PROJECTIVE_METRIC
#define SDF_GEMS_HAS_HI_PROJECTIVE_METRIC

// ---------------------------------------------------------------------------
// 1. Chord-Intersection Hilbert & Funk Distances on a Box K = (-b, b)^3
//    For x, y in K, let q = y - x. The ray x + t q (t >= 0) pierces dK at t_+ = 1 / mu_+,
//    and the backward ray x - t q pierces dK at t_- = 1 / mu_-, where
//      mu_+(x, y) = max_k max( q_k / (b_k - x_k), -q_k / (b_k + x_k)),
//      mu_-(x, y) = max_k max(-q_k / (b_k - x_k),  q_k / (b_k + x_k)).
//    Evaluating reciprocal exit times avoids 0/0 ray-slab divisions on axis-aligned chords.
// ---------------------------------------------------------------------------

// Forward Funk weak distance d_F(x, y) = -ln(1 - mu_+(x, y)) in (-b, b)^3.
float hiFunkDistanceCube( in vec3 x, in vec3 y, in vec3 boxHalf )
{
    vec3 q   = y - x;
    vec3 dp  = max(boxHalf - x, 1e-6);
    vec3 dm  = max(boxHalf + x, 1e-6);
    vec3 mup = max(q / dp, -q / dm);
    float mp = clamp(max(mup.x, max(mup.y, mup.z)), 0.0, 1.0 - 1e-6);
    return -log(1.0 - mp);
}

// Hilbert projective distance d_H(x, y) = 0.5 * ln((1 + mu_-(x, y)) / (1 - mu_+(x, y))) in (-b, b)^3.
float hiChordDistanceCube( in vec3 x, in vec3 y, in vec3 boxHalf )
{
    vec3 q   = y - x;
    vec3 dp  = max(boxHalf - x, 1e-6);
    vec3 dm  = max(boxHalf + x, 1e-6);
    vec3 mup = max(q / dp, -q / dm);
    vec3 mum = max(-q / dp, q / dm);
    float mp = clamp(max(mup.x, max(mup.y, mup.z)), 0.0, 1.0 - 1e-6);
    float mm = max(0.0, max(mum.x, max(mum.y, mum.z)));
    return 0.5 * log((1.0 + mm) / (1.0 - mp));
}

// ---------------------------------------------------------------------------
// 2. Chord-Intersection Hilbert & Funk Distances on a Regular Tetrahedron K
//    K = { p : max_i (v_i . p) < s } with v_i in {(1,1,-1), (1,-1,1), (-1,1,1), (-1,-1,-1)}.
// ---------------------------------------------------------------------------
float hiFunkDistanceTetra( in vec3 x, in vec3 y, in float s )
{
    vec3 q = y - x;
    vec4 hx = max(vec4(s) - vec4(x.x + x.y - x.z, x.x - x.y + x.z, -x.x + x.y + x.z, -x.x - x.y - x.z), 1e-6);
    vec4 dq = vec4(q.x + q.y - q.z, q.x - q.y + q.z, -q.x + q.y + q.z, -q.x - q.y - q.z);
    vec4 r  = dq / hx;
    float mp = clamp(max(max(r.x, r.y), max(r.z, r.w)), 0.0, 1.0 - 1e-6);
    return -log(1.0 - mp);
}

float hiChordDistanceTetra( in vec3 x, in vec3 y, in float s )
{
    vec3 q = y - x;
    vec4 hx = max(vec4(s) - vec4(x.x + x.y - x.z, x.x - x.y + x.z, -x.x + x.y + x.z, -x.x - x.y - x.z), 1e-6);
    vec4 dq = vec4(q.x + q.y - q.z, q.x - q.y + q.z, -q.x + q.y + q.z, -q.x - q.y - q.z);
    vec4 rp = dq / hx;
    vec4 rm = -rp;
    float mp = clamp(max(max(rp.x, rp.y), max(rp.z, rp.w)), 0.0, 1.0 - 1e-6);
    float mm = max(0.0, max(max(rm.x, rm.y), max(rm.z, rm.w)));
    return 0.5 * log((1.0 + mm) / (1.0 - mp));
}

// ---------------------------------------------------------------------------
// 3. Exact Euclidean SDFs of Origin-Centered Hilbert Balls in Symmetric Bodies (Cube & Octahedron)
//    When K = -K, mu_+(p) = mu_-(p) = mu_K(p), so d_H(0, p) = artanh(mu_K(p)) <= RH
//    iff mu_K(p) <= tanh(RH): the Hilbert ball B_H(0, RH) is an exact homothety tanh(RH) * K.
// ---------------------------------------------------------------------------
vec4 sdgHilbertBallCube( in vec3 p, in vec3 boxHalf, in float RH )
{
    vec3  b  = boxHalf * tanh(max(RH, 0.0));
    vec3  s  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3  q  = abs(p) - b;
    vec3  qp = max(q, 0.0);
    float l2 = dot(qp, qp);
    if (l2 > 0.0)
    {
        float invL = inversesqrt(l2);
        return vec4(l2 * invL, s * qp * invL);
    }
    float m = max(q.x, max(q.y, q.z));
    vec3  g = (q.x >= q.y && q.x >= q.z) ? vec3(1.0, 0.0, 0.0)
            : ((q.y >= q.z) ? vec3(0.0, 1.0, 0.0) : vec3(0.0, 0.0, 1.0));
    return vec4(m, s * g);
}

float sdHilbertBallCube( in vec3 p, in vec3 boxHalf, in float RH )
{
    vec3 b = boxHalf * tanh(max(RH, 0.0));
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

vec4 sdgHilbertBallOcta( in vec3 p, in float s, in float RH )
{
    float r  = s * tanh(max(RH, 0.0));
    vec3  s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3  q  = abs(p);
    bool sw1 = (q.x < q.y); if (sw1) q.xy = q.yx;
    bool sw2 = (q.y < q.z); if (sw2) q.yz = q.zy;
    bool sw3 = (q.x < q.y); if (sw3) q.xy = q.yx;

    float m  = q.x + q.y + q.z - r;
    float a  = max(m - 3.0 * q.z, 0.0);
    float b  = max(q.x - q.y - r, 0.0);
    float sm = (m < 0.0) ? -1.0 : 1.0;
    float d2 = (2.0 * m * m + a * a + 3.0 * b * b) * (1.0 / 6.0);
    float d  = sm * sqrt(d2);
    vec3  g  = (d2 > 1e-30)
             ? (vec3(2.0 * m + a + 3.0 * b, 2.0 * m + a - 3.0 * b, 2.0 * m - 2.0 * a) / (6.0 * d))
             : vec3(0.57735027);
    if (sw3) g.xy = g.yx;
    if (sw2) g.yz = g.zy;
    if (sw1) g.xy = g.yx;
    return vec4(d, s0 * g);
}

float sdHilbertBallOcta( in vec3 p, in float s, in float RH )
{
    float r = s * tanh(max(RH, 0.0));
    p = abs(p);
    float m = p.x + p.y + p.z - r, z = min(p.x, min(p.y, p.z));
    float a = max(m - 3.0 * z, 0.0);
    float b = max(2.0 * max(p.x, max(p.y, p.z)) + z - m - 2.0 * r, 0.0);
    return ((m < 0.0) ? -1.0 : 1.0) * sqrt((2.0 * m * m + a * a + 3.0 * b * b) * (1.0 / 6.0));
}

// ---------------------------------------------------------------------------
// 4. Exact Euclidean SDF of the Hilbert Ball in a Regular Tetrahedron (sdHilbertBallTetra)
//    In the T_d Weyl chamber x >= y >= |z|:
//      mu_+(p) = (x + y - z) / s,   mu_-(p) = (x + y + z) / s.
//    Thus d_H(0, p) <= RH  <=>  mu_+(p) + exp(-2 RH) mu_-(p) <= 1 - exp(-2 RH)
//                          <=>  x + y - tau * z <= s * tau,   tau = tanh(RH).
//    Across R^3 this single Weyl-chamber half-space unfolds into a 12-faced Deltoidal Dodecahedron
//    that morphs continuously from the Catalan Rhombic Dodecahedron (RH -> 0) into the Regular
//    Tetrahedron (RH -> inf). Its exact Euclidean SDF needs 0 sqrts in the interior & 12 face prisms,
//    and 1 sqrt on the two active Weyl ridge segments unified by sgn = +/-1, a = 1 - sgn * tau.
// ---------------------------------------------------------------------------
float sdHilbertBallTetra( in vec3 p, in float s, in float RH )
{
    float tau = tanh(max(RH, 0.0));
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if (p.x < p.y) p.xy = p.yx;
    if (p.y < p.z) p.yz = p.zy;
    if (p.x < p.y) p.xy = p.yx;
    p.z *= sig;

    float st = s * tau;
    float n2 = 2.0 + tau * tau;
    float m  = p.x + p.y - tau * p.z - st;
    if (m <= 0.0 || ((1.0 + tau) * m <= n2 * (p.y - p.z) && (1.0 - tau) * m <= n2 * (p.y + p.z)))
        return m * inversesqrt(n2);

    float sgn = ((1.0 + tau) * m > n2 * (p.y - p.z)) ? 1.0 : -1.0;
    float a   = 1.0 - sgn * tau;
    float t   = clamp((p.y + sgn * p.z - a * (p.x - st)) / (2.0 + a * a), 0.0, st / (1.0 + a));
    return length(vec3(p.x - st + a * t, p.y - t, p.z - sgn * t));
}

vec4 sdgHilbertBallTetra( in vec3 p, in float s, in float RH )
{
    float tau = tanh(max(RH, 0.0));
    vec3  sg  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    float sig = sg.x * sg.y * sg.z;
    p = abs(p);
    bool s1 = (p.x < p.y); if (s1) p.xy = p.yx;
    bool s2 = (p.y < p.z); if (s2) p.yz = p.zy;
    bool s3 = (p.x < p.y); if (s3) p.xy = p.yx;
    p.z *= sig;

    float st = s * tau;
    float n2 = 2.0 + tau * tau;
    float m  = p.x + p.y - tau * p.z - st;
    float d;
    vec3  g;
    if (m <= 0.0 || ((1.0 + tau) * m <= n2 * (p.y - p.z) && (1.0 - tau) * m <= n2 * (p.y + p.z)))
    {
        float invN = inversesqrt(n2);
        d = m * invN;
        g = vec3(invN, invN, -tau * invN);
    }
    else
    {
        float sgn  = ((1.0 + tau) * m > n2 * (p.y - p.z)) ? 1.0 : -1.0;
        float a    = 1.0 - sgn * tau;
        float t    = clamp((p.y + sgn * p.z - a * (p.x - st)) / (2.0 + a * a), 0.0, st / (1.0 + a));
        vec3  v    = vec3(p.x - st + a * t, p.y - t, p.z - sgn * t);
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = (v2 > 1e-30) ? (v * invV) : vec3(1.0, 0.0, 0.0);
    }

    g.z *= sig;
    if (s3) g.xy = g.yx;
    if (s2) g.yz = g.zy;
    if (s1) g.xy = g.yx;
    return vec4(d, sg * g);
}

// ---------------------------------------------------------------------------
// 5. Exact Euclidean SDFs of Forward & Symmetrized Max-Funk Balls in a Tetrahedron
//    - Forward Funk ball B_F^+(0, RF) = { p : -ln(1 - mu_+(p)) <= RF } is the homothetic
//      regular tetrahedron of size s * (1 - exp(-RF)).
//    - Symmetrized Max-Funk ball B_F^sym(0, RF) = B_F^+(0, RF) cap B_F^-(0, RF)
//      = { p : mu_+(p) <= 1 - exp(-RF) and mu_-(p) <= exp(RF) - 1 } is an 8-faced Truncated
//      Tetrahedron that morphs continuously from the Regular Octahedron (RF -> 0), through the
//      Archimedean Truncated Tetrahedron (at RF = ln(5/3)), into the Regular Tetrahedron (RF >= ln 3).
// ---------------------------------------------------------------------------
float sdFunkBallTetra( in vec3 p, in float s, in float RF )
{
    float sf  = s * (1.0 - exp(-max(RF, 0.0)));
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if (p.x < p.y) p.xy = p.yx;
    if (p.y < p.z) p.yz = p.zy;
    if (p.x < p.y) p.xy = p.yx;
    p.z *= sig;

    float m = p.x + p.y - p.z - sf;
    if (m <= 1.5 * (p.y - p.z))
        return m * 0.57735027;
    float t = min(0.5 * (p.y + p.z), sf);
    return length(vec3(p.x - sf, p.y - t, p.z - t));
}

vec4 sdgFunkBallTetra( in vec3 p, in float s, in float RF )
{
    float sf  = s * (1.0 - exp(-max(RF, 0.0)));
    vec3  sg  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    float sig = sg.x * sg.y * sg.z;
    p = abs(p);
    bool s1 = (p.x < p.y); if (s1) p.xy = p.yx;
    bool s2 = (p.y < p.z); if (s2) p.yz = p.zy;
    bool s3 = (p.x < p.y); if (s3) p.xy = p.yx;
    p.z *= sig;

    float m = p.x + p.y - p.z - sf;
    float d;
    vec3  g;
    if (m <= 1.5 * (p.y - p.z))
    {
        d = m * 0.57735027;
        g = vec3(0.57735027, 0.57735027, -0.57735027);
    }
    else
    {
        float t    = min(0.5 * (p.y + p.z), sf);
        vec3  v    = vec3(p.x - sf, p.y - t, p.z - t);
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = (v2 > 1e-30) ? (v * invV) : vec3(1.0, 0.0, 0.0);
    }

    g.z *= sig;
    if (s3) g.xy = g.yx;
    if (s2) g.yz = g.zy;
    if (s1) g.xy = g.yx;
    return vec4(d, sg * g);
}

float sdSymFunkBallTetra( in vec3 p, in float s, in float RF )
{
    float rf  = max(RF, 0.0);
    float a   = s * (1.0 - exp(-rf));
    float b   = s * (exp(rf) - 1.0);
    float c   = min(0.5 * (b - a), a);

    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if (p.x < p.y) p.xy = p.yx;
    if (p.y < p.z) p.yz = p.zy;
    if (p.x < p.y) p.xy = p.yx;
    p.z *= sig;

    float d_hex  = (p.x + p.y - p.z - a) * 0.57735027;
    float d_tri  = (p.x + p.y + p.z - (a + 2.0 * c)) * 0.57735027;
    float d_face = max(d_hex, d_tri);
    bool  in_hex = (d_hex >= d_tri) && (p.x + p.y + 2.0 * p.z <= a + 3.0 * c) && (2.0 * p.x - p.y + p.z <= 2.0 * a);
    bool  in_tri = (d_tri >= d_hex) && (p.x + p.y - 2.0 * p.z <= a - c);
    if (d_face <= 0.0 || in_hex || in_tri)
        return d_face;

    float t1 = max(0.5 * (a - c - (p.x - p.y)), 0.0);
    float u2 = min(0.5 * (p.y + p.z), c);
    vec3  v  = ((p.x - p.y - (a - c)) <= (p.y + p.z - 2.0 * c))
             ? vec3(p.x - a + t1, p.y - c - t1, p.z - c)
             : vec3(p.x - a, p.y - u2, p.z - u2);
    return length(v);
}

vec4 sdgSymFunkBallTetra( in vec3 p, in float s, in float RF )
{
    float rf  = max(RF, 0.0);
    float a   = s * (1.0 - exp(-rf));
    float b   = s * (exp(rf) - 1.0);
    float c   = min(0.5 * (b - a), a);

    vec3  sg  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    float sig = sg.x * sg.y * sg.z;
    p = abs(p);
    bool s1 = (p.x < p.y); if (s1) p.xy = p.yx;
    bool s2 = (p.y < p.z); if (s2) p.yz = p.zy;
    bool s3 = (p.x < p.y); if (s3) p.xy = p.yx;
    p.z *= sig;

    float d_hex  = (p.x + p.y - p.z - a) * 0.57735027;
    float d_tri  = (p.x + p.y + p.z - (a + 2.0 * c)) * 0.57735027;
    float d_face = max(d_hex, d_tri);
    bool  in_hex = (d_hex >= d_tri) && (p.x + p.y + 2.0 * p.z <= a + 3.0 * c) && (2.0 * p.x - p.y + p.z <= 2.0 * a);
    bool  in_tri = (d_tri >= d_hex) && (p.x + p.y - 2.0 * p.z <= a - c);
    float d;
    vec3  g;
    if (d_face <= 0.0 || in_hex || in_tri)
    {
        d = d_face;
        g = (d_tri > d_hex) ? vec3(0.57735027) : vec3(0.57735027, 0.57735027, -0.57735027);
    }
    else
    {
        float t1   = max(0.5 * (a - c - (p.x - p.y)), 0.0);
        float u2   = min(0.5 * (p.y + p.z), c);
        vec3  v    = ((p.x - p.y - (a - c)) <= (p.y + p.z - 2.0 * c))
                   ? vec3(p.x - a + t1, p.y - c - t1, p.z - c)
                   : vec3(p.x - a, p.y - u2, p.z - u2);
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = (v2 > 1e-30) ? (v * invV) : vec3(1.0, 0.0, 0.0);
    }

    g.z *= sig;
    if (s3) g.xy = g.yx;
    if (s2) g.yz = g.zy;
    if (s1) g.xy = g.yx;
    return vec4(d, sg * g);
}

// ---------------------------------------------------------------------------
// 6. Off-Center Hilbert Ball in a Box K = (-boxHalf, boxHalf)^3 (sdHilbertBallAt)
//    For an off-center basepoint c in K, B_H(c, RH) = { p : mu_{c,+}(p) + alpha * mu_{c,-}(p) <= 1 - alpha }
//    with alpha = exp(-2 RH) is an 18-faced convex polyhedron (6 axial faces + 12 edge-bevel faces
//    where the forward ray and backward ray from c exit through perpendicular faces of the box).
//    Taking the maximum over the 6 axial planes and 12 bevel planes yields a closed-form 1-Lipschitz
//    Euclidean distance field that is EXACT everywhere inside B_H(c, RH) and on all 18 face prisms,
//    and a conservative Euclidean bound outside the edges and vertices.
// ---------------------------------------------------------------------------
vec4 sdgHilbertBallAt( in vec3 p, in vec3 c, in float RH, in vec3 boxHalf )
{
    float al  = exp(-2.0 * max(RH, 0.0));
    float h   = 1.0 - al;
    vec3  q   = p - c;
    vec3  ifp = 1.0 / max(boxHalf - c, 1e-6);
    vec3  ifm = 1.0 / max(boxHalf + c, 1e-6);
    vec3  ibp = al * ifm;
    vec3  ibm = al * ifp;

    // 6 axial faces (+x, -x, +y, -y, +z, -z):
    vec3  dp = q - h / (ifp + ibp);
    vec3  dm = -q - h / (ifm + ibm);
    float d  = dp.x;
    vec3  g  = vec3(1.0, 0.0, 0.0);
    if (dm.x > d) { d = dm.x; g = vec3(-1.0, 0.0, 0.0); }
    if (dp.y > d) { d = dp.y; g = vec3(0.0, 1.0, 0.0); }
    if (dm.y > d) { d = dm.y; g = vec3(0.0, -1.0, 0.0); }
    if (dp.z > d) { d = dp.z; g = vec3(0.0, 0.0, 1.0); }
    if (dm.z > d) { d = dm.z; g = vec3(0.0, 0.0, -1.0); }

    // 24 pairwise bevel planes (forward on signed axis i, backward on signed axis j != i):
    #define SDF_GEMS_HI_BEVEL4(QI, QJ, IFP_I, IFM_I, IBP_J, IBM_J, VEC_G) \
    { \
        float n1 = IFP_I, m1 = IBP_J, inv1 = inversesqrt(n1 * n1 + m1 * m1); \
        float d1 = ( n1 * QI + m1 * QJ - h) * inv1; \
        if (d1 > d) { float sx = n1 * inv1, sy = m1 * inv1; d = d1; g = VEC_G; } \
        float n2 = IFP_I, m2 = IBM_J, inv2 = inversesqrt(n2 * n2 + m2 * m2); \
        float d2 = ( n2 * QI - m2 * QJ - h) * inv2; \
        if (d2 > d) { float sx = n2 * inv2, sy = -m2 * inv2; d = d2; g = VEC_G; } \
        float n3 = IFM_I, m3 = IBP_J, inv3 = inversesqrt(n3 * n3 + m3 * m3); \
        float d3 = (-n3 * QI + m3 * QJ - h) * inv3; \
        if (d3 > d) { float sx = -n3 * inv3, sy = m3 * inv3; d = d3; g = VEC_G; } \
        float n4 = IFM_I, m4 = IBM_J, inv4 = inversesqrt(n4 * n4 + m4 * m4); \
        float d4 = (-n4 * QI - m4 * QJ - h) * inv4; \
        if (d4 > d) { float sx = -n4 * inv4, sy = -m4 * inv4; d = d4; g = VEC_G; } \
    }
    SDF_GEMS_HI_BEVEL4(q.x, q.y, ifp.x, ifm.x, ibp.y, ibm.y, vec3(sx, sy, 0.0))
    SDF_GEMS_HI_BEVEL4(q.y, q.x, ifp.y, ifm.y, ibp.x, ibm.x, vec3(sy, sx, 0.0))
    SDF_GEMS_HI_BEVEL4(q.x, q.z, ifp.x, ifm.x, ibp.z, ibm.z, vec3(sx, 0.0, sy))
    SDF_GEMS_HI_BEVEL4(q.z, q.x, ifp.z, ifm.z, ibp.x, ibm.x, vec3(sy, 0.0, sx))
    SDF_GEMS_HI_BEVEL4(q.y, q.z, ifp.y, ifm.y, ibp.z, ibm.z, vec3(0.0, sx, sy))
    SDF_GEMS_HI_BEVEL4(q.z, q.y, ifp.z, ifm.z, ibp.y, ibm.y, vec3(0.0, sy, sx))
    #undef SDF_GEMS_HI_BEVEL4

    return vec4(d, g);
}

float sdHilbertBallAt( in vec3 p, in vec3 c, in float RH, in vec3 boxHalf )
{
    return sdgHilbertBallAt(p, c, RH, boxHalf).x;
}

#endif // SDF_GEMS_HAS_HI_PROJECTIVE_METRIC
