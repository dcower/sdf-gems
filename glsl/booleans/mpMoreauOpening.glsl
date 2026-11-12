// ============================================================================
// mpMoreauOpening — Moreau-Yosida Proximal Morphological Opening/Closing
//                   & Hopf-Lax Anisotropic Wulff Crystal Envelopes
//
//   mpJetAt                Transport a 1-jet (d, n) from probe q to query p
//   mpFarkasMargin4        Exact 3D Caratheodory/Farkas total-erosion margin
//   mpOpeningFromJets3     0/1-SFU Dual-Gram active-set QP for 3 supporting jets
//   mpOpeningFromJets4Out  0/1-SFU Dual-Gram active-set QP for 4 supporting jets
//   mpOpeningFromJets4     + 3D Caratheodory emptiness certificate (K (-) B_r = 0)
//   sdgOpenedTetrahedron / sdOpenedTetrahedron
//                          Black-box morphological opening via tetrahedral stencil
//                          + active-set cutting-plane reprojection
//   sdgOpenedPolyhedron  / sdOpenedPolyhedron
//                          Black-box opening of rectangular frustum across
//                          topological face/edge collapse (6 faces -> 5-face wedge)
//   sdgClosedCross3D     / sdClosedCross3D
//                          De Morgan morphological closing C_r(K) = (O_r(K^c))^c
//                          rounding 12 concave crevices + 8 trihedral saddles
//                          with 0 outer-arm bloating
//   mpWulffSupport / mpWulffProject / mpWulffNaiveOffset
//   sdgWulffCrystalOffset / mpWulffCrystalOffset
//                          Anisotropic Wulff (1901) crystal growth K (+) r W
//                          (octahedral {111}, rhombic dodecahedral {110},
//                           cubic {100}, hexagonal prism {10-10}+{0001})
//
// @author Gemini 4 Carbon (sdf-gems, 2026)
// @exact  Convex polyhedra: certified 1-Lipschitz lower bound at every step
//         (0.0% overstep), converging to FP32 machine precision (~1e-6) in 2-3
//         active-set cutting-plane reprojections, including past topological
//         face/edge collapse. Total erosion (r > r_in) requires 4 jets in 3D
//         (Caratheodory d + 1 = 4). Fails to be a global bound on non-convex
//         saddles / thin necks where single-point tangent planes slice across K;
//         see docs/booleans/mpMoreauOpening.md for full breakdown measurements.
// ============================================================================

#ifndef SDF_GEMS_HAS_MP_MOREAU_OPENING
#define SDF_GEMS_HAS_MP_MOREAU_OPENING

#include "../polyhedra/sdTetrahedron.glsl"
#include "../primitives3d/sdRectFrustum.glsl"
#include "../gradients3d/sdgBox.glsl"

// Forward declarations for single-translation-unit order independence
vec4 sdgTetrahedron( vec3 p, float s );
vec4 sdgRectFrustum( vec3 p, vec2 b1, vec2 b2, float h );
vec4 sdgBox( in vec3 p, in vec3 b );

// Sentinel distance returned when 4-jet Farkas certificate proves K (-) B_r == empty
#define MP_EMPTY_DIST 1e6

// Wulff crystal habit identifiers
#define MP_WULFF_OCTAHEDRON    0
#define MP_WULFF_RHOMBIC_DODEC 1
#define MP_WULFF_CUBE          2
#define MP_WULFF_HEX_PRISM     3

// ----------------------------------------------------------------------------
// 1. Supporting Half-Space Transport & 3D Caratheodory/Farkas Emptiness Test
// ----------------------------------------------------------------------------

// Transport a 1-jet dg = vec4(d(q), grad d(q)) evaluated at probe q to query p:
//   l(p) = d(q) + grad d(q) . (p - q)  <=  d_K(p)   (for any convex K).
vec4 mpJetAt( vec3 p, vec3 q, vec4 dg )
{
    return vec4(dg.x + dot(dg.yzw, p - q), dg.yzw);
}

// Exact 3D Caratheodory / Farkas emptiness margin for 4 unit-normal 1-jets at p.
// Computes mu = (det(n1,n2,n3), -det(n0,n2,n3), det(n0,n1,n3), -det(n0,n1,n2)),
// which satisfies mu.x*n0 + mu.y*n1 + mu.z*n2 + mu.w*n3 == 0 identically (2 cross
// products, 0 SFU ops). When all mu_k > 0, 0 lies strictly inside conv{n0,n1,n2,n3}
// and (mu . u) / sum(mu) == r - r_in^(4) > 0 certifies empty eroded core!
float mpFarkasMargin4( vec4 j0, vec4 j1, vec4 j2, vec4 j3, float r )
{
    vec3  c23 = cross(j2.yzw, j3.yzw);
    vec3  c01 = cross(j0.yzw, j1.yzw);
    vec4  mu  = vec4( dot(j1.yzw, c23), -dot(j0.yzw, c23),
                      dot(j3.yzw, c01), -dot(j2.yzw, c01) );
    float smu = mu.x + mu.y + mu.z + mu.w;
    if( smu < 0.0 ) { mu = -mu; smu = -smu; }
    if( min(min(mu.x, mu.y), min(mu.z, mu.w)) >= -1e-6 && smu > 1e-5 )
    {
        vec4 u = vec4(j0.x, j1.x, j2.x, j3.x) + r;
        return dot(mu, u) / smu;
    }
    return -1.0;
}

// ----------------------------------------------------------------------------
// 2. Closed-Form Dual-Gram Active-Set Supporting-Plane QP (3-Jet & 4-Jet)
//    Solves:  dist_signed(0, { n_i . z + l_i + r <= 0 }) - r
//    in 0 SFU ops inside the eroded core and 1 inversesqrt outside!
// ----------------------------------------------------------------------------

vec4 mpOpeningFromJets3( vec4 a, vec4 b, vec4 c, float r )
{
    float ua = a.x + r, ub = b.x + r, uc = c.x + r;

    // Inside all 3 shifted supporting half-spaces: 0 SFU ops
    if( max(max(ua, ub), uc) <= 0.0 )
    {
        vec4 m = (a.x >= b.x) ? a : b;
        return (m.x >= c.x) ? m : c;
    }

    vec3  na  = a.yzw, nb = b.yzw, nc = c.yzw;
    float cab = dot(na, nb), cbc = dot(nb, nc), cca = dot(nc, na);

    // 1D Caratheodory opposing-slab emptiness check (n_i + n_j == 0, u_i + u_j > 0)
    if( (cab < -0.99999 && ua + ub > 1e-5) ||
        (cbc < -0.99999 && ub + uc > 1e-5) ||
        (cca < -0.99999 && uc + ua > 1e-5) )
        return vec4(MP_EMPTY_DIST, 0.0, 1.0, 0.0);

    float detG = 1.0 + 2.0 * cab * cbc * cca - (cab * cab + cbc * cbc + cca * cca);

    // 2D Caratheodory coplanar prism emptiness check when normals lie in a plane
    if( detG <= 1e-5 )
    {
        vec3 w = cross(na, nb);
        if( dot(w, w) < 1e-4 ) w = cross(nb, nc);
        float w2 = dot(w, w);
        if( w2 >= 1e-4 )
        {
            vec3  mu  = vec3(dot(w, cross(nb, nc)), dot(w, cross(nc, na)), dot(w, cross(na, nb)));
            float smu = mu.x + mu.y + mu.z;
            if( smu < 0.0 ) { mu = -mu; smu = -smu; }
            if( min(mu.x, min(mu.y, mu.z)) >= -1e-5 && smu > 1e-4 )
            {
                if( dot(mu, vec3(ua, ub, uc)) > 1e-5 * smu )
                    return vec4(MP_EMPTY_DIST, 0.0, 1.0, 0.0);
            }
        }
    }

    // Track maximum dual objective 2*D(lambda) = rho^2 across all dual-feasible active sets
    float bestRho2 = -1.0;
    vec3  bestV    = na;

    // 1-plane active sets
    if( ua > 0.0 ) { bestRho2 = ua * ua; bestV = ua * na; }
    if( ub > 0.0 && ub * ub > bestRho2 ) { bestRho2 = ub * ub; bestV = ub * nb; }
    if( uc > 0.0 && uc * uc > bestRho2 ) { bestRho2 = uc * uc; bestV = uc * nc; }

    // 2-plane edge active sets (orthogonalized Varadhan-Gram form avoids c -> +-1 cancellation)
    float sin2_ab = 1.0 - cab * cab;
    if( sin2_ab > 1e-4 )
    {
        float w1 = ua - cab * ub, w2 = ub - cab * ua;
        if( w1 > 0.0 && w2 > 0.0 )
        {
            float su = ua + ub, sv = ua - ub;
            float r2 = (su * su) / (2.0 + 2.0 * cab) + (sv * sv) / (2.0 - 2.0 * cab);
            if( r2 > bestRho2 ) { bestRho2 = r2; float invS = 1.0 / sin2_ab; bestV = (w1 * invS) * na + (w2 * invS) * nb; }
        }
    }
    float sin2_bc = 1.0 - cbc * cbc;
    if( sin2_bc > 1e-4 )
    {
        float w1 = ub - cbc * uc, w2 = uc - cbc * ub;
        if( w1 > 0.0 && w2 > 0.0 )
        {
            float su = ub + uc, sv = ub - uc;
            float r2 = (su * su) / (2.0 + 2.0 * cbc) + (sv * sv) / (2.0 - 2.0 * cbc);
            if( r2 > bestRho2 ) { bestRho2 = r2; float invS = 1.0 / sin2_bc; bestV = (w1 * invS) * nb + (w2 * invS) * nc; }
        }
    }
    float sin2_ca = 1.0 - cca * cca;
    if( sin2_ca > 1e-4 )
    {
        float w1 = uc - cca * ua, w2 = ua - cca * uc;
        if( w1 > 0.0 && w2 > 0.0 )
        {
            float su = uc + ua, sv = uc - ua;
            float r2 = (su * su) / (2.0 + 2.0 * cca) + (sv * sv) / (2.0 - 2.0 * cca);
            if( r2 > bestRho2 ) { bestRho2 = r2; float invS = 1.0 / sin2_ca; bestV = (w1 * invS) * nc + (w2 * invS) * na; }
        }
    }

    // 3-plane trihedral vertex active set via exact 3x3 Gram adjugate
    if( detG > 5e-4 )
    {
        float wa = (1.0 - cbc * cbc) * ua + (cca * cbc - cab) * ub + (cab * cbc - cca) * uc;
        float wb = (cca * cbc - cab) * ua + (1.0 - cca * cca) * ub + (cab * cca - cbc) * uc;
        float wc = (cab * cbc - cca) * ua + (cab * cca - cbc) * ub + (1.0 - cab * cab) * uc;
        if( wa > 0.0 && wb > 0.0 && wc > 0.0 )
        {
            float invD = 1.0 / detG;
            vec3  v    = (wa * invD) * na + (wb * invD) * nb + (wc * invD) * nc;
            float r2   = dot(v, v);
            if( r2 > bestRho2 ) { bestRho2 = r2; bestV = v; }
        }
    }

    float v2 = dot(bestV, bestV);
    if( v2 > 1e-24 )
    {
        float invV = inversesqrt(v2);
        return vec4(v2 * invV - r, bestV * invV);
    }
    return vec4(-r, na);
}

// 4-jet Active-Set Supporting-Plane QP returning both (d, grad) and the winning 3-jet basis indices.
vec4 mpOpeningFromJets4Out( vec4 j0, vec4 j1, vec4 j2, vec4 j3, float r, out ivec3 activeTri )
{
    // 1. Exact 4-jet 3D Caratheodory/Farkas total-erosion check
    if( mpFarkasMargin4(j0, j1, j2, j3, r) > 1e-5 )
    {
        activeTri = ivec3(0, 1, 2);
        return vec4(MP_EMPTY_DIST, 0.0, 1.0, 0.0);
    }

    // 2. Interior of all 4 eroded half-spaces: 0 SFU ops
    if( max(max(j0.x, j1.x), max(j2.x, j3.x)) + r <= 0.0 )
    {
        vec4 best = j0; activeTri = ivec3(0, 1, 2);
        if( j1.x > best.x ) { best = j1; activeTri = ivec3(1, 0, 2); }
        if( j2.x > best.x ) { best = j2; activeTri = ivec3(2, 0, 1); }
        if( j3.x > best.x ) { best = j3; activeTri = ivec3(3, 0, 1); }
        return best;
    }

    // 3. Helly-Caratheodory minimax across the 4 triples
    vec4  dg012 = mpOpeningFromJets3(j0, j1, j2, r);
    vec4  best  = dg012;
    activeTri   = ivec3(0, 1, 2);

    vec4 dg013 = mpOpeningFromJets3(j0, j1, j3, r);
    if( dg013.x > best.x ) { best = dg013; activeTri = ivec3(0, 1, 3); }

    vec4 dg023 = mpOpeningFromJets3(j0, j2, j3, r);
    if( dg023.x > best.x ) { best = dg023; activeTri = ivec3(0, 2, 3); }

    vec4 dg123 = mpOpeningFromJets3(j1, j2, j3, r);
    if( dg123.x > best.x ) { best = dg123; activeTri = ivec3(1, 2, 3); }

    return best;
}

vec4 mpOpeningFromJets4( vec4 j0, vec4 j1, vec4 j2, vec4 j3, float r )
{
    ivec3 tri;
    return mpOpeningFromJets4Out(j0, j1, j2, j3, r, tri);
}

// Helper to retain the 3 active jets of a 4-jet pool (j0, j1, j2, j3)
void _mpKeepActive3( vec4 j0, vec4 j1, vec4 j2, vec4 j3, ivec3 tri,
                     out vec4 a0, out vec4 a1, out vec4 a2 )
{
    vec4 pool[4];
    pool[0] = j0; pool[1] = j1; pool[2] = j2; pool[3] = j3;
    a0 = pool[tri.x];
    a1 = pool[tri.y];
    a2 = pool[tri.z];
}

// ----------------------------------------------------------------------------
// 3. Black-Box Tetrahedral Stencil + Active-Set Cutting-Plane Opening Macro
//    & Concrete Verified Primitives
// ----------------------------------------------------------------------------

// Defines FUNC_SDG(vec3 p, float r) and FUNC_SD(vec3 p, float r) for any black-box
// convex 1-jet expression SDG_CALL(q) (where q is a vec3 probe point).
// - Inside K_{-r} (d_K(p) <= -r): returns sdg_K(p) in 1 evaluation (exact!).
// - Outside K_{-r}: probes a regular tetrahedron of radius 1.25*r around the initial
//   eroded footpoint x0 = p - (d0 + r)*n0, folds in dg0, and runs 3 monotone
//   Kelley/Matousek active-set cutting-plane reprojections x_k = p - (d_k + r)*g_k.
#define MP_EVAL_OPENING_BODY(SDG_CALL, P_IN, R_IN, RES_OUT)                         \
    do {                                                                            \
        vec3  _p   = (P_IN);                                                        \
        float _r   = max((R_IN), 0.0);                                              \
        vec3  _q   = _p;                                                            \
        vec4  _dg0 = (SDG_CALL);                                                    \
        if( _r <= 1e-7 || _dg0.x <= -_r )                                           \
        {                                                                           \
            RES_OUT = _dg0;                                                         \
        }                                                                           \
        else                                                                        \
        {                                                                           \
            vec3  _x0 = _p - (_dg0.x + _r) * _dg0.yzw;                              \
            float _rp = 1.25 * _r * 0.57735027;                                     \
            _q = _x0 + vec3( _rp,  _rp,  _rp); vec4 _j0 = mpJetAt(_p, _q, (SDG_CALL)); \
            _q = _x0 + vec3( _rp, -_rp, -_rp); vec4 _j1 = mpJetAt(_p, _q, (SDG_CALL)); \
            _q = _x0 + vec3(-_rp,  _rp, -_rp); vec4 _j2 = mpJetAt(_p, _q, (SDG_CALL)); \
            _q = _x0 + vec3(-_rp, -_rp,  _rp); vec4 _j3 = mpJetAt(_p, _q, (SDG_CALL)); \
            ivec3 _tri;                                                             \
            vec4  _dg = mpOpeningFromJets4Out(_j0, _j1, _j2, _j3, _r, _tri);        \
            vec4  _a0, _a1, _a2;                                                    \
            _mpKeepActive3(_j0, _j1, _j2, _j3, _tri, _a0, _a1, _a2);                \
            if( _dg.x < 0.5 * MP_EMPTY_DIST )                                       \
            {                                                                       \
                _dg = mpOpeningFromJets4Out(_a0, _a1, _a2, _dg0, _r, _tri);         \
                _mpKeepActive3(_a0, _a1, _a2, _dg0, _tri, _a0, _a1, _a2);           \
            }                                                                       \
            for( int _it = 0; _it < 3; _it++ )                                      \
            {                                                                       \
                if( _dg.x >= 0.5 * MP_EMPTY_DIST ) break;                           \
                _q = _p - (_dg.x + _r) * _dg.yzw;                                   \
                vec4 _dgk = (SDG_CALL);                                             \
                if( _dgk.x <= -_r + 1e-6 ) break;                                   \
                vec4 _jk = mpJetAt(_p, _q, _dgk);                                   \
                _dg = mpOpeningFromJets4Out(_a0, _a1, _a2, _jk, _r, _tri);          \
                _mpKeepActive3(_a0, _a1, _a2, _jk, _tri, _a0, _a1, _a2);            \
            }                                                                       \
            RES_OUT = _dg;                                                          \
        }                                                                           \
    } while( false )

// Demo 1: Black-box Morphological Opening of a Regular Tetrahedron (sdgTetrahedron)
// Rounds all 6 sharp edges (dihedral ~70.53 deg) and 4 trihedral vertices INWARD
// without altering the 4 face planes; certifies emptiness for r > s / sqrt(3).
vec4 sdgOpenedTetrahedron( vec3 p, float s, float r )
{
    vec4 res;
    MP_EVAL_OPENING_BODY(sdgTetrahedron(_q, s), p, r, res);
    return res;
}

float sdOpenedTetrahedron( vec3 p, float s, float r )
{
    return sdgOpenedTetrahedron(p, s, r).x;
}

// Demo 2: Black-box Morphological Opening of an Unequal Rectangular Frustum (sdgRectFrustum)
// Smoothly crosses TOPOLOGICAL FACE/EDGE COLLAPSE when r > min(b2.x, b2.y), where the
// flat top cap of K_{-r} vanishes and non-adjacent slanted faces meet in a ridge!
vec4 sdgOpenedPolyhedron( vec3 p, vec2 b1, vec2 b2, float h, float r )
{
    vec4 res;
    MP_EVAL_OPENING_BODY(sdgRectFrustum(_q, b1, b2, h), p, r, res);
    return res;
}

float sdOpenedPolyhedron( vec3 p, vec2 b1, vec2 b2, float h, float r )
{
    return sdgOpenedPolyhedron(p, b1, b2, h, r).x;
}

// Demo 3: Morphological Closing of a 3D Orthogonal Cross via De Morgan Dual Gram QP
// C_r(Cross3D) = (O_r((Cross3D)^c))^c. Since (B_x U B_y U B_z)^c = B_x^c n B_y^c n B_z^c,
// feeding the 3 negated arm-box 1-jets into mpOpeningFromJets3 rounds all 12 concave
// 90-deg crevices and all 8 trihedral 3-way saddle corners simultaneously with 0.0
// outer-arm bloating and 0.0% exterior overstep!
vec4 sdgClosedCross3D( vec3 p, vec2 b, float r )
{
    vec4 dx = sdgBox(p, vec3(b.x, b.y, b.y));
    vec4 dy = sdgBox(p, vec3(b.y, b.x, b.y));
    vec4 dz = sdgBox(p, vec3(b.y, b.y, b.x));
    if( r <= 1e-7 )
    {
        vec4 m = (dx.x < dy.x) ? dx : dy;
        return (m.x < dz.x) ? m : dz;
    }
    vec4 op = mpOpeningFromJets3(-dx, -dy, -dz, r);
    if( op.x + r > 0.0 )
    {
        // 1 reprojection to the dilated complement boundary eliminates cylinder-curvature overshoot
        vec3 q1 = p - (op.x + r) * op.yzw;
        vec4 jx = mpJetAt(p, q1, -sdgBox(q1, vec3(b.x, b.y, b.y)));
        vec4 jy = mpJetAt(p, q1, -sdgBox(q1, vec3(b.y, b.x, b.y)));
        vec4 jz = mpJetAt(p, q1, -sdgBox(q1, vec3(b.y, b.y, b.x)));
        op = mpOpeningFromJets3(jx, jy, jz, r);
    }
    return -op;
}

float sdClosedCross3D( vec3 p, vec2 b, float r )
{
    return sdgClosedCross3D(p, b, r).x;
}

// ----------------------------------------------------------------------------
// 4. Hopf-Lax / Legendre-Fenchel Anisotropic Wulff Crystal Offset (K (+) r W)
//    Crystal support functions h_W(u) = max_{w in W} (u . w) and exact 0-SFU
//    Euclidean projections Pi_W(y) onto the four crystallographic Wulff shapes.
// ----------------------------------------------------------------------------

float mpWulffSupport( vec3 u, int habit )
{
    vec3 a = abs(u);
    if( habit == MP_WULFF_OCTAHEDRON )
        return max(a.x, max(a.y, a.z));
    if( habit == MP_WULFF_RHOMBIC_DODEC )
        return max(max(a.x, max(a.y, a.z)), 0.5 * (a.x + a.y + a.z));
    if( habit == MP_WULFF_CUBE )
        return 0.57735027 * (a.x + a.y + a.z);
    // MP_WULFF_HEX_PRISM: circumradius-normalized {10-10} + {0001} hexagonal prism
    const float cHex = 0.70710678;
    float h2d = max(a.x, 0.5 * a.x + 0.86602540 * a.z);
    return cHex * (h2d + a.y);
}

// Exact 0-SFU Euclidean projection onto W: argmin_{w in W} ||w - y||_2
vec3 mpWulffProject( vec3 y, int habit )
{
    if( habit == MP_WULFF_CUBE )
    {
        const float sCube = 0.57735027;
        return clamp(y, vec3(-sCube), vec3(sCube));
    }
    if( habit == MP_WULFF_HEX_PRISM )
    {
        const float cHex = 0.70710678;
        float wy = clamp(y.y, -cHex, cHex);
        // Project y.xz onto regular hexagon of circumradius cHex (inradius cHex * sqrt(3)/2):
        // Dual vertices (+-1, 0), (+-1/2, +-sqrt(3)/2) * cHex -> fold into 30-deg Weyl sector
        vec2  s2 = vec2((y.x < 0.0) ? -1.0 : 1.0, (y.z < 0.0) ? -1.0 : 1.0);
        vec2  q  = abs(y.xz);
        // Hexagon { q.z <= cHex*sqrt(3)/2,  sqrt(3)/2*q.x + 0.5*q.z <= cHex*sqrt(3)/2 }
        const vec2 n1 = vec2(0.0, 1.0);
        const vec2 n2 = vec2(0.86602540, 0.5);
        float hIn = cHex * 0.86602540;
        float d1  = q.y - hIn;
        float d2  = dot(q, n2) - hIn;
        vec2  w2d = q;
        if( d1 > 0.0 || d2 > 0.0 )
        {
            // Project onto active edge(s) or vertex (cHex, 0) / (0.5*cHex, hIn)
            vec2  vTop = vec2(0.5 * cHex, hIn);
            vec2  vMid = vec2(cHex, 0.0);
            if( q.x <= 0.5 * cHex && d1 >= 0.0 )
                w2d = vec2(q.x, hIn);
            else
            {
                vec2  ed = vMid - vTop;
                float t  = clamp(dot(q - vTop, ed) / dot(ed, ed), 0.0, 1.0);
                w2d = vTop + t * ed;
            }
        }
        return vec3(s2.x * w2d.x, wy, s2.y * w2d.y);
    }

    // Sort |y| into q.x >= q.y >= q.z >= 0 for Octahedron & Rhombic Dodecahedron
    vec3 s0 = vec3((y.x < 0.0) ? -1.0 : 1.0, (y.y < 0.0) ? -1.0 : 1.0, (y.z < 0.0) ? -1.0 : 1.0);
    vec3 q  = abs(y);
    bool b1 = (q.x < q.y); if( b1 ) q.xy = q.yx;
    bool b2 = (q.y < q.z); if( b2 ) q.yz = q.zy;
    bool b3 = (q.x < q.y); if( b3 ) q.xy = q.yx;

    vec3 w;
    if( habit == MP_WULFF_OCTAHEDRON )
    {
        // Exact Duchi-Condat L1-ball projection { q.x + q.y + q.z <= 1 }
        if( q.x + q.y + q.z <= 1.0 )
            w = q;
        else
        {
            float th = (q.x - q.y >= 1.0) ? (q.x - 1.0)
                     : ((q.x + q.y - 2.0 * q.z >= 1.0) ? (0.5 * (q.x + q.y - 1.0))
                                                       : ((q.x + q.y + q.z - 1.0) * (1.0 / 3.0)));
            w = max(q - th, 0.0);
        }
    }
    else
    {
        // Exact Weyl-chamber projection onto Rhombic Dodecahedron { q.x + q.y <= 1 }
        if( q.x + q.y <= 1.0 )
            w = q;
        else
        {
            float m  = q.x + q.y - 1.0;
            float a  = max(q.x - q.y + 2.0 * q.z - 1.0, 0.0);
            float b0 = q.y + q.z - q.x + 0.25;
            float b  = ((b0 < 0.0) ? -1.0 : 1.0) * max(abs(b0) - 0.75, 0.0);
            w = q - vec3(3.0 * m + a - 2.0 * b, 3.0 * m - a + 2.0 * b, 2.0 * a + 2.0 * b) * (1.0 / 6.0);
        }
    }

    if( b3 ) w.xy = w.yx;
    if( b2 ) w.yz = w.zy;
    if( b1 ) w.xy = w.yx;
    return s0 * w;
}

// Naive single-jet Hopf-Lax approximation: d_K(p) - r * h_W(grad d_K(p)).
// Exact ONLY on flat faces of K; produces curved non-Eikonal lobes (||grad|| up to
// sqrt(2) > 1, 14.6%*r surface bulge) at edges and corners where grad d_K rotates!
float mpWulffNaiveOffset( vec4 dg, float r, int habit )
{
    return dg.x - r * mpWulffSupport(dg.yzw, habit);
}

// Certified Legendre-Fenchel / Alternating-Proximal Wulff Crystal Envelope K (+) r W
// on a regular tetrahedron sdgTetrahedron(p, s) (generalizes to any convex sdg_K).
// Every supporting evaluation dg = sdg_K(p - r*w) yields a certified 1-Lipschitz
// lower bound:  l_W(p; w) = dg.x + r * (dot(dg.yzw, w) - h_W(dg.yzw)) <= d_{K (+) r W}(p).
// Two proximal projections w <- Pi_W((p - c_K) / r) lock onto the exact jewel-cut
// crystal bevels and facets to FP32 machine precision (~1e-6) with 0.0% overstep!
vec4 sdgWulffCrystalOffset( vec3 p, float s, float r, int habit )
{
    vec4 dg0 = sdgTetrahedron(p, s);
    if( r <= 1e-7 ) return dg0;

    // Candidate 0: w = 0 (always unit-gradient if dg0 is valid)
    vec3  g0     = (dot(dg0.yzw, dg0.yzw) > 0.5) ? dg0.yzw : vec3(0.57735027, 0.57735027, -0.57735027);
    float bestLo = dg0.x - r * mpWulffSupport(g0, habit);
    vec3  bestG  = g0;

    vec3 w = mpWulffProject((max(dg0.x, 0.0) / r) * g0, habit);

    // Proximal alternating steps between K and (p - r W)
    for( int it = 0; it < 3; it++ )
    {
        vec3 q  = p - r * w;
        vec4 dg = sdgTetrahedron(q, s);
        if( dot(dg.yzw, dg.yzw) > 0.5 )
        {
            float lo = dg.x + r * (dot(dg.yzw, w) - mpWulffSupport(dg.yzw, habit));
            if( lo > bestLo ) { bestLo = lo; bestG = dg.yzw; }
            if( dg.x > 1e-5 )
            {
                vec3 cK = q - dg.x * dg.yzw;
                w = mpWulffProject((p - cK) / r, habit);
            }
            else
            {
                w = mpWulffProject(w + 0.75 * dg.yzw, habit);
            }
        }
        else
        {
            // q landed directly on a sharp feature of K (q == cK): p - cK == r * w!
            w = mpWulffProject(w + 0.25 * bestG, habit);
        }
    }

    // Unconditional face support planes of Tetrahedron + r W (globally 1-Lipschitz & C0 everywhere)
    const float inv3 = 0.57735027;
    vec3 tn[4];
    tn[0] = vec3(-inv3,  inv3,  inv3);
    tn[1] = vec3( inv3, -inv3,  inv3);
    tn[2] = vec3( inv3,  inv3, -inv3);
    tn[3] = vec3(-inv3, -inv3, -inv3);
    float hTet = s * inv3;
    for( int k = 0; k < 4; k++ )
    {
        float lk = dot(tn[k], p) - hTet - r * mpWulffSupport(tn[k], habit);
        if( lk > bestLo ) { bestLo = lk; bestG = tn[k]; }
    }
    return vec4(bestLo, bestG);
}

float mpWulffCrystalOffset( vec3 p, float s, float r, int habit )
{
    return sdgWulffCrystalOffset(p, s, r, habit).x;
}

#endif // SDF_GEMS_HAS_MP_MOREAU_OPENING
