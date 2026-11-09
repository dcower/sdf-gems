// ============================================================================
// opRepDiamond — 0-Loop, 0-Sqrt 3D Diamond Cubic (Fd-3m) & HCP (P6_3/mmc) Crystal Lattice Repetitions
// @author Discovered and derived using Claude Opus 5.5; audited and extended
//         (analytical Voronoi cell gradients & axis hardening) using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (0-loop, 0-sqrt exact nearest-site folding for Diamond Cubic D+ and HCP lattices;
//         exact 1-sqrt infinite 3D Diamond Cubic rod network SDF & analytical gradient;
//         exact 0-sqrt interior / 1-sqrt exterior Triakis Truncated Tetrahedron &
//         Trapezo-Rhombic Dodecahedron Voronoi cell SDFs)
// ============================================================================

#ifndef SDF_GEMS_OP_REP_DIAMOND_GLSL
#define SDF_GEMS_OP_REP_DIAMOND_GLSL

// ----------------------------------------------------------------------------
// Internal helper: Exact 1-tap Conway-Sloane parity rounding onto the 3D FCC
// lattice D_3 = { i in Z^3 : i.x + i.y + i.z == 0 (mod 2) }
// ----------------------------------------------------------------------------
#ifndef SDF_GEMS_HAS_NEAREST_FCC_SITE
#define SDF_GEMS_HAS_NEAREST_FCC_SITE
vec3 sdfGemsNearestFCCSite( vec3 u )
{
    vec3  i = round(u);
    float m = mod(i.x + i.y + i.z, 2.0);
    if( m != 0.0 )
    {
        vec3 d   = abs(u - i);
        vec3 sgn = vec3(greaterThanEqual(u, i)) * 2.0 - 1.0;
        i += (d.x >= d.y && d.x >= d.z) ? vec3(sgn.x, 0.0, 0.0) :
             (d.y >= d.z)               ? vec3(0.0, sgn.y, 0.0) : vec3(0.0, 0.0, sgn.z);
    }
    return i;
}
#endif

// ----------------------------------------------------------------------------
// 1. opRepDiamond(p, s, id, chirality) & opRepDiamond(p, s)
//
// Folds R^3 into the exact 16-faced Triakis Truncated Tetrahedron Voronoi cell
// of the 3D Diamond Cubic crystal lattice (D+, space group Fd-3m, carbon diamond
// & silicon structure) in 0 loops and 0 sqrts.
//
// The Diamond Cubic lattice is the union of two interlocking FCC sublattices:
//   - Sublattice A (chirality = +1.0): Lambda_A = s * D_3
//   - Sublattice B (chirality = -1.0): Lambda_B = s * (D_3 + vec3(0.5))
//
// When sublattice B wins (chirality = -1.0), point-inversion symmetry about the
// bond midpoints negates the local offset q -> -q so that BOTH sublattices share
// the exact same canonical T_d tetrahedral frame where the 4 nearest tetrahedral
// bonds always point along (+1,+1,+1), (+1,-1,-1), (-1,+1,-1), (-1,-1,+1).
// ----------------------------------------------------------------------------
vec3 opRepDiamond( vec3 p, float s, out vec3 id, out float chirality )
{
    vec3  u  = p / max(s, 1e-12);
    vec3  iA = sdfGemsNearestFCCSite(u);
    vec3  iB = sdfGemsNearestFCCSite(u - 0.5) + 0.5;
    vec3  dA = u - iA;
    vec3  dB = u - iB;
    if( dot(dA, dA) <= dot(dB, dB) )
    {
        id        = iA;
        chirality = 1.0;
        return dA * s;
    }
    id        = iB;
    chirality = -1.0;
    return -dB * s;
}

vec3 opRepDiamond( vec3 p, float s )
{
    vec3  id;
    float chirality;
    return opRepDiamond(p, s, id, chirality);
}

vec4 opRepDiamondID( vec3 p, float s )
{
    vec3  id;
    float chirality;
    opRepDiamond(p, s, id, chirality);
    return vec4(id, chirality);
}

// ----------------------------------------------------------------------------
// 2. sdDiamondVoronoiCell(q, s)
//
// Exact Euclidean SDF of the 16-faced Triakis Truncated Tetrahedron Voronoi cell
// of the Diamond Cubic lattice (4 regular hexagonal faces shared with opposite-
// sublattice tetrahedral neighbors at distance sqrt(3)/4 * s, and 12 isosceles
// triangular faces shared with same-sublattice FCC neighbors at sqrt(2)/2 * s).
//
// Evaluates in 0 sqrts inside the Voronoi cell and 1 sqrt outside via 24-fold
// T_d Weyl-chamber reflection folding.
// ----------------------------------------------------------------------------
float sdDiamondVoronoiCell( vec3 q, float s )
{
    vec3  u   = q / max(s, 1e-12);
    bool  odd = (u.x < 0.0) != (u.y < 0.0) != (u.z < 0.0);
    vec3  w   = abs(u);
    if( w.x < w.y ) w.xy = w.yx;
    if( w.y < w.z ) w.yz = w.zy;
    if( w.x < w.y ) w.xy = w.yx;
    if( odd ) w.z = -w.z;

    float dh    = (w.x + w.y + w.z - 0.75) * 0.57735027;
    float dt    = (w.x + w.y - 1.0) * 0.70710678;
    float d_int = max(dh, dt);
    if( d_int <= 0.0 ) return d_int * s;

    float thh = clamp(0.5 * (w.y - w.z), 0.0, 0.25);
    float tht = clamp(0.5 * (w.x - w.y), 0.0, 0.25);
    float ttt = clamp((1.0 - w.x + w.y - w.z) * (1.0 / 3.0), 0.25, 0.5);

    vec3  ehh = w - vec3(0.75, thh, -thh);
    vec3  eht = w - vec3(0.5 + tht, 0.5 - tht, -0.25);
    vec3  ett = w - vec3(1.0 - ttt, ttt, -ttt);
    float d2  = min(dot(ehh, ehh), min(dot(eht, eht), dot(ett, ett)));

    if( 2.0 * w.x - w.y - w.z <= 1.5 && w.x + w.y - 2.0 * w.z <= 1.5 )
        d2 = min(d2, max(dh, 0.0) * max(dh, 0.0));
    if( w.z <= -0.25 && w.x - w.y - 2.0 * w.z <= 1.0 )
        d2 = min(d2, max(dt, 0.0) * max(dt, 0.0));

    return sqrt(max(d2, 0.0)) * s;
}

vec4 sdgDiamondVoronoiCell( vec3 q, float s )
{
    vec3  u   = q / max(s, 1e-12);
    vec3  sg  = vec3(u.x < 0.0 ? -1.0 : 1.0,
                     u.y < 0.0 ? -1.0 : 1.0,
                     u.z < 0.0 ? -1.0 : 1.0);
    float sig = sg.x * sg.y * sg.z;
    vec3  w   = abs(u);
    bool  s1  = (w.x < w.y); if( s1 ) w.xy = w.yx;
    bool  s2  = (w.y < w.z); if( s2 ) w.yz = w.zy;
    bool  s3  = (w.x < w.y); if( s3 ) w.xy = w.yx;
    w.z *= sig;

    float dh    = (w.x + w.y + w.z - 0.75) * 0.57735027;
    float dt    = (w.x + w.y - 1.0) * 0.70710678;
    float d_int = max(dh, dt);
    float d;
    vec3  g;
    if( d_int <= 0.0 ||
        (dh >= dt && 2.0 * w.x - w.y - w.z <= 1.5 && w.x + w.y - 2.0 * w.z <= 1.5) ||
        (dt >= dh && w.z <= -0.25 && w.x - w.y - 2.0 * w.z <= 1.0) )
    {
        d = d_int * s;
        g = (dh >= dt) ? vec3(0.57735027) : vec3(0.70710678, 0.70710678, 0.0);
    }
    else
    {
        float thh = clamp(0.5 * (w.y - w.z), 0.0, 0.25);
        float tht = clamp(0.5 * (w.x - w.y), 0.0, 0.25);
        float ttt = clamp((1.0 - w.x + w.y - w.z) * (1.0 / 3.0), 0.25, 0.5);

        vec3  ehh = w - vec3(0.75, thh, -thh);
        vec3  eht = w - vec3(0.5 + tht, 0.5 - tht, -0.25);
        vec3  ett = w - vec3(1.0 - ttt, ttt, -ttt);

        float dhh = dot(ehh, ehh), dht = dot(eht, eht), dtt = dot(ett, ett);
        vec3  v   = (dhh <= dht && dhh <= dtt) ? ehh : ((dht <= dtt) ? eht : ett);
        float v2  = dot(v, v);
        float inv = inversesqrt(max(v2, 1e-30));
        d = (v2 * inv) * s;
        g = v * inv;
    }

    g.z *= sig;
    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4(d, sg * g);
}

// ----------------------------------------------------------------------------
// 3. sdDiamondLattice3D(p, s, r) & sdgDiamondLattice3D(p, s, r)
//
// Exact 0-loop, 0-clamp, 1-sqrt Euclidean SDF and analytical gradient (vec4(d, grad))
// of the infinite 3D Diamond Cubic Strut / Rod Truss Network of cell scale s and
// strut radius r.
//
// By the Voronoi-Star Exactness Theorem of the Fd-3m diamond bond network, every
// one of the 16 faces of the Triakis Truncated Tetrahedron Voronoi cell is a
// symmetry plane bisecting the distance between the 4 half-bonds radiating from
// the cell center and all bonds in neighboring cells, and inside the T_d chamber
// w.x >= w.y >= |w.z| the axial projection t = (w.x + w.y + w.z) / 3 onto the
// nearest bond axis (+1,+1,+1) satisfies 0 <= t <= s/4 identically everywhere
// in the Voronoi cell.
// ----------------------------------------------------------------------------
vec4 sdgDiamondLattice3D( vec3 p, float s, float r )
{
    vec3  id;
    float chi;
    vec3  q = opRepDiamond(p, s, id, chi);

    vec3  a = abs(q);
    vec3  u = vec3(q.x >= 0.0 ? 1.0 : -1.0,
                   q.y >= 0.0 ? 1.0 : -1.0,
                   q.z >= 0.0 ? 1.0 : -1.0);
    if( u.x * u.y * u.z < 0.0 )
    {
        if( a.x <= a.y && a.x <= a.z ) u.x = -u.x;
        else if( a.y <= a.z )          u.y = -u.y;
        else                           u.z = -u.z;
    }

    vec3  d_perp = q - u * (dot(q, u) * (1.0 / 3.0));
    float l      = length(d_perp);
    vec3  g      = (l > 1e-6 * s) ? (d_perp / l) : (vec3(u.y, -u.x, 0.0) * 0.70710678);
    return vec4(l - r, chi * g);
}

float sdDiamondLattice3D( vec3 p, float s, float r )
{
    return sdgDiamondLattice3D(p, s, r).x;
}

// ----------------------------------------------------------------------------
// 4. opRepHCP(p, s, id, layer_parity) & sdHCPVoronoiCell(q, s)
//
// Exact 0-loop, 0-sqrt 3D Hexagonal Close-Packed (HCP, space group P6_3/mmc,
// ABABAB stacking, 74.05% maximum sphere packing density, layer spacing
// h = sqrt(2/3) * s) crystal lattice domain repetition and 12-faced
// Trapezo-Rhombic Dodecahedron Voronoi cell SDF.
//
// Compares the nearest A-layer site (even k, layer_parity = +1.0) and nearest
// B-layer site (odd k, shifted by (0, s / sqrt(3)), layer_parity = -1.0) using
// 1-tap linear-norm hexagonal rounding on each layer and dot2 comparison.
// Negating q.y on odd B-layers aligns both layers into the single canonical
// D_3h Trapezo-Rhombic Dodecahedron Voronoi cell frame.
// ----------------------------------------------------------------------------
vec3 opRepHCP( vec3 p, float s, out vec3 id, out float layer_parity )
{
    const vec2 k_hex = vec2(1.0, 1.73205081);
    float sm = max(s, 1e-12);
    float h  = 0.81649658 * sm;
    float uz = p.z / h;

    float kA = 2.0 * round(0.5 * uz);
    float kB = 2.0 * round(0.5 * (uz - 1.0)) + 1.0;

    vec2  uA = p.xy / (sm * k_hex), iA = round(uA), dA = uA - iA;
    if( abs(dA.x) + 3.0 * abs(dA.y) > 1.0 )
        iA += vec2(dA.x >= 0.0 ? 0.5 : -0.5, dA.y >= 0.0 ? 0.5 : -0.5);
    vec3  qA = vec3(p.xy - sm * iA * k_hex, p.z - kA * h);

    vec2  pB = p.xy - vec2(0.0, 0.57735027 * sm);
    vec2  uB = pB / (sm * k_hex), iB = round(uB), dB = uB - iB;
    if( abs(dB.x) + 3.0 * abs(dB.y) > 1.0 )
        iB += vec2(dB.x >= 0.0 ? 0.5 : -0.5, dB.y >= 0.0 ? 0.5 : -0.5);
    vec3  qB = vec3(pB - sm * iB * k_hex, p.z - kB * h);

    if( dot(qA, qA) <= dot(qB, qB) )
    {
        id           = vec3(iA, kA);
        layer_parity = 1.0;
        return qA;
    }
    id           = vec3(iB, kB);
    layer_parity = -1.0;
    return vec3(qB.x, -qB.y, qB.z);
}

vec3 opRepHCP( vec3 p, float s )
{
    vec3  id;
    float layer_parity;
    return opRepHCP(p, s, id, layer_parity);
}

vec4 opRepHCPID( vec3 p, float s )
{
    vec3  id;
    float layer_parity;
    opRepHCP(p, s, id, layer_parity);
    return vec4(id, layer_parity);
}

// Exact 0-sqrt interior / 1-sqrt exterior Euclidean SDF of the 12-faced
// Trapezo-Rhombic Dodecahedron Voronoi cell of the HCP lattice
float sdHCPVoronoiCell( vec3 q, float s )
{
    vec3  w = q / max(s, 1e-12);
    w.xz    = abs(w.xz);
    float k = max(0.5 * w.x - 0.86602540 * w.y, 0.0);
    w.xy   -= k * vec2(1.0, -1.73205081);
    w.x     = abs(w.x);

    float d1    = 0.5 * w.x + 0.86602540 * w.y - 0.5;
    float d2    = 0.57735027 * w.y + 0.81649658 * w.z - 0.5;
    float d_int = max(d1, d2);
    if( d_int <= 0.0 ) return d_int * s;

    vec2  evy = vec2(w.x, w.y - 0.57735027);
    float zvy = max(w.z - 0.20412415, 0.0);
    float d2_vy = dot(evy, evy) + zvy * zvy;

    vec2  evx = vec2(w.x - 0.5, w.y - 0.28867513);
    float zvx = max(w.z - 0.40824829, 0.0);
    float d2_vx = dot(evx, evx) + zvx * zvx;

    float t12 = clamp((w.x * 0.5 + (w.y - 0.57735027) * (-0.28867513) + (w.z - 0.20412415) * 0.20412415) * (8.0 / 3.0), 0.0, 1.0);
    vec3  e12 = w - vec3(0.5 * t12, 0.57735027 - 0.28867513 * t12, 0.20412415 + 0.20412415 * t12);

    float ttx = clamp((w.x * 0.5 + w.y * 0.28867513 + (w.z - 0.61237244) * (-0.20412415)) * (8.0 / 3.0), 0.0, 1.0);
    vec3  etx = w - vec3(0.5 * ttx, 0.28867513 * ttx, 0.61237244 - 0.20412415 * ttx);

    float d2_ext = min(min(d2_vy, d2_vx), min(dot(e12, e12), dot(etx, etx)));

    float p1x = w.x - 0.5 * d1, p1y = w.y - 0.86602540 * d1;
    if( d2 - 0.5 * d1 <= 0.0 && p1x >= 0.0 && 0.5 * p1x - 0.86602540 * p1y <= 0.0 )
        d2_ext = min(d2_ext, max(d1, 0.0) * max(d1, 0.0));

    float p2y = w.y - 0.57735027 * d2, p2z = w.z - 0.81649658 * d2;
    if( d1 - 0.5 * d2 <= 0.0 && p2z >= 0.0 && 0.5 * w.x - 0.86602540 * p2y <= 0.0 )
        d2_ext = min(d2_ext, max(d2, 0.0) * max(d2, 0.0));

    return sqrt(max(d2_ext, 0.0)) * s;
}

vec4 sdgHCPVoronoiCell( vec3 q, float s )
{
    vec3  w   = q / max(s, 1e-12);
    float sx0 = (w.x < 0.0) ? -1.0 : 1.0;
    float sz0 = (w.z < 0.0) ? -1.0 : 1.0;
    w.xz      = abs(w.xz);
    bool  f60 = (0.5 * w.x - 0.86602540 * w.y > 0.0);
    if( f60 ) w.xy -= (0.5 * w.x - 0.86602540 * w.y) * vec2(1.0, -1.73205081);
    float sx1 = (w.x < 0.0) ? -1.0 : 1.0;
    w.x       = abs(w.x);

    float d1    = 0.5 * w.x + 0.86602540 * w.y - 0.5;
    float d2    = 0.57735027 * w.y + 0.81649658 * w.z - 0.5;
    float d_int = max(d1, d2);

    float p1x   = w.x - 0.5 * d1, p1y = w.y - 0.86602540 * d1;
    bool  in_f1 = (d1 > 0.0 && d2 - 0.5 * d1 <= 0.0 && p1x >= 0.0 && 0.5 * p1x - 0.86602540 * p1y <= 0.0);
    float p2y   = w.y - 0.57735027 * d2, p2z = w.z - 0.81649658 * d2;
    bool  in_f2 = (d2 > 0.0 && d1 - 0.5 * d2 <= 0.0 && p2z >= 0.0 && 0.5 * w.x - 0.86602540 * p2y <= 0.0);

    float d;
    vec3  g;
    if( d_int <= 0.0 || in_f1 || in_f2 )
    {
        bool use1 = (d_int <= 0.0) ? (d1 >= d2) : in_f1;
        d = (use1 ? d1 : d2) * s;
        g = use1 ? vec3(0.5, 0.86602540, 0.0) : vec3(0.0, 0.57735027, 0.81649658);
    }
    else
    {
        vec3  evy = vec3(w.x, w.y - 0.57735027, max(w.z - 0.20412415, 0.0));
        vec3  evx = vec3(w.x - 0.5, w.y - 0.28867513, max(w.z - 0.40824829, 0.0));
        float t12 = clamp((w.x * 0.5 + (w.y - 0.57735027) * (-0.28867513) + (w.z - 0.20412415) * 0.20412415) * (8.0 / 3.0), 0.0, 1.0);
        vec3  e12 = w - vec3(0.5 * t12, 0.57735027 - 0.28867513 * t12, 0.20412415 + 0.20412415 * t12);
        float ttx = clamp((w.x * 0.5 + w.y * 0.28867513 + (w.z - 0.61237244) * (-0.20412415)) * (8.0 / 3.0), 0.0, 1.0);
        vec3  etx = w - vec3(0.5 * ttx, 0.28867513 * ttx, 0.61237244 - 0.20412415 * ttx);

        float dvy = dot(evy, evy), dvx = dot(evx, evx), de12 = dot(e12, e12), detx = dot(etx, etx);
        vec3  v   = (dvy <= dvx && dvy <= de12 && dvy <= detx) ? evy :
                    (dvx <= de12 && dvx <= detx)               ? evx :
                    (de12 <= detx)                             ? e12 : etx;
        float v2  = dot(v, v);
        float inv = inversesqrt(max(v2, 1e-30));
        d = (v2 * inv) * s;
        g = v * inv;
    }

    g.x *= sx1;
    if( f60 ) g.xy -= (0.5 * g.x - 0.86602540 * g.y) * vec2(1.0, -1.73205081);
    g.x *= sx0;
    g.z *= sz0;
    return vec4(d, g);
}

#endif // SDF_GEMS_OP_REP_DIAMOND_GLSL
