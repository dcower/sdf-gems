// ============================================================================
// opRepQuasicrystal — Cut-and-Project Domain Repetition (Icosahedral Quasicrystals)
//
//   opRepQuasi6, opQuasi6Project, hashQuasi6,
//   sdQuasicrystalBalls        (slice of the 6D ball field, 1-Lipschitz bound)
//   sdQuasicrystalBallsTight   (7 candidate atoms + certificate, conservative bound)
//   sdQuasicrystalBallsCorners (64 candidate atoms + certificate)
//   sdQuasicrystalCells        (slice of the 6D box: eroded cells of the icosahedral 6-grid)
//   sdQuasicrystalRods         (slice of the 6D axis-line field: prolate ellipsoids)
//   opRepAmmannBeenker2D, sdAmmannBeenkerDiscs2D   (Z^4 -> R^2, 8-fold)
//   opRepPenrose2D, sdPenroseDiscs2D               (Z^5 -> R^2, 10-fold, pentagrid ids)
//
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  sdQuasicrystalBalls: conservative 1-Lipschitz bound with the exact zero set and sign of the
//         union of balls (centres E^T n, radii sqrt(r^2 - h_n^2), over n in Z^6 with h_n < r).
//         sdQuasicrystalBallsTight: conservative bound, continuous (a candidate enters / leaves the set
//         exactly when its 6D distance equals the certificate radius), exact on 73% of exterior points at r = 1/2.
//         sdQuasicrystalBallsCorners: conservative bound for every r; for r >= 0.5 it matched the brute-force
//         exterior SDF of the union on every sampled point (measured, no proof; 99.998% at r = 0.49).
//         sdQuasicrystalCells: exact inside each (convex) cell and near faces; conservative bound outside.
//         sdQuasicrystalRods, *Discs2D: conservative 1-Lipschitz bounds.
//
// Idea: p in R^3 is embedded isometrically into R^6 as P = E p, where the six rows of E are the six
// 5-fold axes of the icosahedron scaled by 1/sqrt(2) (so E^T E = I_3). The residual Q = P - round(P)
// to the nearest point of Z^6 is a 6D "cell coordinate". Any 1-Lipschitz function of Q is 1-Lipschitz
// in p, and its zero set is an aperiodic, icosahedrally symmetric structure with no cell seams.
// All distances are returned in world units; `scale` is the 6D lattice constant, `r` is in lattice units.
// ============================================================================

#ifndef SDF_GEMS_OP_REP_QUASICRYSTAL_GLSL
#define SDF_GEMS_OP_REP_QUASICRYSTAL_GLSL

// k = 1 / sqrt(2 (1 + tau^2)),  k * tau;  k^2 + (k tau)^2 = 1/2.
#define SDF_GEMS_QC_K1 0.37174803446018456
#define SDF_GEMS_QC_KT 0.60150095500754570

// Rows of E (the six 5-fold axes / sqrt(2)); row i is the parallel projection of the 6D basis vector e_i.
const vec3 SDF_GEMS_QC_AXIS[6] = vec3[6](
    vec3( 0.0,             SDF_GEMS_QC_K1, SDF_GEMS_QC_KT),
    vec3( 0.0,            -SDF_GEMS_QC_K1, SDF_GEMS_QC_KT),
    vec3( SDF_GEMS_QC_K1,  SDF_GEMS_QC_KT, 0.0),
    vec3(-SDF_GEMS_QC_K1,  SDF_GEMS_QC_KT, 0.0),
    vec3( SDF_GEMS_QC_KT,  0.0,            SDF_GEMS_QC_K1),
    vec3( SDF_GEMS_QC_KT,  0.0,           -SDF_GEMS_QC_K1));

// 1. Parallel projection E^T (R^6 -> R^3) of a 6D vector stored as two vec3 (components 0..2 and 3..5).
//    E^T applied to an integer cell id gives the atom / tiling-vertex position (in lattice units).
vec3 opQuasi6Project( in vec3 qa, in vec3 qb )
{
    return vec3(SDF_GEMS_QC_K1 * (qa.z - qb.x) + SDF_GEMS_QC_KT * (qb.y + qb.z),
                SDF_GEMS_QC_K1 * (qa.x - qa.y) + SDF_GEMS_QC_KT * (qa.z + qb.x),
                SDF_GEMS_QC_KT * (qa.x + qa.y) + SDF_GEMS_QC_K1 * (qb.y - qb.z));
}

// 2. Cut-and-project repetition: residual Q = E p - round(E p) in [-1/2, 1/2]^6 and the integer cell id.
//    p is in lattice units (divide by the lattice constant first). |Q| is the distance from E p to Z^6.
void opRepQuasi6( in vec3 p, out vec3 qa, out vec3 qb, out ivec3 na, out ivec3 nb )
{
    vec3 a = vec3( SDF_GEMS_QC_K1 * p.y + SDF_GEMS_QC_KT * p.z,
                  -SDF_GEMS_QC_K1 * p.y + SDF_GEMS_QC_KT * p.z,
                   SDF_GEMS_QC_K1 * p.x + SDF_GEMS_QC_KT * p.y);
    vec3 b = vec3(-SDF_GEMS_QC_K1 * p.x + SDF_GEMS_QC_KT * p.y,
                   SDF_GEMS_QC_KT * p.x + SDF_GEMS_QC_K1 * p.z,
                   SDF_GEMS_QC_KT * p.x - SDF_GEMS_QC_K1 * p.z);
    vec3 ra = floor(a + 0.5);
    vec3 rb = floor(b + 0.5);
    qa = a - ra;
    qb = b - rb;
    na = ivec3(ra);
    nb = ivec3(rb);
}

void opRepQuasi6( in vec3 p, out vec3 qa, out vec3 qb )
{
    ivec3 na, nb;
    opRepQuasi6(p, qa, qb, na, nb);
}

// 3. Hash of a 6D integer cell id to [0, 1) for per-atom / per-cell colouring.
float hashQuasi6( in ivec3 na, in ivec3 nb )
{
    uvec3 a = uvec3(na), b = uvec3(nb);
    uint h = a.x * 0x9E3779B1u + a.y * 0x85EBCA77u + a.z * 0xC2B2AE3Du
           + b.x * 0x27D4EB2Fu + b.y * 0x165667B1u + b.z * 0xD3A2646Cu;
    h ^= h >> 15; h *= 0x2C1B3C6Du; h ^= h >> 12; h *= 0x297A2D39u; h ^= h >> 15;
    return float(h >> 8) * (1.0 / 16777216.0);
}

// 4. Slice of the 6D ball field: d = scale * (dist(E p / scale, Z^6) - r).
//    Zero set and sign are exactly those of the union of balls centred at scale * E^T n with radius
//    scale * sqrt(r^2 - h_n^2), over all n in Z^6 whose distance h_n to the slice is < r.
//    1-Lipschitz, |d| <= true distance. The balls are pairwise disjoint iff r <= 1/2.
//    On the surface the exact normal is normalize(opQuasi6Project(qa, qb)).
float sdQuasicrystalBalls( in vec3 p, in float scale, in float r )
{
    vec3 qa, qb;
    opRepQuasi6(p / scale, qa, qb);
    return scale * (sqrt(dot(qa, qa) + dot(qb, qb)) - r);
}

// 5. Tighter conservative bound: exact distances to the balls of the 7 lattice points N and N + s_i e_i
//    (s_i = sign(Q_i)); every other lattice point is at 6D distance >= R, hence its ball is at distance
//    >= R - r, with R^2 = |Q|^2 + min(2 - 2 (a1 + a2), 1 + 2 a6), a1 >= a2 >= ... >= a6 the sorted |Q_i|.
//    Never overestimates the exterior distance. Inside a ball it returns min over candidates of
//    (d_n - rho_n), which is the exact interior distance when the balls are disjoint (r <= 1/2).
float sdQuasicrystalBallsTight( in vec3 p, in float scale, in float r )
{
    vec3 qa, qb;
    opRepQuasi6(p / scale, qa, qb);
    float q2 = dot(qa, qa) + dot(qb, qb);
    float r2 = r * r;
    float f  = sqrt(q2) - r;
    vec3  qp = opQuasi6Project(qa, qb);
    float d2 = dot(qp, qp);
    float h2 = max(q2 - d2, 0.0);
    float c  = (h2 < r2) ? sqrt(d2) - sqrt(r2 - h2) : 1e20;
    float m1 = 0.0, m2 = 0.0, mn = 1.0;
    for( int i = 0; i < 6; i++ )
    {
        float qi = (i < 3) ? qa[i] : qb[i - 3];
        float ai = abs(qi);
        float si = (qi < 0.0) ? -1.0 : 1.0;
        m2 = max(m2, min(m1, ai));
        m1 = max(m1, ai);
        mn = min(mn, ai);
        vec3  cp  = qp - si * SDF_GEMS_QC_AXIS[i];
        float cd2 = dot(cp, cp);
        float ch2 = max(q2 + 1.0 - 2.0 * ai - cd2, 0.0);
        if( ch2 < r2 ) c = min(c, sqrt(cd2) - sqrt(r2 - ch2));
    }
    if( f < 0.0 ) return scale * ((c < 1e19) ? c : f);
    float lb = max(sqrt(q2 + min(2.0 - 2.0 * (m1 + m2), 1.0 + 2.0 * mn)) - r, f);
    return scale * min(c, lb);
}

// 6. 64 candidates: all corners N + (s_i b_i) of the unit 6-cube that contains E p, pruned with
//    |d_n - rho_n| >= |P - n| - r. Lattice points outside the cube are at 6D distance
//    >= R = sqrt(|Q|^2 + 1 + 2 min|Q_i|), which certifies the result as a conservative bound for every r.
//    Measured: equal to the brute-force exterior SDF of the union of balls for r >= 1/2 (no proof);
//    for r < 1/2 the nearest atom can lie outside the cube and the certificate term takes over.
float sdQuasicrystalBallsCorners( in vec3 p, in float scale, in float r )
{
    vec3 qa, qb;
    opRepQuasi6(p / scale, qa, qb);
    float q2 = dot(qa, qa) + dot(qb, qb);
    float r2 = r * r;
    float f  = sqrt(q2) - r;
    vec3  sa = vec3(qa.x < 0.0 ? -1.0 : 1.0, qa.y < 0.0 ? -1.0 : 1.0, qa.z < 0.0 ? -1.0 : 1.0);
    vec3  sb = vec3(qb.x < 0.0 ? -1.0 : 1.0, qb.y < 0.0 ? -1.0 : 1.0, qb.z < 0.0 ? -1.0 : 1.0);
    vec3  wa = 1.0 - 2.0 * abs(qa);
    vec3  wb = 1.0 - 2.0 * abs(qb);
    float c  = 1e20;
    // Candidates farther than sqrt(lim) in 6D cannot win: outside, d_n - rho_n >= |P - n| - r >= c;
    // inside, only balls that contain p matter and those have |P - n| < r. The nearest point N (m = 0)
    // is always evaluated, so that rounding of f on the surface cannot prune every candidate.
    float lim = 1e30;
    for( int m = 0; m < 64; m++ )
    {
        vec3  ba = vec3(float(m & 1), float((m >> 1) & 1), float((m >> 2) & 1));
        vec3  bb = vec3(float((m >> 3) & 1), float((m >> 4) & 1), float((m >> 5) & 1));
        float n2 = q2 + dot(ba, wa) + dot(bb, wb);          // |P - n|^2
        if( n2 >= lim ) continue;
        vec3  cp  = opQuasi6Project(qa - sa * ba, qb - sb * bb);
        float cd2 = dot(cp, cp);
        float ch2 = max(n2 - cd2, 0.0);
        if( ch2 < r2 )
        {
            c = min(c, sqrt(cd2) - sqrt(r2 - ch2));
            float t = max(c, 0.0) + r;
            if( f >= 0.0 ) lim = t * t;
        }
        if( f < 0.0 ) lim = r2;
    }
    if( f < 0.0 ) return scale * ((c < 1e19) ? c : f);
    vec3  wm = max(wa, wb);
    float mx = max(max(wm.x, wm.y), wm.z);                 // 1 - 2 min|Q_i|
    float lb = max(sqrt(q2 + 2.0 - mx) - r, f);            // R^2 = |Q|^2 + 1 + 2 min|Q_i|
    return scale * min(c, lb);
}

// 7. Slice of the 6D box field, rescaled by sqrt(2): d = scale * sqrt(2) * (max_i |Q_i| - r), 0 < r < 1/2.
//    Each cell is the convex polytope { |E_i . p - n_i| < r for all i } bounded by planes normal to the
//    5-fold axes (a regular dodecahedron of inradius sqrt(2) r scale when n lies on the slice). The cells
//    are the meshes of the icosahedral 6-grid eroded by sqrt(2) (1/2 - r); the negated function is the exact
//    exterior distance to the union of the 6 families of plates of thickness sqrt(2) (1 - 2 r).
float sdQuasicrystalCells( in vec3 p, in float scale, in float r )
{
    vec3 qa, qb;
    opRepQuasi6(p / scale, qa, qb);
    vec3 m = max(abs(qa), abs(qb));
    return scale * 1.4142135623730951 * (max(max(m.x, m.y), m.z) - r);
}

// 8. Slice of the 6D field of axis-parallel lattice lines: d = scale * (sqrt(|Q|^2 - max_i Q_i^2) - r).
//    Each 6D line n + t e_i within distance r of the slice appears as a prolate ellipsoid (aspect sqrt(2))
//    whose long axis is the 5-fold axis i. 1-Lipschitz conservative bound.
float sdQuasicrystalRods( in vec3 p, in float scale, in float r )
{
    vec3 qa, qb;
    opRepQuasi6(p / scale, qa, qb);
    vec3 a2 = qa * qa, b2 = qb * qb;
    vec3 m  = max(a2, b2);
    float s = a2.x + a2.y + a2.z + b2.x + b2.y + b2.z - max(max(m.x, m.y), m.z);
    return scale * (sqrt(max(s, 0.0)) - r);
}

// 9. 2D, 8-fold (Ammann-Beenker type): Z^4 -> R^2 with rows (cos(k pi/4), sin(k pi/4)) / sqrt(2), k = 0..3.
//    Returns the residual Q in [-1/2, 1/2]^4 and the integer cell id. p in lattice units.
vec4 opRepAmmannBeenker2D( in vec2 p, out ivec4 n )
{
    vec4 P = vec4(0.7071067811865476 * p.x, 0.5 * (p.x + p.y), 0.7071067811865476 * p.y, 0.5 * (p.y - p.x));
    vec4 R = floor(P + 0.5);
    n = ivec4(R);
    return P - R;
}

// Discs at the projections of the lattice points within r of the slice (disc window), radii sqrt(r^2 - h^2).
float sdAmmannBeenkerDiscs2D( in vec2 p, in float scale, in float r )
{
    ivec4 n;
    return scale * (length(opRepAmmannBeenker2D(p / scale, n)) - r);
}

// 10. 2D, 10-fold (Penrose type): Z^5 -> R^2 with rows sqrt(2/5) (cos(2 pi k/5), sin(2 pi k/5)), k = 0..4,
//     and a common phase gamma added to the five coordinates. n[] is constant on each mesh of the pentagrid
//     { P_k + gamma = integer + 1/2 }; with 5 gamma = 1/2 (mod 1), e.g. gamma = 0.1, the points
//     sum_k n[k] e_k are the vertices of a Penrose rhomb tiling (de Bruijn 1981).
//     (1/2 - max_k |q[k]|) * sqrt(5/2) is the exact distance to the nearest pentagrid line.
void opRepPenrose2D( in vec2 p, in float gamma, out float q[5], out int n[5] )
{
    const vec2 e[5] = vec2[5](
        vec2( 0.6324555320336759,  0.0),
        vec2( 0.1954395076113904,  0.6015009550075457),
        vec2(-0.5116672736282283,  0.3717480344601846),
        vec2(-0.5116672736282283, -0.3717480344601846),
        vec2( 0.1954395076113904, -0.6015009550075457));
    for( int k = 0; k < 5; k++ )
    {
        float P = dot(p, e[k]) + gamma;
        float R = floor(P + 0.5);
        q[k] = P - R;
        n[k] = int(R);
    }
}

// Discs from the 5D ball field (gamma = 0): a 10-fold cut-and-project set with a ball window.
// NOT the Penrose vertex set (that needs the pentagrid ids above).
float sdPenroseDiscs2D( in vec2 p, in float scale, in float r )
{
    float q[5]; int n[5];
    opRepPenrose2D(p / scale, 0.0, q, n);
    return scale * (sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3] + q[4] * q[4]) - r);
}

#endif // SDF_GEMS_OP_REP_QUASICRYSTAL_GLSL
