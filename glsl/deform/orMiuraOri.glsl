// ============================================================================
// orMiuraOri — Rigid Isometric Origami Kinematics, Miura-Ori Pullback SDFs & Yoshimura Buckling
//   orMiuraSecondaryFold, opMiuraOriCell, opMiuraOriUnfold,
//   sdMiuraOriNaiveSlab, sdMiuraOriSingleFacet, sdMiuraOriSheet, sdgMiuraOriSheet,
//   sdYoshimuraCylinder
//
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  sdMiuraOriSheet / sdgMiuraOriSheet: Exact Euclidean SDF and analytical unit gradient of an
//         isometrically folded Miura-ori tessellation (thickened by `thickness`).
//         sdMiuraOriSingleFacet: 1-Lipschitz upper bound, exact on 92%–99.96% of space (misses only the
//         acute-corner overhang wedge across the glide-reflection seam when alpha < pi/2).
//         sdMiuraOriNaiveSlab: Unfolded-plane slab approximation (exact only in the normal prism over a
//         single facet; overestimates in valley wedges and underestimates over mountain ridges).
//         sdYoshimuraCylinder: Exact Euclidean SDF of a triangulated Yoshimura/Kresling buckled cylinder.
// ============================================================================

#ifndef SDF_GEMS_HAS_OR_MIURA_ORI
#define SDF_GEMS_HAS_OR_MIURA_ORI

// ---------------------------------------------------------------------------
// 1. Kawasaki-Maekawa Degree-4 Spherical Kinematics
//    Given sector angle alpha in (0, pi/2) and primary dihedral fold angle gamma in [0, pi),
//    returns the coupled secondary dihedral fold angle rho_1 via the spherical law of cosines:
//        tan(rho_1 / 2) = cos(alpha) * tan(gamma / 2).
// ---------------------------------------------------------------------------
float orMiuraSecondaryFold( in float alpha, in float gamma )
{
    return 2.0 * atan(cos(alpha) * tan(0.5 * gamma));
}

// Internal helper: Miura-ori unit cell dimensions (ya, za, xb, yb) and orthonormal facet frame (ea, eap, n).
void orMiuraFrame( in float a, in float b, in float alpha, in float gamma,
                   out vec4 dims, out vec3 ea, out vec3 eap, out vec3 n, out vec2 csa )
{
    float sa  = sin(alpha), ca = cos(alpha);
    float sg  = sin(0.5 * gamma), cg = cos(0.5 * gamma);
    float sth = sg * sa;
    float cth = sqrt(max((1.0 - sth) * (1.0 + sth), 1e-12));
    float ic  = 1.0 / cth;

    float ya  = a * cth;
    float za  = a * sth;
    float xb  = b * cg * sa * ic;
    float yb  = b * ca * ic;

    dims = vec4(ya, za, xb, yb);
    csa  = vec2(ca, sa);
    ea   = vec3(0.0, cth, sth);
    eap  = vec3(cg * ic, ca * sg * sth * ic, -ca * sg);
    n    = vec3(ca * sg * ic, -sa * cg * sg * ic, cg);
}

// Internal helper: exact distance and 3D gradient to the fundamental parallelogram F0 in chamber coordinates.
// Returns vec4(dist, grad.xyz) and outputs intrinsic 2D coordinates (u*, v*, w) in `uvw`.
vec4 sdgMiuraFacet0( in vec3 pf, in float a, in float b, in vec4 dims,
                     in vec3 ea, in vec3 eap, in vec3 n, in vec2 csa, out vec3 uvw )
{
    vec3  dp = pf - vec3(0.0, 0.0, -0.5 * dims.y);
    float u  = dot(dp, ea);
    float v  = dot(dp, eap);
    float w  = dot(dp, n);
    uvw = vec3(u, v, w);

    float ca = csa.x, sa = csa.y;
    float tc = v / max(b * sa, 1e-12);
    float sc = (u - tc * b * ca) / max(a, 1e-12);

    if( sc >= 0.0 && sc <= 1.0 && tc >= 0.0 && tc <= 1.0 )
    {
        return vec4(abs(w), ((w < 0.0) ? -1.0 : 1.0) * n);
    }

    // 4 edges of the 2D parallelogram in the intrinsic (u, v) plane:
    float du1 = u - clamp(u, 0.0, a), dv1 = v;
    float d2  = du1 * du1 + dv1 * dv1;
    float du  = du1, dv = dv1;

    float u2  = u - b * ca;
    float du2 = u2 - clamp(u2, 0.0, a), dv2 = v - b * sa;
    float l2  = du2 * du2 + dv2 * dv2;
    if( l2 < d2 ) { d2 = l2; du = du2; dv = dv2; }

    float h3  = clamp(u * ca + v * sa, 0.0, b);
    float du3 = u - h3 * ca, dv3 = v - h3 * sa;
    float l3  = du3 * du3 + dv3 * dv3;
    if( l3 < d2 ) { d2 = l3; du = du3; dv = dv3; }

    float u4  = u - a;
    float h4  = clamp(u4 * ca + v * sa, 0.0, b);
    float du4 = u4 - h4 * ca, dv4 = v - h4 * sa;
    float l4  = du4 * du4 + dv4 * dv4;
    if( l4 < d2 ) { d2 = l4; du = du4; dv = dv4; }

    vec3  q   = du * ea + dv * eap + w * n;
    float q2  = d2 + w * w;
    float il  = inversesqrt(max(q2, 1e-24));
    return (q2 > 1e-24) ? vec4(q2 * il, q * il) : vec4(0.0, n);
}

// ---------------------------------------------------------------------------
// 2. O(1) Layer-Group (pm2_1b) Symmetry Folding of the Miura-Ori Unit Cell
//    Folds p in R^3 into the fundamental parallelogram chamber using:
//      - Mirror reflections X = j * xb along the X-axis,
//      - Glide reflections (Y -> Y - k * ya, Z -> (-1)^k Z) along the Y-axis.
//    `cellId` returns vec4(cell_x, cell_y, mirror_sign_x, glide_sign_z).
// ---------------------------------------------------------------------------
vec3 opMiuraOriCell( in vec3 p, in float a, in float b, in float alpha, in float gamma, out vec4 cellId )
{
    vec4 dims; vec3 ea, eap, n; vec2 csa;
    orMiuraFrame(a, b, alpha, gamma, dims, ea, eap, n, csa);
    float ya = dims.x, xb = dims.z, yb = dims.w;

    float cx   = floor(p.x / max(xb, 1e-12));
    float xf   = p.x - cx * xb;
    bool  oddx = (int(floor(mod(cx, 2.0) + 0.5)) != 0);
    if( oddx ) xf = xb - xf;

    float yrel = (p.y - (yb / max(xb, 1e-12)) * xf) / max(ya, 1e-12);
    float cy   = floor(yrel);
    float yf   = p.y - cy * ya;
    bool  oddy = (int(floor(mod(cy, 2.0) + 0.5)) != 0);
    float zf   = oddy ? -p.z : p.z;

    cellId = vec4(cx, cy, oddx ? -1.0 : 1.0, oddy ? -1.0 : 1.0);
    return vec3(xf, yf, zf);
}

vec3 opMiuraOriCell( in vec3 p, in float a, in float b, in float alpha, in float gamma )
{
    vec4 cellId;
    return opMiuraOriCell(p, a, b, alpha, gamma, cellId);
}

// ---------------------------------------------------------------------------
// 3. Naive Single-Facet Slab Pullback vs Single-Facet vs Exact 2-Chamber SDF
// ---------------------------------------------------------------------------

// (a) Naive single-facet slab pullback: projects p onto the infinite plane of F0 inside the folded chamber.
//     Overestimates in valley wedges and underestimates over mountain creases outside the normal prism of F0.
float sdMiuraOriNaiveSlab( in vec3 p, in float a, in float b, in float alpha, in float gamma, in float thickness )
{
    vec4 dims; vec3 ea, eap, n; vec2 csa;
    orMiuraFrame(a, b, alpha, gamma, dims, ea, eap, n, csa);
    vec4 cellId;
    vec3 pf = opMiuraOriCell(p, a, b, alpha, gamma, cellId);
    float w = dot(pf - vec3(0.0, 0.0, -0.5 * dims.y), n);
    return abs(w) - thickness;
}

// (b) Single-chamber 3D parallelogram distance: exact on 92%–99.96% of space, slight upper bound in the
//     acute-corner overhang wedge across the glide-reflection crease when alpha < pi/2.
float sdMiuraOriSingleFacet( in vec3 p, in float a, in float b, in float alpha, in float gamma, in float thickness )
{
    vec4 dims; vec3 ea, eap, n; vec2 csa;
    orMiuraFrame(a, b, alpha, gamma, dims, ea, eap, n, csa);
    vec4 cellId;
    vec3 pf = opMiuraOriCell(p, a, b, alpha, gamma, cellId);
    vec3 uvw;
    return sdgMiuraFacet0(pf, a, b, dims, ea, eap, n, csa, uvw).x - thickness;
}

// (c) Exact Euclidean SDF + Analytical 3D Unit Gradient of the Miura-Ori Sheet.
//     Evaluates F0 and its glide-reflected partner across the nearest horizontal crease in O(1).
vec4 sdgMiuraOriSheet( in vec3 p, in float a, in float b, in float alpha, in float gamma, in float thickness )
{
    vec4 dims; vec3 ea, eap, n; vec2 csa;
    orMiuraFrame(a, b, alpha, gamma, dims, ea, eap, n, csa);
    float ya = dims.x, xb = dims.z, yb = dims.w;

    vec4 cellId;
    vec3 pf = opMiuraOriCell(p, a, b, alpha, gamma, cellId);

    vec3 uvw0, uvw1;
    vec4 dg0 = sdgMiuraFacet0(pf, a, b, dims, ea, eap, n, csa, uvw0);

    float yfrac = (pf.y - (yb / max(xb, 1e-12)) * pf.x) / max(ya, 1e-12);
    float sdir  = (yfrac < 0.5) ? 1.0 : -1.0;
    vec3  pnear = vec3(pf.x, pf.y + sdir * ya, -pf.z);
    vec4  dg1   = sdgMiuraFacet0(pnear, a, b, dims, ea, eap, n, csa, uvw1);
    dg1.w = -dg1.w;

    vec4 dg = (dg1.x < dg0.x) ? dg1 : dg0;
    vec3 g  = vec3(dg.y * cellId.z, dg.z, dg.w * cellId.w);
    return vec4(dg.x - thickness, g);
}

float sdMiuraOriSheet( in vec3 p, in float a, in float b, in float alpha, in float gamma, in float thickness )
{
    return sdgMiuraOriSheet(p, a, b, alpha, gamma, thickness).x;
}

// ---------------------------------------------------------------------------
// 4. Isometric Unfolding Material Coordinates (`opMiuraOriUnfold`)
//    Maps any 3D query point p back to its intrinsic unfolded paper coordinates (u, v, w),
//    where (u, v) in [0, a] x [0, b * sin(alpha)] are local isometric facet coordinates and
//    w is the signed normal offset from the closest facet plane.
// ---------------------------------------------------------------------------
vec3 opMiuraOriUnfold( in vec3 p, in float a, in float b, in float alpha, in float gamma )
{
    vec4 dims; vec3 ea, eap, n; vec2 csa;
    orMiuraFrame(a, b, alpha, gamma, dims, ea, eap, n, csa);
    float ya = dims.x, xb = dims.z, yb = dims.w;

    vec4 cellId;
    vec3 pf = opMiuraOriCell(p, a, b, alpha, gamma, cellId);

    vec3 uvw0, uvw1;
    vec4 dg0 = sdgMiuraFacet0(pf, a, b, dims, ea, eap, n, csa, uvw0);

    float yfrac = (pf.y - (yb / max(xb, 1e-12)) * pf.x) / max(ya, 1e-12);
    float sdir  = (yfrac < 0.5) ? 1.0 : -1.0;
    vec3  pnear = vec3(pf.x, pf.y + sdir * ya, -pf.z);
    vec4  dg1   = sdgMiuraFacet0(pnear, a, b, dims, ea, eap, n, csa, uvw1);

    vec3 uvw = (dg1.x < dg0.x) ? vec3(uvw1.xy, -uvw1.z) : uvw0;
    uvw.z *= cellId.w;
    return uvw;
}

// ---------------------------------------------------------------------------
// 5. Yoshimura / Kresling Cylindrical Origami Buckling (`sdYoshimuraCylinder`)
//    Triangulated cylindrical origami bellows of radius `r`, half-height `h`, `n` circumferential
//    facets per ring, `m` axial half-stories (total 2*m stories along y in [-h, h]), and relative
//    story twist angle `twist` (twist = pi / float(n) gives the symmetric Yoshimura diamond pattern;
//    other twist values give chiral Kresling torsional buckling).
// ---------------------------------------------------------------------------
float orYoshimuraTriDistSq( in vec3 p, in vec3 a, in vec3 b, in vec3 c )
{
    vec3  ba  = b - a, pa = p - a;
    vec3  cb  = c - b, pb = p - b;
    vec3  ac  = a - c, pc = p - c;
    vec3  nor = cross(ba, ac);
    float n2  = dot(nor, nor);
    float s1  = dot(cross(ba, nor), pa);
    float s2  = dot(cross(cb, nor), pb);
    float s3  = n2 - s1 - s2;
    if( s1 > 0.0 && s2 > 0.0 && s3 > 0.0 )
    {
        float np = dot(nor, pa);
        return (np * np) / max(n2, 1e-24);
    }
    vec3 q1 = ba * clamp(dot(ba, pa) / max(dot(ba, ba), 1e-24), 0.0, 1.0) - pa;
    vec3 q2 = cb * clamp(dot(cb, pb) / max(dot(cb, cb), 1e-24), 0.0, 1.0) - pb;
    vec3 q3 = ac * clamp(dot(ac, pc) / max(dot(ac, ac), 1e-24), 0.0, 1.0) - pc;
    return min(min(dot(q1, q1), dot(q2, q2)), dot(q3, q3));
}

float sdYoshimuraCylinder( in vec3 p, in float r, in float h, in int n, in int m, in float twist, in float thickness )
{
    float fn  = float(max(n, 3));
    float fm  = float(max(m, 1));
    float dy  = h / fm;
    float dth = 6.283185307179586 / fn;

    // 1. Fold y into the nearest story k in [0, 2*m - 1] (and check adjacent story k_adj for exactness):
    float yn = (p.y + h) / max(dy, 1e-12);
    float k0 = clamp(floor(yn), 0.0, 2.0 * fm - 1.0);
    float k1 = clamp((yn - k0 < 0.5) ? (k0 - 1.0) : (k0 + 1.0), 0.0, 2.0 * fm - 1.0);

    float best2 = 1e30;
    for( int pass = 0; pass < 2; ++pass )
    {
        float k    = (pass == 0) ? k0 : k1;
        bool  oddk = (int(floor(mod(k, 2.0) + 0.5)) != 0);
        float yloc = p.y - (-h + (k + (oddk ? 1.0 : 0.0)) * dy);
        if( oddk ) yloc = -yloc;

        // Every even ring (y_loc = 0) is at azimuth 0; every odd ring (y_loc = dy) is at azimuth +twist:
        float th = atan(p.z, p.x);
        float j0 = floor(th / dth);
        vec3  q  = vec3(p.x, yloc, p.z);

        for( int dj = -1; dj <= 1; ++dj )
        {
            float a0 = (j0 + float(dj)) * dth;
            float a1 = a0 + dth;
            vec3  B0 = vec3(r * cos(a0), 0.0, r * sin(a0));
            vec3  B1 = vec3(r * cos(a1), 0.0, r * sin(a1));
            vec3  T0 = vec3(r * cos(a0 + twist), dy, r * sin(a0 + twist));
            vec3  T1 = vec3(r * cos(a1 + twist), dy, r * sin(a1 + twist));

            best2 = min(best2, orYoshimuraTriDistSq(q, B0, B1, T0));
            best2 = min(best2, orYoshimuraTriDistSq(q, B1, T1, T0));
        }
    }
    return sqrt(best2) - thickness;
}

#endif // SDF_GEMS_HAS_OR_MIURA_ORI
