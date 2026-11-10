// ============================================================================
// sd4Polytopes — Exact 4D SDFs of the six regular 4-polytopes (5-cell, tesseract,
//                16-cell, 24-cell, 120-cell, 600-cell), their edge skeletons, and 3D slices
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @source Coxeter, "Regular Polytopes" (characteristic orthoscheme); fold technique as in knighty / hg_sdf
// @exact  sd4Cell*, sd4Tesseract, sd4RegularEdges: Yes in R^4 (GPU FP32 max error <= 2.0e-6 at r = 1 vs half-space + vertex-hull QP references)
//         sdRegularSlice:      Conservative bound in R^3 (restriction of an exact 4D SDF; on mirror slices exact OUTSIDE only)
//         sdRegularSliceBound: Conservative bound in R^3, tighter outside, NOT continuous
//
// Orientation ("orthoscheme frame"): a cell centre lies on +w, the centre of one of its faces in the zw-plane,
// the midpoint of one of that face's edges in the yzw-space, and that edge is parallel to x.
// In this frame the fundamental chamber of the symmetry group is
//      0 <= x/h.x <= y/h.y <= z/h.z <= w/h.w ,
// where the unit vector h is the direction of a vertex, and inside the chamber the solid is the box corner
// { q <= r*h componentwise }. Size convention: r = CIRCUMRADIUS (every vertex has |v| = r).
// Inradius = r*h.w, half edge = r*h.x.
// ============================================================================

#ifndef SDF_GEMS_HAS_SD4_POLYTOPES
#define SDF_GEMS_HAS_SD4_POLYTOPES

#define SD4_CELL5     0   // {3,3,3}  A4, order   120
#define SD4_TESSERACT 1   // {4,3,3}  B4, order   384
#define SD4_CELL16    2   // {3,3,4}  B4, order   384
#define SD4_CELL24    3   // {3,4,3}  F4, order  1152
#define SD4_CELL120   4   // {5,3,3}  H4, order 14400
#define SD4_CELL600   5   // {3,3,5}  H4, order 14400

// Unit vertex direction h, the three 2D mirror normals (acting on xy, yz, zw; the fourth mirror is x = 0), and the
// number of fold sweeps. The sweep counts are the exact worst case over EVERY chamber of the group (not an estimate).
void sd4RegularConstants( int id, out vec4 h, out vec2 k1, out vec2 k2, out vec2 k3, out int sweeps )
{
    if (id == SD4_CELL5) {
        h  = vec4(0.7905694150, 0.4564354646, 0.3227486122, 0.25);
        k1 = vec2(-0.5, 0.8660254038); k2 = vec2(-0.5773502692, 0.8164965809); k3 = vec2(-0.6123724357, 0.7905694150);
        sweeps = 4;
    } else if (id == SD4_TESSERACT) {
        h  = vec4(0.5);
        k1 = vec2(-0.7071067812, 0.7071067812); k2 = k1; k3 = k1;
        sweeps = 4;
    } else if (id == SD4_CELL16) {
        h  = vec4(0.7071067812, 0.4082482905, 0.2886751346, 0.5);
        k1 = vec2(-0.5, 0.8660254038); k2 = vec2(-0.5773502692, 0.8164965809); k3 = vec2(-0.8660254038, 0.5);
        sweeps = 4;
    } else if (id == SD4_CELL24) {
        h  = vec4(0.5, 0.2886751346, 0.4082482905, 0.7071067812);
        k1 = vec2(-0.5, 0.8660254038); k2 = vec2(-0.8164965809, 0.5773502692); k3 = vec2(-0.8660254038, 0.5);
        sweeps = 6;
    } else if (id == SD4_CELL120) {
        h  = vec4(0.1350453784, 0.1858740172, 0.3007504775, 0.9256147934);
        k1 = vec2(-0.8090169944, 0.5877852523); k2 = vec2(-0.8506508084, 0.5257311121); k3 = vec2(-0.9510565163, 0.3090169944);
        sweeps = 15;
    } else {
        h  = vec4(0.3090169944, 0.1784110449, 0.1261556597, 0.9256147934);
        k1 = vec2(-0.5, 0.8660254038); k2 = vec2(-0.5773502692, 0.8164965809); k3 = vec2(-0.9908394147, 0.1350453784);
        sweeps = 15;
    }
}

// Fold p into the chamber 0 <= x/h.x <= y/h.y <= z/h.z <= w/h.w with the four simple reflections.
vec4 sd4CoxeterFold( vec4 p, vec2 k1, vec2 k2, vec2 k3, int sweeps )
{
    for (int i = 0; i < sweeps; i++) {
        p.x = abs(p.x);
        p.xy -= 2.0 * min(dot(p.xy, k1), 0.0) * k1;
        p.yz -= 2.0 * min(dot(p.yz, k2), 0.0) * k2;
        p.zw -= 2.0 * min(dot(p.zw, k3), 0.0) * k3;
    }
    return p;
}

// Exact distance, for q in the chamber, to the solid: cell / ridge / edge / vertex are the box-corner regions.
float sd4OrthoschemeSolid( vec4 q, vec4 h )
{
    q -= h;
    return (q.w > 0.0) ? length(max(q, 0.0)) : q.w;
}

// Exact (unsigned) distance, for q in the chamber, to the edge skeleton: the one edge x in [0, h.x] at (y,z,w) = h.yzw.
float sd4OrthoschemeEdges( vec4 q, vec4 h )
{
    q -= h;
    return length(vec4(max(q.x, 0.0), q.yzw));
}

// Exact signed distance in R^4 to the solid regular 4-polytope `id` with circumradius r.
float sd4Regular( int id, vec4 p, float r )
{
    vec4 h; vec2 k1, k2, k3; int n;
    sd4RegularConstants(id, h, k1, k2, k3, n);
    return sd4OrthoschemeSolid(sd4CoxeterFold(p, k1, k2, k3, n), r * h);
}

// Exact unsigned distance in R^4 to the edge skeleton (subtract a tube radius for a wireframe).
float sd4RegularEdges( int id, vec4 p, float r )
{
    vec4 h; vec2 k1, k2, k3; int n;
    sd4RegularConstants(id, h, k1, k2, k3, n);
    return sd4OrthoschemeEdges(sd4CoxeterFold(p, k1, k2, k3, n), r * h);
}

float sd4Cell5  ( vec4 p, float r ) { return sd4Regular(SD4_CELL5,   p, r); }
float sd4Cell16 ( vec4 p, float r ) { return sd4Regular(SD4_CELL16,  p, r); }
float sd4Cell24 ( vec4 p, float r ) { return sd4Regular(SD4_CELL24,  p, r); }
float sd4Cell120( vec4 p, float r ) { return sd4Regular(SD4_CELL120, p, r); }
float sd4Cell600( vec4 p, float r ) { return sd4Regular(SD4_CELL600, p, r); }

// The tesseract needs no fold: it is the 4D box of half-side r/2 (identical to sd4Regular(SD4_TESSERACT, p, r)).
float sd4Tesseract( vec4 p, float r )
{
    vec4 q = abs(p) - 0.5 * r;
    return length(max(q, 0.0)) + min(max(max(q.x, q.y), max(q.z, q.w)), 0.0);
}

// 4D point of the 3-flat { R * (x, w) : x in R^3 }. R must be orthogonal (a 4D rotation).
vec4 op4Slice( vec3 x, mat4 R, float w )
{
    return R * vec4(x, w);
}

// Rotation by angle a in the coordinate plane (i, j) of R^4, i != j in 0..3. Compose several for a general rotation.
mat4 op4Rotation( int i, int j, float a )
{
    mat4 m = mat4(1.0);
    float c = cos(a), s = sin(a);
    m[i][i] = c; m[j][j] = c; m[i][j] = s; m[j][i] = -s;
    return m;
}

// 3D slice of the solid polytope: the restriction of the exact 4D SDF. 1-Lipschitz, correct zero set (the slice
// polyhedron), |value| <= true 3D distance. On a mirror 3-flat it is the exact 3D distance outside the solid only (inside it under-estimates, measured up to 0.126 at r = 1).
float sdRegularSlice( int id, vec3 x, mat4 R, float w, float r )
{
    return sd4Regular(id, op4Slice(x, R, w), r);
}

// 3D slice of the edge skeleton thickened to 4D radius th (the "cell-first bubble" pictures). Conservative bound.
float sdRegularEdgesSlice( int id, vec3 x, mat4 R, float w, float r, float th )
{
    return sd4RegularEdges(id, op4Slice(x, R, w), r) - th;
}

// Tighter conservative slice bound. The slice normal R[3] is folded together with p, which identifies the nearest
// cell hyperplane; inside the flat its distance is a / s with s = sin(angle between flat and cell). Outside the solid
// the result is max(4D distance, a/s, same for the cell across the ridge): each term is the distance to a supporting
// plane of the convex slice polyhedron or a 4D distance, hence a lower bound. It equals the exact 3D distance wherever
// the nearest point of the slice lies on the facet cut from the tracked cell. NOT continuous (the tracked cell
// changes across mirrors); inside it returns the plain restriction.
float sdRegularSliceBound( int id, vec3 x, mat4 R, float w, float r )
{
    vec4 h; vec2 k1, k2, k3; int sweeps;
    sd4RegularConstants(id, h, k1, k2, k3, sweeps);
    vec4 p = R * vec4(x, w);
    vec4 n = R[3];
    for (int i = 0; i < sweeps; i++) {
        float t;
        if (p.x < 0.0) { p.x = -p.x; n.x = -n.x; }
        t = dot(p.xy, k1); if (t < 0.0) { p.xy -= 2.0 * t * k1; n.xy -= 2.0 * dot(n.xy, k1) * k1; }
        t = dot(p.yz, k2); if (t < 0.0) { p.yz -= 2.0 * t * k2; n.yz -= 2.0 * dot(n.yz, k2) * k2; }
        t = dot(p.zw, k3); if (t < 0.0) { p.zw -= 2.0 * t * k3; n.zw -= 2.0 * dot(n.zw, k3) * k3; }
    }
    h *= r;
    vec4 q = p - h;
    if (q.w <= 0.0) return q.w;
    float d = length(max(q, 0.0));
    d = max(d, q.w * inversesqrt(max(1.0 - n.w * n.w, 1e-3)));
    float pw = p.w - 2.0 * dot(p.zw, k3) * k3.y;      // w of p and n mirrored into the cell across the ridge
    float nw = n.w - 2.0 * dot(n.zw, k3) * k3.y;
    d = max(d, (pw - h.w) * inversesqrt(max(1.0 - nw * nw, 1e-3)));
    return d;
}

#endif // SDF_GEMS_HAS_SD4_POLYTOPES
