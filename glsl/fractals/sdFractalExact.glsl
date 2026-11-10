// ============================================================================
// sdSierpinski2D, sdSierpinskiTetra, sdMengerExact, sdMengerGreedy (+ LOD variants)
// — Euclidean SDFs of level-n self-similar sets by symmetry-pruned descent of the IFS tree
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  Measured numbers and proofs: docs/fractals/sdFractalExact.md
//         sdSierpinski2D    : Exact (proved; nearest-vertex descent, one triangle SDF in total).
//         sdSierpinskiTetra : Exact (proved; nearest-vertex descent, one tetrahedron SDF in total).
//         sdMengerExact     : Exact for n <= 5 with the default beam of 16 (the folded tree has at most
//                             2^(n-1) nodes per level, so the beam cannot overflow). For larger n it is
//                             a conservative bound (outside: d <= distance; inside: |d| <= depth).
//         sdMengerGreedy    : Two box distances per level. Agreed with sdMengerExact at every sampled
//                             point (see doc); NOT proved, so treat it as approximate.
//         *Lod variants     : Conservative (never above the level-n distance); within eps of it wherever
//                             the result is positive.
//
// A level-n approximant is a finite union of similar copies of one primitive, so its distance is a min
// over copies, and the copies form a tree. Each node is symmetric under a reflection group, so the nearest
// point of the node to p lies in the closed chamber of that group containing p. For the Sierpinski sets
// the chamber meets one child only (the one at the nearest vertex): the descent is a single path. For the
// Menger sponge the chamber x >= y >= z >= 0 meets two kept children (corner and edge) and two removed
// ones (centre and face centre): the tree has branching factor 2 instead of 20.
// ============================================================================

#ifndef SDF_GEMS_SD_FRACTAL_EXACT_GLSL
#define SDF_GEMS_SD_FRACTAL_EXACT_GLSL

// Candidates kept per level by sdMengerExact. 2^(n-1) guarantees exactness at level n.
#ifndef SDF_FRACTAL_MENGER_BEAM
#define SDF_FRACTAL_MENGER_BEAM 16
#endif

// ---------------------------------------------------------------------------
// Primitives
// ---------------------------------------------------------------------------

// Exact equilateral triangle, centroid at the origin, one vertex at (0, R).
float sdFractalEqTri( in vec2 p, in float R )
{
    const float k = 1.7320508075688772;
    float r = 0.8660254037844386 * R;
    p.x = abs(p.x) - r;
    p.y = p.y + r / k;
    if( p.x + k * p.y > 0.0 ) p = vec2(p.x - k * p.y, -k * p.x - p.y) * 0.5;
    p.x -= clamp(p.x, -2.0 * r, 0.0);
    return (p.y > 0.0) ? -length(p) : length(p);
}

// Exact regular tetrahedron with vertices s*(1,1,1), s*(1,-1,-1), s*(-1,1,-1), s*(-1,-1,1)
// (parity-sorted Weyl chamber of Td: one face plane, one edge, one vertex).
float sdFractalTetra( in vec3 p, in float s )
{
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;
    p.z *= sig;
    float m = p.x + p.y - p.z - s;
    if( m <= 1.5 * (p.y - p.z) ) return m * 0.5773502691896258;
    float t = min(0.5 * (p.y + p.z), s);
    return length(vec3(p.x - s, p.y - t, p.z - t));
}

// Unsigned distance to the cube of half-size s.
float udFractalBox( in vec3 p, in float s )
{
    return length(max(abs(p) - s, 0.0));
}

// Fold into the chamber x >= y >= z >= 0 of the cube group.
vec3 opFractalCubeFold( in vec3 p )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;
    return p;
}

// Level of detail: the smallest level k <= n whose hull distance is within eps of the level-n distance,
// i.e. h * ratio^k <= eps (h = sup over the root hull of the distance to the attractor).
int sdFractalLodLevel( in int n, in float h, in float invRatio, in float eps )
{
    if( eps <= 0.0 ) return n;
    return clamp(int(ceil(log(max(h / eps, 1.0)) / log(invRatio))), 0, n);
}

// ---------------------------------------------------------------------------
// 1. Sierpinski triangle (2D), level n: union of 3^n filled triangles.
//    Root triangle: circumradius 1, vertices (0, 1), (+-sqrt(3)/2, -1/2).
// ---------------------------------------------------------------------------
float sdSierpinski2D( in vec2 p, in int n )
{
    float sc = 1.0;
    for( int k = 0; k < n; k++ )
    {
        // nearest vertex = largest dot product
        vec2 v = (p.y > 0.5773502691896258 * abs(p.x)) ? vec2(0.0, 1.0)
               : vec2((p.x > 0.0) ? 0.8660254037844386 : -0.8660254037844386, -0.5);
        p   = 2.0 * p - v;
        sc *= 0.5;
    }
    return sdFractalEqTri(p, 1.0) * sc;
}

// ---------------------------------------------------------------------------
// 2. Sierpinski tetrahedron, level n: union of 4^n tetrahedra.
//    Root: vertices (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1) (circumradius sqrt(3)).
// ---------------------------------------------------------------------------
float sdSierpinskiTetra( in vec3 p, in int n )
{
    float sc = 1.0;
    for( int k = 0; k < n; k++ )
    {
        float d0 =  p.x + p.y + p.z;
        float d1 =  p.x - p.y - p.z;
        float d2 = -p.x + p.y - p.z;
        float d3 = -p.x - p.y + p.z;
        vec3  v  = vec3(1.0, 1.0, 1.0);
        float dm = d0;
        if( d1 > dm ) { dm = d1; v = vec3( 1.0, -1.0, -1.0); }
        if( d2 > dm ) { dm = d2; v = vec3(-1.0,  1.0, -1.0); }
        if( d3 > dm ) { dm = d3; v = vec3(-1.0, -1.0,  1.0); }
        p   = 2.0 * p - v;
        sc *= 0.5;
    }
    return sdFractalTetra(p, 1.0) * sc;
}

// Conservative LOD: result <= sdSierpinskiTetra(p, n) everywhere, and the gap is <= eps wherever the
// result is positive (inside the coarse hull the result is negative: the coarse solid stands in).
float sdSierpinskiTetraLod( in vec3 p, in int n, in float eps )
{
    return sdSierpinskiTetra(p, sdFractalLodLevel(n, 0.6123724356957945, 2.0, eps)); // h = sqrt(3/8)
}

// ---------------------------------------------------------------------------
// 3. Menger sponge, level n, in the cube [-1, 1]^3: union of 20^n cubes.
//    Outside the solid: distance to the union of cubes.
//    Inside the solid : minus the distance to the complement, i.e. to the nearest removed cube of any
//                       level 1..n or to the outside of the root cube.
//    Each candidate stores the query point in its node's frame (cq) and the node's box distance (cl).
// ---------------------------------------------------------------------------
float sdFractalMengerOuter( in vec3 p, in int n, in float eps )
{
    const int   K = SDF_FRACTAL_MENGER_BEAM;
    const float H = 0.4714045207910317; // sqrt(2)/3: sup over a cube of half-size s of the distance to the sponge / s
    vec3  cq[K];
    float cl[K];
    cq[0] = p;
    cl[0] = udFractalBox(p, 1.0);
    if( n <= 0 || H <= eps ) return cl[0];
    int   cnt    = 1;
    float s      = 1.0;
    float ub     = cl[0] + H;
    float pruned = 1e30;
    for( int k = 1; k <= n; k++ )
    {
        float sc   = s / 3.0;
        bool  leaf = (k == n) || (H * sc <= eps);
        float hh   = leaf ? 0.0 : H * sc;
        vec3  nq[K];
        float nl[K];
        int   m    = 0;
        float best = 1e30;
        for( int i = 0; i < cnt; i++ )
        {
            if( cl[i] >= ub ) continue;
            vec3 a = opFractalCubeFold(cq[i]);
            for( int j = 0; j < 2; j++ )
            {
                vec3  q  = a - (2.0 * sc) * vec3(1.0, 1.0, (j == 0) ? 1.0 : 0.0); // corner, edge
                float lb = udFractalBox(q, sc);
                if( leaf ) { best = min(best, lb); continue; }
                ub = min(ub, lb + hh);
                if( lb >= ub ) continue;
                int pos;
                if( m < K ) { pos = m; m++; }
                else if( lb < nl[K - 1] ) { pruned = min(pruned, nl[K - 1]); pos = K - 1; }
                else { pruned = min(pruned, lb); continue; }
                while( pos > 0 && nl[pos - 1] > lb ) { nl[pos] = nl[pos - 1]; nq[pos] = nq[pos - 1]; pos--; }
                nl[pos] = lb; nq[pos] = q;
            }
        }
        if( leaf ) return min(best, pruned);
        for( int i = 0; i < m; i++ ) { cq[i] = nq[i]; cl[i] = nl[i]; }
        cnt = m;
        s   = sc;
    }
    return pruned;
}

// Distance from a point of the solid to its complement (>= 0).
float sdFractalMengerInner( in vec3 p, in int n )
{
    const int K = SDF_FRACTAL_MENGER_BEAM;
    vec3  cq[K];
    float cl[K];
    cq[0] = p;
    cl[0] = 0.0;
    int   cnt    = 1;
    float s      = 1.0;
    vec3  w      = 1.0 - abs(p);
    float best   = max(min(w.x, min(w.y, w.z)), 0.0);
    float pruned = 1e30;
    for( int k = 1; k <= n; k++ )
    {
        float sc = s / 3.0;
        // removed children (centre, face centre) first: they tighten `best` before the kept ones are filtered
        for( int i = 0; i < cnt; i++ )
        {
            if( cl[i] >= best ) continue;
            vec3 a = opFractalCubeFold(cq[i]);
            best = min(best, udFractalBox(a, sc));
            best = min(best, udFractalBox(a - vec3(2.0 * sc, 0.0, 0.0), sc));
        }
        if( k == n ) break;
        vec3  nq[K];
        float nl[K];
        int   m = 0;
        for( int i = 0; i < cnt; i++ )
        {
            if( cl[i] >= best ) continue;
            vec3 a = opFractalCubeFold(cq[i]);
            for( int j = 0; j < 2; j++ )
            {
                vec3  q  = a - (2.0 * sc) * vec3(1.0, 1.0, (j == 0) ? 1.0 : 0.0);
                float lb = udFractalBox(q, sc);
                if( lb >= best ) continue;
                int pos;
                if( m < K ) { pos = m; m++; }
                else if( lb < nl[K - 1] ) { pruned = min(pruned, nl[K - 1]); pos = K - 1; }
                else { pruned = min(pruned, lb); continue; }
                while( pos > 0 && nl[pos - 1] > lb ) { nl[pos] = nl[pos - 1]; nq[pos] = nq[pos - 1]; pos--; }
                nl[pos] = lb; nq[pos] = q;
            }
        }
        for( int i = 0; i < m; i++ ) { cq[i] = nq[i]; cl[i] = nl[i]; }
        cnt = m;
        s   = sc;
    }
    return min(best, pruned);
}

float sdMengerExactLod( in vec3 p, in int n, in float eps )
{
    float d = sdFractalMengerOuter(p, n, eps);
    // 2e-6: a point of the solid can come out a few ulp outside its leaf cube after n folds in FP32
    if( d > 2e-6 ) return d;
    return -sdFractalMengerInner(p, n);
}

float sdMengerExact( in vec3 p, in int n )
{
    return sdMengerExactLod(p, n, 0.0);
}

// Single-path variant: always follow the nearer of the two kept children (2 box distances per level
// outside, 4 inside). Equal to sdMengerExact at every point sampled so far; no proof.
float sdMengerGreedy( in vec3 p, in int n )
{
    vec3  q  = p;
    float s  = 1.0;
    float d  = udFractalBox(p, 1.0);
    for( int k = 0; k < n; k++ )
    {
        float sc = s / 3.0;
        vec3  a  = opFractalCubeFold(q);
        vec3  q0 = a - vec3(2.0 * sc);
        vec3  q1 = a - vec3(2.0 * sc, 2.0 * sc, 0.0);
        float l0 = udFractalBox(q0, sc);
        float l1 = udFractalBox(q1, sc);
        if( l0 <= l1 ) { q = q0; d = l0; } else { q = q1; d = l1; }
        s = sc;
    }
    if( d > 2e-6 ) return d;
    vec3  w    = 1.0 - abs(p);
    float best = max(min(w.x, min(w.y, w.z)), 0.0);
    q = p;
    s = 1.0;
    for( int k = 1; k <= n; k++ )
    {
        float sc = s / 3.0;
        vec3  a  = opFractalCubeFold(q);
        best = min(best, udFractalBox(a, sc));
        best = min(best, udFractalBox(a - vec3(2.0 * sc, 0.0, 0.0), sc));
        vec3  q0 = a - vec3(2.0 * sc);
        vec3  q1 = a - vec3(2.0 * sc, 2.0 * sc, 0.0);
        q = (udFractalBox(q0, sc) <= udFractalBox(q1, sc)) ? q0 : q1;
        s = sc;
    }
    return -best;
}

#endif // SDF_GEMS_SD_FRACTAL_EXACT_GLSL
