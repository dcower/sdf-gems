// ============================================================================
// sdSnubCube & sdgSnubCube — Exact 0/1-Sqrt 3D Chiral Archimedean Snub Cube (38 Faces, 24 Vertices)
//
// @author  Claude Opus 5.5 (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & all 38 face prisms,
//          1 sqrt in the edge/vertex cones)
// @novelty Exact Euclidean SDF and analytical gradient for the chiral snub cube (6 squares,
//          8 + 24 equilateral triangles, 60 edges, 24 vertices; rotation group O only, no
//          mirror planes) via two results:
//            1. Chiral Parity Fold: fold with all 48 O_h reflections (abs + 3-swap sort) while
//               tracking the reflection-count parity; odd parity is undone by one fixed wall
//               reflection (swap y,z), landing in a double chamber D = C u s_yz(C) that is a
//               fundamental domain of the rotation group O. Inside D only the 5 faces of the
//               single vertex star v0 = (1, 1/t^2, 1/t) r are ever the nearest face plane.
//            2. Nearest-Vertex Tangent-Cone Lemma: for every exterior point, the closest
//               surface point lies in the star of the closest vertex, so the exact exterior
//               distance equals the distance to the 5-face tangent cone at that vertex. In D
//               the closest vertex is one of 3 (v0, v16, v8), each mapped to v0 by a swizzle.
//
// Chirality: this is the snub cube whose vertices are the even permutations of
//   (+-1, +-1/t, +-t) with an even number of plus signs together with the odd permutations
//   with an odd number of plus signs (scaled by r/t; t = tribonacci constant 1.83928676).
//   The mirror-image enantiomer is sdSnubCube(p * vec3(-1.0, 1.0, 1.0), r).
//
// Parameters:
//   p : 3D query point
//   r : inradius to the 6 axial square faces (triangle inradius = 1.06191270 * r,
//       edge length = 0.87518657 * r, circumradius = 1.17599990 * r)
// ============================================================================

#ifndef SDF_GEMS_SD_SNUB_CUBE_GLSL
#define SDF_GEMS_SD_SNUB_CUBE_GLSL

float sdSnubCube( vec3 p, float r )
{
    // 1. Chiral parity fold into the rotation-group fundamental domain D = C u s_yz(C).
    bool odd = (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    if( p.x < p.y ) { p.xy = p.yx; odd = !odd; }
    if( p.y < p.z ) { p.yz = p.zy; odd = !odd; }
    if( p.x < p.y ) { p.xy = p.yx; odd = !odd; }
    if( odd ) p.yz = p.zy;

    // 2. Interior (exact) / plane bound: square + axial triangle + 3 snub triangles at v0.
    const float K = 1.06191270; // triangle inradius / square inradius
    float dm = max(p.x - r,
                   max(max(dot(p, vec3(0.79846143, 0.57735027, 0.17066344)),
                           dot(p, vec3(0.57735027, 0.57735027, 0.57735027))),
                       max(dot(p, vec3(0.57735027, 0.17066344, 0.79846143)),
                           dot(p, vec3(0.79846143, -0.17066344, 0.57735027)))) - K * r);
    if( dm <= 0.0 ) return dm;

    // 3. Nearest vertex among v0 = (1, a, b), v16 = (b, 1, a), v8 = (1, b, -a); rotate it to v0.
    float s16 = dot(p, vec3(-0.45631099, 0.70440226, -0.24809127));
    float s8  = dot(p, vec3( 0.0,        0.24809127, -0.83928675));
    if( max(s16, s8) > 0.0 ) p = (s16 >= s8) ? p.yzx : vec3(p.x, -p.z, p.y);

    // 4. Exact distance to the 5-face tangent cone at v0.
    //    Face tests accept h > -1e-5 r: on the surface, FP32 rounding can give dm > 0 above
    //    but h < 0 here; a point inside a face prism within 1e-5 r of its plane is within
    //    1e-5 r of the face, so this only affects a rounding shell (error < 1e-5 r).
    vec3  q  = p - vec3(1.0, 0.29559774, 0.54368901) * r;
    float h0 = q.x;
    float h1 = dot(q, vec3(0.79846143,  0.57735027, 0.17066344));
    float h2 = dot(q, vec3(0.57735027,  0.57735027, 0.57735027));
    float h3 = dot(q, vec3(0.57735027,  0.17066344, 0.79846143));
    float h4 = dot(q, vec3(0.79846143, -0.17066344, 0.57735027));
    if( h0 > -1e-5 * r && dot(q.yz, vec2(-0.28347244, 0.95898038)) <= 0.0
                 && dot(q.yz, vec2( 0.95898038, 0.28347244)) <= 0.0 ) return h0;
    if( h1 > -1e-5 * r && dot(q, vec3( 0.60204596, -0.76570885, -0.22634181)) <= 0.0
                 && dot(q, vec3(-0.30102298,  0.13736009,  0.94367228)) <= 0.0 ) return h1;
    if( h2 > -1e-5 * r && dot(q, vec3( 0.62834876,  0.13736009, -0.76570885)) <= 0.0
                 && dot(q, vec3( 0.13736009, -0.76570885,  0.62834876)) <= 0.0 ) return h2;
    if( h3 > -1e-5 * r && dot(q, vec3( 0.13736009,  0.94367228, -0.30102298)) <= 0.0
                 && dot(q, vec3( 0.62834876, -0.71733047, -0.30102298)) <= 0.0 ) return h3;
    if( h4 > -1e-5 * r && dot(q, vec3(-0.30102298,  0.71733047,  0.62834876)) <= 0.0
                 && dot(q, vec3( 0.60204596,  0.22634181, -0.76570885)) <= 0.0 ) return h4;

    // Edge rays (5 edges at v0) or the vertex itself: project off the best edge direction.
    vec3  e0 = vec3( 0.0,        -0.95898038, -0.28347244);
    vec3  e1 = vec3( 0.0,         0.28347244, -0.95898038);
    vec3  e2 = vec3(-0.52138710,  0.80485954, -0.28347244);
    vec3  e3 = vec3(-0.80485954,  0.28347244,  0.52138710);
    vec3  e4 = vec3(-0.52138710, -0.67550795,  0.52138710);
    float u = 0.0;
    vec3  w = vec3(0.0);
    float t;
    t = dot(q, e0); if( t > u ) { u = t; w = e0; }
    t = dot(q, e1); if( t > u ) { u = t; w = e1; }
    t = dot(q, e2); if( t > u ) { u = t; w = e2; }
    t = dot(q, e3); if( t > u ) { u = t; w = e3; }
    t = dot(q, e4); if( t > u ) { u = t; w = e4; }
    return length(q - u * w);
}

vec4 sdgSnubCube( vec3 p, float r )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    bool odd = (p.x < 0.0) ^^ (p.y < 0.0) ^^ (p.z < 0.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) { p.xy = p.yx; odd = !odd; }
    bool s2 = (p.y < p.z); if( s2 ) { p.yz = p.zy; odd = !odd; }
    bool s3 = (p.x < p.y); if( s3 ) { p.xy = p.yx; odd = !odd; }
    if( odd ) p.yz = p.zy;

    const float K = 1.06191270;
    const vec3 n1 = vec3(0.79846143,  0.57735027, 0.17066344);
    const vec3 n2 = vec3(0.57735027,  0.57735027, 0.57735027);
    const vec3 n3 = vec3(0.57735027,  0.17066344, 0.79846143);
    const vec3 n4 = vec3(0.79846143, -0.17066344, 0.57735027);

    float d;
    vec3  g;
    float a0 = p.x - r;
    float a1 = dot(p, n1) - K * r;
    float a2 = dot(p, n2) - K * r;
    float a3 = dot(p, n3) - K * r;
    float a4 = dot(p, n4) - K * r;
    float dm = max(max(max(a0, a1), max(a2, a3)), a4);
    int   rot = 0;
    if( dm <= 0.0 )
    {
        d = dm;
        g = (a0 >= dm) ? vec3(1.0, 0.0, 0.0) : (a1 >= dm) ? n1 : (a2 >= dm) ? n2 : (a3 >= dm) ? n3 : n4;
    }
    else
    {
        float s16 = dot(p, vec3(-0.45631099, 0.70440226, -0.24809127));
        float s8  = dot(p, vec3( 0.0,        0.24809127, -0.83928675));
        if( max(s16, s8) > 0.0 )
        {
            if( s16 >= s8 ) { p = p.yzx;                 rot = 1; }
            else            { p = vec3(p.x, -p.z, p.y);  rot = 2; }
        }
        vec3  q  = p - vec3(1.0, 0.29559774, 0.54368901) * r;
        float h0 = q.x;
        float h1 = dot(q, n1);
        float h2 = dot(q, n2);
        float h3 = dot(q, n3);
        float h4 = dot(q, n4);
        if( h0 > -1e-5 * r && dot(q.yz, vec2(-0.28347244, 0.95898038)) <= 0.0
                     && dot(q.yz, vec2( 0.95898038, 0.28347244)) <= 0.0 )      { d = h0; g = vec3(1.0, 0.0, 0.0); }
        else if( h1 > -1e-5 * r && dot(q, vec3( 0.60204596, -0.76570885, -0.22634181)) <= 0.0
                          && dot(q, vec3(-0.30102298,  0.13736009,  0.94367228)) <= 0.0 ) { d = h1; g = n1; }
        else if( h2 > -1e-5 * r && dot(q, vec3( 0.62834876,  0.13736009, -0.76570885)) <= 0.0
                          && dot(q, vec3( 0.13736009, -0.76570885,  0.62834876)) <= 0.0 ) { d = h2; g = n2; }
        else if( h3 > -1e-5 * r && dot(q, vec3( 0.13736009,  0.94367228, -0.30102298)) <= 0.0
                          && dot(q, vec3( 0.62834876, -0.71733047, -0.30102298)) <= 0.0 ) { d = h3; g = n3; }
        else if( h4 > -1e-5 * r && dot(q, vec3(-0.30102298,  0.71733047,  0.62834876)) <= 0.0
                          && dot(q, vec3( 0.60204596,  0.22634181, -0.76570885)) <= 0.0 ) { d = h4; g = n4; }
        else
        {
            vec3  e0 = vec3( 0.0,        -0.95898038, -0.28347244);
            vec3  e1 = vec3( 0.0,         0.28347244, -0.95898038);
            vec3  e2 = vec3(-0.52138710,  0.80485954, -0.28347244);
            vec3  e3 = vec3(-0.80485954,  0.28347244,  0.52138710);
            vec3  e4 = vec3(-0.52138710, -0.67550795,  0.52138710);
            float u = 0.0;
            vec3  w = vec3(0.0);
            float t;
            t = dot(q, e0); if( t > u ) { u = t; w = e0; }
            t = dot(q, e1); if( t > u ) { u = t; w = e1; }
            t = dot(q, e2); if( t > u ) { u = t; w = e2; }
            t = dot(q, e3); if( t > u ) { u = t; w = e3; }
            t = dot(q, e4); if( t > u ) { u = t; w = e4; }
            vec3  v    = q - u * w;
            float v2   = dot(v, v);
            float invV = inversesqrt(max(v2, 1e-30));
            d = v2 * invV;
            g = (v2 > 1e-30) ? v * invV : vec3(0.57735027);
        }
    }

    // Undo the rotation, the parity swap, the sort and the sign flips (transpose of each map).
    if( rot == 1 ) g = g.zxy;
    if( rot == 2 ) g = vec3(g.x, g.z, -g.y);
    if( odd ) g.yz = g.zy;
    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4(d, sg * g);
}

#endif // SDF_GEMS_SD_SNUB_CUBE_GLSL
