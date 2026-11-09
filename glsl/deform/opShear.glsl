// ============================================================================
// opShear, opShearBound, sdgShearFinalize, opShearBox, sdgShearBox,
// sdRhombohedron & sdgRhombohedron — Unimodular Shear Singular-Value Bounds,
// Eikonal Chain-Rule Gradients, Exact 0/1-Sqrt Sheared Box & D_3d Rhombohedron
//
// @author  Discovered and derived using Claude Opus 5.5; audited and completed
//          (acute rhombohedron polar edge e3 for r/h < sqrt(2)/3) using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   opShearBound is the exact 1-sqrt 0-div unimodular singular-value
//          Lipschitz bound; sdgShearFinalize is 1st-order Eikonal-exact on the
//          surface; opShearBox, sdgShearBox, sdRhombohedron, and sdgRhombohedron
//          are exact Euclidean SDFs (`||grad d|| = 1` everywhere in R^3; 0 sqrts
//          in the interior and all 6 face prisms, 1 sqrt on edges/vertices)
// ============================================================================

#ifndef SDF_GEMS_OP_SHEAR_GLSL
#define SDF_GEMS_OP_SHEAR_GLSL

// 1. Standard 3D Linear Shear along x (or xz) proportional to y:
//    Maps p -> q = vec3(p.x - k * p.y, p.y, p.z).
vec3 opShear( in vec3 p, in float k )
{
    return vec3(p.x - k * p.y, p.y, p.z);
}

vec3 opShear( in vec3 p, in vec2 k )
{
    return vec3(p.x - k.x * p.y, p.y, p.z - k.y * p.y);
}

vec3 opShearXZ( in vec3 p, in vec2 k )
{
    return vec3(p.x - k.x * p.y, p.y, p.z - k.y * p.y);
}

// 2. Exact Unimodular Shear Singular-Value Bound (1 sqrt, 0 divisions):
//    The shear Jacobian J = [1, -k_x, 0; 0, 1, 0; 0, -k_z, 1] has det(J) = 1
//    and restricts on span(k, y) to M = [1, -||k||; 0, 1] with
//      sigma_max * sigma_min = 1,   sigma_max^2 + sigma_min^2 = 2 + ||k||^2.
//    Notice that naive column-norm scaling 1 / sqrt(1 + ||k||^2) STILL OVERSTEPS
//    by up to +14.4% at ||k|| = 1 because sigma_max(J) > sqrt(1 + ||k||^2).
//    Completing the square gives (sigma_max +/- sigma_min)^2 = {4 + ||k||^2, ||k||^2},
//    yielding the exact reciprocal maximum singular value with 0 divisions:
//      1 / sigma_max(J) = sigma_min(J) = 0.5 * (sqrt(4 + ||k||^2) - ||k||).
//    Usage: vec4 sp = opShearBound(p, k); float d = sdPrimitive(sp.xyz) * sp.w;
vec4 opShearBound( in vec3 p, in float k )
{
    return vec4(p.x - k * p.y, p.y, p.z, 0.5 * (sqrt(4.0 + k * k) - abs(k)));
}

vec4 opShearXZBound( in vec3 p, in vec2 k )
{
    float k2 = dot(k, k);
    return vec4(p.x - k.x * p.y, p.y, p.z - k.y * p.y, 0.5 * (sqrt(4.0 + k2) - sqrt(k2)));
}

vec4 opShearBound( in vec3 p, in vec2 k )
{
    return opShearXZBound(p, k);
}

// 3. First-Order Eikonal-Exact Chain-Rule Pullback (sdgShearFinalize):
//    Given q = opShear(p, k) and dg = sdgPrimitive(q), the exact pullback
//    gradient in p-space is g = J^T g_q = (g_qx, g_qy - k_x g_qx - k_z g_qz, g_qz).
//    Normalizing dg by inversesqrt(dot(g, g)) restores ||grad d|| = 1 on the
//    sheared surface and eliminates both overstepping and face understepping.
vec4 sdgShearFinalize( in float k, in vec4 dg )
{
    vec3  g    = vec3(dg.y, dg.z - k * dg.y, dg.w);
    float invG = inversesqrt(max(dot(g, g), 1e-30));
    return vec4(dg.x * invG, g * invG);
}

vec4 sdgShearFinalize( in vec3 p, in float k, in vec4 dg )
{
    return sdgShearFinalize(k, dg);
}

vec4 sdgShearFinalize( in vec2 k, in vec4 dg )
{
    vec3  g    = vec3(dg.y, dg.z - k.x * dg.y - k.y * dg.w, dg.w);
    float invG = inversesqrt(max(dot(g, g), 1e-30));
    return vec4(dg.x * invG, g * invG);
}

vec4 sdgShearFinalize( in vec3 p, in vec2 k, in vec4 dg )
{
    return sdgShearFinalize(k, dg);
}

// 4. Exact 0/1-Sqrt 3D Sheared Box (half-extents b, horizontal shear rate k):
//    Exact Euclidean SDF of { p in R^3 : |p.x - k * p.y| <= b.x, |p.y| <= b.y, |p.z| <= b.z }.
//    Evaluates 0 sqrts in the interior and all 6 exterior face Voronoi prisms,
//    and 1 sqrt on the 12 edges and 8 vertex cones.
float opShearBox( in vec3 p, in vec3 b, in float k )
{
    float sy   = (p.y < 0.0) ? -1.0 : 1.0;
    float ux0  = p.x - k * p.y;
    float sx   = (ux0 < 0.0) ? -1.0 : 1.0;
    float qx   = abs(ux0) - b.x;
    float dy   = abs(p.y) - b.y;
    float dz   = abs(p.z) - b.z;
    float l2   = 1.0 + k * k;
    float dx   = qx * inversesqrt(l2);
    float cx   = p.x - k * sy * b.y;
    float proj = p.y + k * (p.x - sx * b.x);

    if( (qx <= 0.0 && dy <= 0.0) ||
        (dz <= 0.0 && ((dy > 0.0 && abs(cx) <= b.x) || (qx > 0.0 && abs(proj) <= l2 * b.y))) )
        return max(max(dx, dy), dz);

    vec2  w  = vec2(cx - clamp(cx, -b.x, b.x), dy);
    float ty = clamp(proj / l2, -b.y, b.y);
    vec2  v  = vec2(p.x - (sx * b.x + k * ty), p.y - ty);
    float mz = max(dz, 0.0);
    return sqrt(min(dot(w, w), dot(v, v)) + mz * mz);
}

// 5. Exact 0/1-Sqrt 3D Sheared Box with Analytical Unit Gradient vec4(d, grad.xyz)
vec4 sdgShearBox( in vec3 p, in vec3 b, in float k )
{
    float sy   = (p.y < 0.0) ? -1.0 : 1.0;
    float sz   = (p.z < 0.0) ? -1.0 : 1.0;
    float ux0  = p.x - k * p.y;
    float sx   = (ux0 < 0.0) ? -1.0 : 1.0;
    float qx   = abs(ux0) - b.x;
    float dy   = abs(p.y) - b.y;
    float dz   = abs(p.z) - b.z;
    float l2   = 1.0 + k * k;
    float invL = inversesqrt(l2);
    float dx   = qx * invL;
    float cx   = p.x - k * sy * b.y;
    float proj = p.y + k * (p.x - sx * b.x);

    bool in_int  = (qx <= 0.0 && dy <= 0.0 && dz <= 0.0);
    bool in_zcap = (dz > 0.0 && qx <= 0.0 && dy <= 0.0);
    bool in_ycap = (dz <= 0.0 && dy > 0.0 && abs(cx) <= b.x);
    bool in_xcap = (dz <= 0.0 && qx > 0.0 && abs(proj) <= l2 * b.y);

    if( in_int || in_zcap || in_ycap || in_xcap )
    {
        float d     = max(max(dx, dy), dz);
        bool  win_z = in_zcap || (in_int && dz >= dx && dz >= dy);
        bool  win_y = in_ycap || (in_int && !win_z && dy >= dx);
        vec3  g     = win_z ? vec3(0.0, 0.0, sz)
                            : (win_y ? vec3(0.0, sy, 0.0) : vec3(sx * invL, -k * sx * invL, 0.0));
        return vec4(d, g);
    }

    vec2  w    = vec2(cx - clamp(cx, -b.x, b.x), sy * dy);
    float ty   = clamp(proj / l2, -b.y, b.y);
    vec2  v    = vec2(p.x - (sx * b.x + k * ty), p.y - ty);
    float w2   = dot(w, w), v2 = dot(v, v);
    vec3  e    = vec3((w2 <= v2) ? w : v, sz * max(dz, 0.0));
    float e2   = min(w2, v2) + e.z * e.z;
    float invE = inversesqrt(max(e2, 1e-30));
    vec3  g    = (e2 > 1e-30) ? (e * invE) : vec3(0.0, sy, 0.0);
    return vec4(e2 * invE, g);
}

// 6. Exact 0-Trig 0/1-Sqrt 3D Trigonal Rhombohedron (equatorial radius r, polar half-height h):
//    Calcite crystal / body-diagonal-stretched cube with D_3d antiprismatic symmetry (order 12),
//    polar apices at (0, +/-h, 0), and 6 equatorial zig-zag vertices at height +/-h/3 and radius r.
//    Reduces to an exact regular cube tilted along its (1,1,1) body diagonal at r = (2*sqrt(2)/3)*h.
//    Folds p into the 12-fold D_3d fundamental chamber where all 6 congruent rhombic faces collapse
//    to 1 plane and all 12 equal-length edges collapse to 3 edges sharing denominator r^2 + (4/9)*h^2.
float sdRhombohedron( in vec3 p, in float r, in float h )
{
    const float c = 0.5, s = 0.86602540;
    float x  = p.x, y = p.y, z = abs(p.z);
    float d3 = max(-s * x + c * z, 0.0);
    x += 2.0 * d3 * s;
    z  = abs(z - 2.0 * d3 * c);

    float u = x * c + z * s, v = x * s - z * c;
    float a = 1.33333333 * h;
    if( a * (x - u) + 2.0 * r * y < 0.0 )
    {
        float tx = x, tz = z;
        x = u; z = v; u = tx; v = tz; y = -y;
    }

    float n2 = a * a + r * r;
    float m  = a * x + r * (y - h);
    float w2 = a * (x - u) + 2.0 * r * y;
    if( a * s * m <= n2 * v && (0.5 * a * a + 2.0 * r * r) * m <= n2 * w2 )
        return m * inversesqrt(n2);

    float invLe2 = 1.0 / (r * r + 0.25 * a * a);
    vec3  e1  = vec3(0.5 * r, -0.5 * a, s * r);
    vec3  q1  = vec3(x, y - h, z);
    vec3  d1  = q1 - clamp(dot(q1, e1) * invLe2, 0.0, 1.0) * e1;

    vec3  e2  = vec3(0.5 * r, -0.5 * a, -s * r);
    vec3  q2  = vec3(x - 0.5 * r, y - 0.33333333 * h, z - s * r);
    vec3  d2  = q2 - clamp(dot(q2, e2) * invLe2, 0.0, 1.0) * e2;

    vec3  e3  = vec3(r, 0.5 * a, 0.0);
    vec3  q3  = vec3(x, y + h, z);
    vec3  d3e = q3 - clamp(dot(q3, e3) * invLe2, 0.0, 1.0) * e3;

    return sqrt(min(min(dot(d1, d1), dot(d2, d2)), dot(d3e, d3e)));
}

// 7. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) for D_3d Trigonal Rhombohedron
vec4 sdgRhombohedron( in vec3 p, in float r, in float h )
{
    const float c = 0.5, s = 0.86602540;
    float sz0 = (p.z < 0.0) ? -1.0 : 1.0;
    float x   = p.x, y = p.y, z = abs(p.z);
    float d3  = -s * x + c * z;
    bool  r0  = (d3 > 0.0);
    if( r0 ) { x += 2.0 * d3 * s; z -= 2.0 * d3 * c; }
    float sz1 = (z < 0.0) ? -1.0 : 1.0;
    z = abs(z);

    float u = x * c + z * s, v = x * s - z * c;
    float a = 1.33333333 * h;
    bool  swap = (a * (x - u) + 2.0 * r * y < 0.0);
    if( swap )
    {
        float tx = x, tz = z;
        x = u; z = v; u = tx; v = tz; y = -y;
    }

    float n2   = a * a + r * r;
    float invN = inversesqrt(n2);
    float m    = a * x + r * (y - h);
    float w2   = a * (x - u) + 2.0 * r * y;
    float d;
    vec3  g;
    if( a * s * m <= n2 * v && (0.5 * a * a + 2.0 * r * r) * m <= n2 * w2 )
    {
        d = m * invN;
        g = vec3(a * invN, r * invN, 0.0);
    }
    else
    {
        float invLe2 = 1.0 / (r * r + 0.25 * a * a);
        vec3  e1  = vec3(0.5 * r, -0.5 * a, s * r);
        vec3  q1  = vec3(x, y - h, z);
        vec3  d1  = q1 - clamp(dot(q1, e1) * invLe2, 0.0, 1.0) * e1;

        vec3  e2  = vec3(0.5 * r, -0.5 * a, -s * r);
        vec3  q2  = vec3(x - 0.5 * r, y - 0.33333333 * h, z - s * r);
        vec3  d2  = q2 - clamp(dot(q2, e2) * invLe2, 0.0, 1.0) * e2;

        vec3  e3  = vec3(r, 0.5 * a, 0.0);
        vec3  q3  = vec3(x, y + h, z);
        vec3  d3e = q3 - clamp(dot(q3, e3) * invLe2, 0.0, 1.0) * e3;

        float v1_2 = dot(d1, d1), v2_2 = dot(d2, d2), v3_2 = dot(d3e, d3e);
        vec3  ve   = (v1_2 <= v2_2) ? d1 : d2;
        float v2m  = min(v1_2, v2_2);
        if( v3_2 < v2m ) { ve = d3e; v2m = v3_2; }
        float invV = inversesqrt(max(v2m, 1e-30));
        d = v2m * invV;
        g = (v2m > 1e-30) ? (ve * invV) : vec3(a * invN, r * invN, 0.0);
    }

    if( swap ) g = vec3(g.x * c + g.z * s, -g.y, g.x * s - g.z * c);
    g.z *= sz1;
    if( r0 )
    {
        float gd = -s * g.x + c * g.z;
        g.x += 2.0 * gd * s;
        g.z -= 2.0 * gd * c;
    }
    g.z *= sz0;
    return vec4(d, g);
}

#endif // SDF_GEMS_OP_SHEAR_GLSL
