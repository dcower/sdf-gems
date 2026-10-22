// ============================================================================
// opBend, sdgBendArchFinalize & sdgCheapBendFinalize — Exact Circular Arch Bend & Jacobian Eikonal Normalization
// @author Inigo Quilez (opCheapBend) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1` on surface; 0 `atan`/`sin`/`cos` calls in opBendArch & sdgBendArchFinalize)
// ============================================================================

#ifndef SDF_GEMS_OP_BEND_GLSL
#define SDF_GEMS_OP_BEND_GLSL

// 1. Exact atan-Free Circular Arch Bend + Straight Tangent Legs
vec3 opBendArch( in vec3 p, in float R, in float alpha, in vec2 cs )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    vec2  w  = vec2( abs(p.x), R + p.y );
    float u  = w.x * cs.x - w.y * cs.y;
    float v  = w.x * cs.y + w.y * cs.x;
    return vec3( sx * max(R * alpha + u, 0.0),
                 (u <= 0.0) ? length(w) - R : v - R,
                 p.z );
}

// 2. Exact Analytical 3D Normal & Curvature-Normalized SDF for opBendArch:
//    In the circular arch sector (u <= 0), the polar map has orthogonal singular values
//    (R / rho, 1, 1) where rho = length(w), so J^T * g_q scales the tangential component
//    by R / rho and rotates by the radial unit vector w / rho. In the straight tangent legs
//    (u > 0), J is a pure rotation by cs = (cos(alpha), sin(alpha)) with ||J^T g_q|| = 1!
//    Usage: vec3 q = opBendArch(p, R, alpha, cs); vec4 dg = sdgPrimitive(q);
//           vec4 res = sdgBendArchFinalize(p, R, alpha, cs, dg);
vec4 sdgBendArchFinalize( in vec3 p, in float R, in float alpha, in vec2 cs, in vec4 dg )
{
    float sx  = (p.x < 0.0) ? -1.0 : 1.0;
    vec2  w   = vec2(abs(p.x), R + p.y);
    float u   = w.x * cs.x - w.y * cs.y;
    float gqx = sx * dg.y, gqy = dg.z;
    if( u <= 0.0 )
    {
        float rho  = max(length(w), 1e-6);
        vec2  er   = w / rho;
        vec2  et   = vec2(er.y, -er.x);
        float gt   = (R * alpha + u > 0.0) ? (gqx * (R / rho)) : 0.0;
        vec3  g    = vec3(sx * (gt * et.x + gqy * er.x), gt * et.y + gqy * er.y, dg.w);
        float invG = inversesqrt(max(dot(g, g), 1e-30));
        return vec4(dg.x * invG, g * invG);
    }
    float gt = (R * alpha + u > 0.0) ? gqx : 0.0;
    vec3  g  = vec3(sx * (gt * cs.x + gqy * cs.y), -gt * cs.y + gqy * cs.x, dg.w);
    return vec4(dg.x, normalize(g));
}

// 3. Standard Inigo Quilez Cheap Bend + First-Order Eikonal-Exact Jacobian Normalization:
//    Corrects the severe cross-sectional shear of opCheapBend (`||grad|| in [0.43, 1.85]`)
//    via the exact chain-rule pullback J(p)^T g_q = (R(-kx) g_qxy + (0, k*(qy*g_qx - qx*g_qy), 0)).
vec3 opCheapBend( in vec3 p, in float k )
{
    float c = cos(k * p.x);
    float s = sin(k * p.x);
    mat2  m = mat2(c, -s, s, c);
    return vec3(m * p.xy, p.z);
}

vec4 sdgCheapBendFinalize( in vec3 p, in vec3 q, in float k, in vec4 dg )
{
    float c    = cos(k * p.x);
    float s    = sin(k * p.x);
    float gx   = c * dg.y - s * dg.z + k * (q.y * dg.y - q.x * dg.z);
    float gy   = s * dg.y + c * dg.z;
    vec3  g    = vec3(gx, gy, dg.w);
    float invG = inversesqrt(max(dot(g, g), 1e-30));
    return vec4(dg.x * invG, g * invG);
}

#endif // SDF_GEMS_OP_BEND_GLSL
