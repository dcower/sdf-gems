#ifndef SDF_GEMS_SD_LIMACON_2D_GLSL
#define SDF_GEMS_SD_LIMACON_2D_GLSL

#include "../common/common.glsl"

vec3 sdgLimacon2D( in vec2 p, in float a, in float b )
{
    a = max(abs(a), 1e-6);
    b = max(abs(b), 0.0);
    float sy = (p.y >= 0.0) ? 1.0 : -1.0;
    float px = p.x, ay = abs(p.y);
    const float PI = 3.14159265358979323846;

    // Symmetry-axis endpoints theta = 0: (b + a, 0) and theta = pi: (a - b, 0)
    float best_th = 0.0;
    vec2  best_q  = vec2(b + a, 0.0);
    float best_d2 = (px - best_q.x) * (px - best_q.x) + ay * ay;

    float qx_pi = a - b;
    float d2_pi = (px - qx_pi) * (px - qx_pi) + ay * ay;
    if( d2_pi < best_d2 )
    {
        best_d2 = d2_pi;
        best_th = PI;
        best_q  = vec2(qx_pi, 0.0);
    }

    // Vectorized 4-seed SIMD Newton solver on [0, pi] (covers outer + inner loops)
    vec4 th = vec4(0.22, 1.00, 2.20, 2.85);
    for( int i = 0; i < 7; i++ )
    {
        vec4 c    = cos(th), s = sin(th);
        vec4 r    = b + a * c;
        vec4 qx   = r * c, qy = r * s;
        vec4 dx   = qx - px;
        vec4 sy_c = step(vec4(0.0), qy) * 2.0 - 1.0;
        vec4 dy   = qy - ay * sy_c;
        vec4 c2s2 = c * c - s * s, sc = s * c;
        vec4 qpx  = -b * s - (2.0 * a) * sc;
        vec4 qpy  =  b * c + a * c2s2;
        vec4 qppx = -b * c - (2.0 * a) * c2s2;
        vec4 qppy = -b * s - (4.0 * a) * sc;
        vec4 g    = dx * qpx + dy * qpy;
        vec4 gp   = qpx * qpx + qpy * qpy + dx * qppx + dy * qppy;
        vec4 step = vec4(
            (gp.x > 0.0) ? (g.x / max(gp.x, 1e-12)) : ((g.x >= 0.0 ? 1.0 : -1.0) * 0.25),
            (gp.y > 0.0) ? (g.y / max(gp.y, 1e-12)) : ((g.y >= 0.0 ? 1.0 : -1.0) * 0.25),
            (gp.z > 0.0) ? (g.z / max(gp.z, 1e-12)) : ((g.z >= 0.0 ? 1.0 : -1.0) * 0.25),
            (gp.w > 0.0) ? (g.w / max(gp.w, 1e-12)) : ((g.w >= 0.0 ? 1.0 : -1.0) * 0.25)
        );
        th = clamp(th - clamp(step, vec4(-0.25), vec4(0.25)), vec4(0.0), vec4(PI));
    }

    for( int j = 0; j < 4; j++ )
    {
        float tj = th[j];
        float c  = cos(tj), s = sin(tj);
        float r  = b + a * c;
        vec2  q  = vec2(r * c, abs(r * s));
        vec2  dv = vec2(px, ay) - q;
        float d2 = dot(dv, dv);
        if( d2 < best_d2 )
        {
            best_d2 = d2;
            best_th = tj;
            best_q  = q;
        }
    }

    vec2  diff = vec2(px, ay) - best_q;
    float d    = sqrt(max(best_d2, 0.0));

    // Exact 0-trig, 0-sqrt factored quartic algebraic sign test:
    // F(x, y) = (x^2 + y^2 - a*x)^2 - b^2*(x^2 + y^2)
    float r2_p = px * px + ay * ay;
    float u_p  = r2_p - a * px;
    float F    = u_p * u_p - (b * b) * r2_p;
    float sgn  = (F < 0.0 || (b > a && r2_p < 1e-12)) ? -1.0 : 1.0;

    vec2 g;
    if( d > 1e-6 * (a + b) )
    {
        g = sgn * normalize(vec2(diff.x, sy * diff.y));
    }
    else
    {
        float c0  = cos(best_th), s0 = sin(best_th);
        float r0  = b + a * c0;
        vec2  qp  = vec2(-b * s0 - 2.0 * a * s0 * c0, (r0 >= 0.0 ? 1.0 : -1.0) * (b * c0 + a * (c0 * c0 - s0 * s0)));
        g = normalize(vec2(qp.y, -sy * qp.x));
    }
    return vec3(sgn * d, g);
}

float sdLimacon2D( in vec2 p, in float a, in float b )
{
    return sdgLimacon2D(p, a, b).x;
}

vec3 sdgTrisectrix2D( in vec2 p, in float a )
{
    return sdgLimacon2D(p, a, 0.5 * a);
}

float sdTrisectrix2D( in vec2 p, in float a )
{
    return sdgLimacon2D(p, a, 0.5 * a).x;
}

vec4 sdgLimacon3D( in vec3 p, in float a, in float b )
{
    float rho = length(p.xz);
    vec2  er  = (rho > 1e-12) ? (p.xz / rho) : vec2(1.0, 0.0);
    vec3  dg  = sdgLimacon2D(vec2(p.y, rho), a, b);
    return vec4(dg.x, dg.z * er.x, dg.y, dg.z * er.y);
}

float sdLimacon3D( in vec3 p, in float a, in float b )
{
    return sdgLimacon2D(vec2(p.y, length(p.xz)), a, b).x;
}

#endif // SDF_GEMS_SD_LIMACON_2D_GLSL