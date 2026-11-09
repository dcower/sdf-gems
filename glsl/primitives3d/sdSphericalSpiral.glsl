#ifndef SDF_GEMS_SD_SPHERICAL_SPIRAL_GLSL
#define SDF_GEMS_SD_SPHERICAL_SPIRAL_GLSL

#include "../common/common.glsl"

vec4 sdgSphericalSpiral( in vec3 p, in float R, in float a, in float rTube )
{
    R = max(R, 1e-6);
    a = clamp(a, 1e-3, 2.0);
    const float PI         = 3.14159265358979323846;
    const float TWO_PI     = 6.28318530717958647692;
    const float INV_TWO_PI = 0.15915494309189533577;

    float rho    = length(p.xy);
    float phi0   = atan(p.y, p.x);
    float psi0   = atan(rho, p.z); // spherical colatitude in [0, pi]
    float th_max = PI / a;

    // Polar endpoints theta = 0 (North Pole) and theta = pi/a (South Pole)
    vec3  best_q  = vec3(0.0, 0.0, R);
    vec3  dv0     = p - best_q;
    float best_d2 = dot(dv0, dv0);

    vec3  q_south  = vec3(0.0, 0.0, -R);
    vec3  dv_south = p - q_south;
    float d2_south = dot(dv_south, dv_south);
    if( d2_south < best_d2 )
    {
        best_d2 = d2_south;
        best_q  = q_south;
    }

    // O(1) Colatitude Turn Bracketing on front, back, and polar meridians
    float k0      = round((psi0 / a - phi0) * INV_TWO_PI);
    float kpi     = round((psi0 / a - (phi0 + PI)) * INV_TWO_PI);
    float k_north = ceil((-phi0 - PI) * INV_TWO_PI);
    float k_south = floor((th_max - phi0 - PI) * INV_TWO_PI);

    vec4 thA = clamp(vec4(
        phi0 + TWO_PI * (k0 - 1.0),
        phi0 + TWO_PI * k0,
        phi0 + TWO_PI * (k0 + 1.0),
        phi0 + PI + TWO_PI * (kpi - 1.0)
    ), vec4(0.0), vec4(th_max));

    vec4 thB = clamp(vec4(
        phi0 + PI + TWO_PI * kpi,
        phi0 + PI + TWO_PI * (kpi + 1.0),
        phi0 + PI + TWO_PI * k_north,
        phi0 + PI + TWO_PI * k_south
    ), vec4(0.0), vec4(th_max));

    for( int batch = 0; batch < 2; batch++ )
    {
        vec4 th = (batch == 0) ? thA : thB;
        for( int it = 0; it < 6; it++ )
        {
            vec4 u   = a * th;
            vec4 v   = th - phi0;
            vec4 cu  = cos(u), su = sin(u);
            vec4 cv  = cos(v), sv = sin(v);
            vec4 Fp  = rho * (a * cu * cv - su * sv) - (a * p.z) * su;
            vec4 Fpp = -rho * ((1.0 + a * a) * su * cv + (2.0 * a) * cu * sv) - (a * a * p.z) * cu;
            vec4 step = vec4(
                (Fpp.x < -1e-12) ? (Fp.x / Fpp.x) : ((Fp.x >= 0.0 ? -1.0 : 1.0) * 0.3),
                (Fpp.y < -1e-12) ? (Fp.y / Fpp.y) : ((Fp.y >= 0.0 ? -1.0 : 1.0) * 0.3),
                (Fpp.z < -1e-12) ? (Fp.z / Fpp.z) : ((Fp.z >= 0.0 ? -1.0 : 1.0) * 0.3),
                (Fpp.w < -1e-12) ? (Fp.w / Fpp.w) : ((Fp.w >= 0.0 ? -1.0 : 1.0) * 0.3)
            );
            th = clamp(th - clamp(step, vec4(-0.6), vec4(0.6)), vec4(0.0), vec4(th_max));
        }
        for( int j = 0; j < 4; j++ )
        {
            float tj = th[j];
            float u  = a * tj;
            float su = sin(u), cu = cos(u);
            vec3  q  = R * vec3(su * cos(tj), su * sin(tj), cu);
            vec3  dv = p - q;
            float d2 = dot(dv, dv);
            if( d2 < best_d2 )
            {
                best_d2 = d2;
                best_q  = q;
            }
        }
    }

    float d = sqrt(max(best_d2, 1e-30));
    return vec4(d - rTube, (p - best_q) / d);
}

float sdSphericalSpiral( in vec3 p, in float R, in float a, in float rTube )
{
    return sdgSphericalSpiral(p, R, a, rTube).x;
}

#endif // SDF_GEMS_SD_SPHERICAL_SPIRAL_GLSL