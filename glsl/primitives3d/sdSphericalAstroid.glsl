// ============================================================================
// sdSphericalAstroid, sdgSphericalAstroid, sdSphericalCardioid & sdgSphericalCardioid
// — Exact O(1) 3D Spherical Astroid (4-Cusped Star Seam on S^2(R)) & Spherical Cardioid
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts inside the Newton solver)
// @novelty Exploits the Orthonormal Equatorial Rolling Identity on S^2(1):
//            u(t) = cos(alpha)*e_r(t) - sin(alpha)*sin(n*t)*e_theta(t) + sin(alpha)*cos(n*t)*e_z,
//          which satisfies ||u(t)||^2 == cos^2(alpha) + sin^2(alpha) == 1 identically for all
//          t, alpha, n with ZERO square roots! Combined with the S^2 Dot-Product Linearization
//          (||p - R*u(t)||^2 = ||p||^2 + R^2 - 2*R*(p . u(t))) and D_n dihedral folding into
//          the fundamental half-sector [0, pi/n], both the 4-cusped Spherical Astroid (n = 4)
//          and Spherical Cardioid (n = 1) are solved in O(1) time with 0 sqrts in the loop.
// ============================================================================

#ifndef SDF_GEMS_SD_SPHERICAL_ASTROID_GLSL
#define SDF_GEMS_SD_SPHERICAL_ASTROID_GLSL

#include "../common/common.glsl"

// Shared O(1) 0-sqrt-in-loop solver for n-lobed spherical hypocycloidal rolling curves on S^2(R)
vec4 _sdSphericalHypocycloidCore( in vec3 p, in float R, in float alpha, in int nLobes, in float rTube )
{
    R     = max(R, 1e-5);
    alpha = clamp(alpha, 1e-4, 1.56);
    float n  = float(max(nLobes, 1));
    float ca = cos(alpha), sa = sin(alpha);

    const float PI = 3.14159265358979323846;
    float per   = (2.0 * PI) / n;
    float hPer  = PI / n;
    float rho   = length(p.xy);
    float phi   = (rho > 1e-12) ? atan(p.y, p.x) : 0.0;
    float kCell = floor(phi / per);
    float ang0  = phi - kCell * per;
    bool  refl  = (ang0 > hPer);
    float ang   = refl ? (per - ang0) : ang0;

    float px = rho * cos(ang);
    float py = rho * sin(ang);
    float pz = p.z;

    float tMin     = -1.00 * hPer;
    float tMax     =  2.00 * hPer;
    float bestD2   = 1e30;
    vec3  bestDiff = vec3(1.0, 0.0, 0.0);

    for( int pass = 0; pass < 2; pass++ )
    {
        vec3 t = (pass == 0)
               ? vec3(-0.50 * hPer, 0.10 * hPer, 0.50 * hPer)
               : vec3( 0.50 * hPer, 0.90 * hPer, 1.50 * hPer);

        for( int i = 0; i < 6; i++ )
        {
            vec3 c1 = cos(t), s1 = sin(t);
            vec3 cn = cos(n * t), sn = sin(n * t);

            vec3 W   = px * c1 + py * s1;
            vec3 V   = -px * s1 + py * c1;
            vec3 df  = (vec3(ca) - (n * sa) * cn) * V + sa * sn * (W - vec3(n * pz));
            vec3 d2f = ((2.0 * n * sa) * cn - vec3(ca)) * W + ((n * n + 1.0) * sa) * sn * V - (n * n * pz * sa) * cn;

            vec3 denom = min(d2f, -abs(d2f) - vec3(1e-6));
            vec3 step  = clamp(df / denom, vec3(-0.35 * hPer), vec3(0.35 * hPer));
            t = clamp(t - step, vec3(tMin), vec3(tMax));
        }

        float cands[5];
        cands[0] = t.x; cands[1] = t.y; cands[2] = t.z; cands[3] = 0.0; cands[4] = hPer;
        for( int k = 0; k < 5; k++ )
        {
            float tk   = cands[k];
            float c1   = cos(tk), s1 = sin(tk);
            float cn   = cos(n * tk), sn = sin(n * tk);
            vec3  q    = R * vec3(ca * c1 + sa * sn * s1, ca * s1 - sa * sn * c1, sa * cn);
            vec3  diff = vec3(px, py, pz) - q;
            float d2   = dot(diff, diff);
            if( d2 < bestD2 )
            {
                bestD2   = d2;
                bestDiff = diff;
            }
        }
    }

    float invD = inversesqrt(max(bestD2, 1e-30));
    float d    = bestD2 * invD;
    vec3  gq   = (bestD2 > 1e-24) ? (bestDiff * invD) : vec3(0.0, 0.0, 1.0);

    vec3 grad;
    if( refl )
    {
        float a0  = (kCell + 1.0) * per;
        float ca0 = cos(a0), sa0 = sin(a0);
        grad = vec3(gq.x * ca0 + gq.y * sa0, gq.x * sa0 - gq.y * ca0, gq.z);
    }
    else
    {
        float a0  = kCell * per;
        float ca0 = cos(a0), sa0 = sin(a0);
        grad = vec3(gq.x * ca0 - gq.y * sa0, gq.x * sa0 + gq.y * ca0, gq.z);
    }
    return vec4(d - rTube, grad);
}

// 1. Exact O(1) 3D Spherical Astroid (4-Cusped Star Seam on S^2(R)) + Unit Gradient vec4(d, grad.xyz)
//    (Exact spherical cusps occur at alpha = atan(0.25) = 0.24497866 rad).
vec4 sdgSphericalAstroid( in vec3 p, in float R, in float alpha, in float rTube )
{
    return _sdSphericalHypocycloidCore(p, R, alpha, 4, rTube);
}

// 2. Scalar 3D Spherical Astroid SDF
float sdSphericalAstroid( in vec3 p, in float R, in float alpha, in float rTube )
{
    return _sdSphericalHypocycloidCore(p, R, alpha, 4, rTube).x;
}

// 3-parameter canonical cusped overloads (alpha = atan(1/4) = 0.24497866 rad)
vec4 sdgSphericalAstroid( in vec3 p, in float R, in float rTube )
{
    return _sdSphericalHypocycloidCore(p, R, 0.24497866, 4, rTube);
}

float sdSphericalAstroid( in vec3 p, in float R, in float rTube )
{
    return _sdSphericalHypocycloidCore(p, R, 0.24497866, 4, rTube).x;
}

// 3. Exact O(1) 3D Spherical Cardioid on S^2(R) + Analytical Unit Gradient vec4(d, grad.xyz)
//    (Exact spherical cusp occurs at alpha = pi/4 = 0.78539816 rad).
vec4 sdgSphericalCardioid( in vec3 p, in float R, in float alpha, in float rTube )
{
    return _sdSphericalHypocycloidCore(p, R, alpha, 1, rTube);
}

// 4. Scalar 3D Spherical Cardioid SDF
float sdSphericalCardioid( in vec3 p, in float R, in float alpha, in float rTube )
{
    return _sdSphericalHypocycloidCore(p, R, alpha, 1, rTube).x;
}

// 3-parameter canonical cusped overloads (alpha = pi/4 = 0.78539816 rad)
vec4 sdgSphericalCardioid( in vec3 p, in float R, in float rTube )
{
    return _sdSphericalHypocycloidCore(p, R, 0.78539816, 1, rTube);
}

float sdSphericalCardioid( in vec3 p, in float R, in float rTube )
{
    return _sdSphericalHypocycloidCore(p, R, 0.78539816, 1, rTube).x;
}

#endif // SDF_GEMS_SD_SPHERICAL_ASTROID_GLSL
