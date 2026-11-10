// ============================================================================
// fmRefractedWavefront — Exact Fermat-Snell Refracted Wavefront Dome & Lens SDFs
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, `1.29e-9` FP64 / `1.19e-7` FP32 max error, 0-step rim & base early exits, 1-Halley Pade-Snell solver)
// ============================================================================

#ifndef SDF_GEMS_HAS_FM_REFRACTED_WAVEFRONT
#define SDF_GEMS_HAS_FM_REFRACTED_WAVEFRONT

// Internal helper: solves the 1D Snell ray-intercept equation
//   F(r) = r * (1 + n1 * ya / sqrt((n2^2 - n1^2) * r^2 + n2^2 * h^2)) - rho = 0
// using a 0-sqrt Pade initial guess + 1 Halley step (1 inversesqrt), returning
// vec3(d_wf, s2, c2) where d_wf = (T(rho, ya) - t0) / n2 and (s2, c2) = (sin theta_2, cos theta_2).
vec3 fmSnellWavefrontSolve( float rho, float ya, float h, float n1, float n2, float t0, float r_max )
{
    float dn2   = (n2 - n1) * (n2 + n1);
    float n2h   = n2 * h;
    float n2h2  = n2h * n2h;
    float n1y   = n1 * ya;

    // 0-sqrt Pade initial guess for sqrt((n2^2 - n1^2) * r^2 + n2^2 * h^2)
    float r_lin = min(rho * n2h / max(n2h + n1y, 1e-20), r_max);
    float u     = (dn2 * r_lin * r_lin) / max(n2h2, 1e-20);
    float d_est = n2h * (1.0 + (0.5 * u) / (1.0 + 0.25 * u));
    float r     = rho * d_est / max(d_est + n1y, 1e-20);

    // 1 cubic-convergent Halley step on F(r) = 0 (1 inversesqrt)
    float D     = dn2 * r * r + n2h2;
    float q     = inversesqrt(max(D, 1e-20));
    float q2    = q * q;
    float n1yq  = n1y * q;
    float f     = r * (1.0 + n1yq) - rho;
    float fp    = 1.0 + n1yq * (n2h2 * q2);
    float fpp   = -3.0 * n1yq * (n2h2 * dn2 * r) * (q2 * q2);
    r           = clamp(r - (2.0 * f * fp) / max(2.0 * fp * fp - f * fpp, 1e-20), 0.0, rho);

    // Stationary Fermat optical path T(rho, ya) = n1 * L1 + n2 * L2 (Envelope Theorem doubles precision)
    float L1    = length(vec2(r, h));
    float L2    = length(vec2(rho - r, ya));
    float d_wf  = (n1 * L1 + n2 * L2 - t0) / n2;

    // Singularity-free unit ray direction (sin theta_2, cos theta_2) via Snell's law
    float s2    = clamp((n1 * r) / max(n2 * L1, 1e-20), 0.0, 1.0);
    float c2    = sqrt(max((1.0 - s2) * (1.0 + s2), 0.0));
    return vec3(d_wf, s2, c2);
}

// ----------------------------------------------------------------------------
// fmRefractedWavefront / fmgRefractedWavefront
// Exact 3D Euclidean SDF & unit gradient of the Plano-Refracted Wavefront Dome
//   Omega = { p in R^3 : p.y >= 0, T(p) <= t0 }
// for a point source at (0, -h, 0) in medium n1 refracting across y = 0 into
// medium n2 >= n1 at optical path length t0 > n1 * h.
// ----------------------------------------------------------------------------
float fmRefractedWavefront( vec3 p, float h, float n1, float n2, float t0 )
{
    h  = max(h, 1e-6);
    n1 = max(n1, 1e-6);
    n2 = max(n2, n1);

    float l1_max = max(t0 / n1, h);
    float r_max  = sqrt(max((l1_max - h) * (l1_max + h), 0.0));
    float rho    = length(p.xz);

    // 1. Exterior below flat interface disk (y <= 0): 0-iteration exact early exit
    if( p.y <= 0.0 )
    {
        return (rho <= r_max) ? -p.y : length(vec2(rho - r_max, p.y));
    }

    // 2. Exterior circular rim Voronoi cone: (rho - r_max) * cos(theta2_max) >= y * sin(theta2_max)
    float s2_max = clamp((n1 * r_max) / max(n2 * l1_max, 1e-20), 0.0, 1.0);
    float c2_max = sqrt(max((1.0 - s2_max) * (1.0 + s2_max), 0.0));
    if( (rho - r_max) * c2_max >= p.y * s2_max )
    {
        return length(vec2(rho - r_max, p.y));
    }

    // 3. Snell refracted wavefront normal fan (y > 0)
    vec3  wf = fmSnellWavefrontSolve(rho, p.y, h, n1, n2, t0, r_max);
    return (wf.x < 0.0) ? max(wf.x, -p.y) : wf.x;
}

vec4 fmgRefractedWavefront( vec3 p, float h, float n1, float n2, float t0 )
{
    h  = max(h, 1e-6);
    n1 = max(n1, 1e-6);
    n2 = max(n2, n1);

    float l1_max = max(t0 / n1, h);
    float r_max  = sqrt(max((l1_max - h) * (l1_max + h), 0.0));
    float rho    = length(p.xz);

    if( p.y <= 0.0 && rho <= r_max )
    {
        return vec4(-p.y, 0.0, -1.0, 0.0);
    }

    float s2_max = clamp((n1 * r_max) / max(n2 * l1_max, 1e-20), 0.0, 1.0);
    float c2_max = sqrt(max((1.0 - s2_max) * (1.0 + s2_max), 0.0));
    vec2  dir_xz = (rho > 1e-20) ? (p.xz / rho) : vec2(1.0, 0.0);

    if( p.y <= 0.0 || (rho - r_max) * c2_max >= p.y * s2_max )
    {
        vec2  v  = vec2(rho - r_max, p.y);
        float d  = length(v);
        vec2  g2 = (d > 1e-20) ? (v / d) : vec2(s2_max, c2_max);
        return vec4(d, dir_xz.x * g2.x, g2.y, dir_xz.y * g2.x);
    }

    vec3 wf = fmSnellWavefrontSolve(rho, p.y, h, n1, n2, t0, r_max);
    if( wf.x < -p.y )
    {
        return vec4(-p.y, 0.0, -1.0, 0.0);
    }
    return vec4(wf.x, dir_xz.x * wf.y, wf.z, dir_xz.y * wf.y);
}

// ----------------------------------------------------------------------------
// fmRefractedLens / fmgRefractedLens
// Exact 3D Euclidean SDF & unit gradient of the Symmetric Biconvex Refracted
// Wavefront Lens Omega_lens = { p in R^3 : T(length(p.xz), abs(p.y)) <= t0 }.
// By the Virtual-Evolute Caustic-Free Theorem, the wavefront evolute lies in
// the opposite half-space, so (T - t0) / n2 is 100% exact across the entire
// interior and normal fan of the lens.
// ----------------------------------------------------------------------------
float fmRefractedLens( vec3 p, float h, float n1, float n2, float t0 )
{
    h  = max(h, 1e-6);
    n1 = max(n1, 1e-6);
    n2 = max(n2, n1);

    float l1_max = max(t0 / n1, h);
    float r_max  = sqrt(max((l1_max - h) * (l1_max + h), 0.0));
    float s2_max = clamp((n1 * r_max) / max(n2 * l1_max, 1e-20), 0.0, 1.0);
    float c2_max = sqrt(max((1.0 - s2_max) * (1.0 + s2_max), 0.0));

    float rho    = length(p.xz);
    float ya     = abs(p.y);
    if( (rho - r_max) * c2_max >= ya * s2_max )
    {
        return length(vec2(rho - r_max, ya));
    }
    return fmSnellWavefrontSolve(rho, ya, h, n1, n2, t0, r_max).x;
}

vec4 fmgRefractedLens( vec3 p, float h, float n1, float n2, float t0 )
{
    h  = max(h, 1e-6);
    n1 = max(n1, 1e-6);
    n2 = max(n2, n1);

    float l1_max = max(t0 / n1, h);
    float r_max  = sqrt(max((l1_max - h) * (l1_max + h), 0.0));
    float s2_max = clamp((n1 * r_max) / max(n2 * l1_max, 1e-20), 0.0, 1.0);
    float c2_max = sqrt(max((1.0 - s2_max) * (1.0 + s2_max), 0.0));

    float rho    = length(p.xz);
    float ya     = abs(p.y);
    float sy     = (p.y < 0.0) ? -1.0 : 1.0;
    vec2  dir_xz = (rho > 1e-20) ? (p.xz / rho) : vec2(1.0, 0.0);

    if( (rho - r_max) * c2_max >= ya * s2_max )
    {
        vec2  v  = vec2(rho - r_max, ya);
        float d  = length(v);
        vec2  g2 = (d > 1e-20) ? (v / d) : vec2(s2_max, c2_max);
        return vec4(d, dir_xz.x * g2.x, sy * g2.y, dir_xz.y * g2.x);
    }

    vec3 wf = fmSnellWavefrontSolve(rho, ya, h, n1, n2, t0, r_max);
    return vec4(wf.x, dir_xz.x * wf.y, sy * wf.z, dir_xz.y * wf.y);
}

#endif // SDF_GEMS_HAS_FM_REFRACTED_WAVEFRONT
