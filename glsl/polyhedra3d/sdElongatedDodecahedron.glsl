// ============================================================================
// sdElongatedDodecahedron — Exact 0-Loop, 0-Trig, 1-Sqrt 3D Elongated Dodecahedron (Fedorov Space-Filling Parallelohedron)
// @author Discovered and derived using Gemini 4 Argon & Gemini 3.8 Flash (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, exact 3D Euclidean distance & gradient)
// ============================================================================

#ifndef SDF_GEMS_SD_ELONGATED_DODECAHEDRON_GLSL
#define SDF_GEMS_SD_ELONGATED_DODECAHEDRON_GLSL

#include "../common/common.glsl"

vec4 sdgElongatedDodecahedron( in vec3 p, in float s, in float h )
{
    s = max(s, 1e-6);
    h = max(h, 0.0);
    vec3 sgn_p = vec3((p.x >= 0.0) ? 1.0 : -1.0,
                      (p.y >= 0.0) ? 1.0 : -1.0,
                      (p.z >= 0.0) ? 1.0 : -1.0);
    vec3 a  = abs(p);
    bool sw = (a.x < a.y);
    vec3 q  = sw ? vec3(a.y, a.x, a.z) : a;

    // Active half-spaces in D_4h chamber q.x >= q.y >= 0, q.z >= 0:
    // F1: q.x + q.y <= s       (vertical hexagonal face, normal = (1, 1, 0)/sqrt(2))
    // F2: q.x + q.z <= s + h   (sloped rhombic cap face, normal = (1, 0, 1)/sqrt(2))
    float m1 = q.x + q.y - s;
    float m2 = q.x + q.z - (s + h);
    const float INV_SQRT2 = 0.70710678118654752440;

    float d;
    vec3  gq;
    if( m1 <= 0.0 && m2 <= 0.0 )
    {
        if( m1 >= m2 )
        {
            d  = m1 * INV_SQRT2;
            gq = vec3(INV_SQRT2, INV_SQRT2, 0.0);
        }
        else
        {
            d  = m2 * INV_SQRT2;
            gq = vec3(INV_SQRT2, 0.0, INV_SQRT2);
        }
    }
    else
    {
        // 1. Clamped projection onto vertical hexagonal face F1 (x + y = s)
        float x1   = clamp(0.5 * (q.x - q.y + s), 0.5 * s, s);
        vec3  c1   = vec3(x1, s - x1, min(q.z, s + h - x1));
        vec3  df1  = q - c1;
        float best_d2 = dot(df1, df1);
        vec3  best_df = df1;

        // 2. Clamped projection onto sloped rhombic cap face F2 (x + z = s + h)
        float x2  = clamp(0.5 * (q.x - q.z + s + h), 0.0, s);
        vec3  c2  = vec3(x2, min(q.y, min(x2, s - x2)), s + h - x2);
        vec3  df2 = q - c2;
        float d2_2 = dot(df2, df2);
        if( d2_2 < best_d2 ) { best_d2 = d2_2; best_df = df2; }

        // 3. Diagonal shared edge E12 = F1 cap F2: (s/2 + t, s/2 - t, s/2 + h - t), t in [0, s/2]
        float t12  = clamp(((q.x - 0.5 * s) - (q.y - 0.5 * s) - (q.z - 0.5 * s - h)) * (1.0 / 3.0), 0.0, 0.5 * s);
        vec3  c12  = vec3(0.5 * s + t12, 0.5 * s - t12, 0.5 * s + h - t12);
        vec3  df12 = q - c12;
        float d2_12 = dot(df12, df12);
        if( d2_12 < best_d2 ) { best_d2 = d2_12; best_df = df12; }

        // 4. Top ridge edge Etop = F2 cap sigma_xy(F2): (t, t, s + h - t), t in [0, s/2]
        float ttop  = clamp((q.x + q.y - (q.z - s - h)) * (1.0 / 3.0), 0.0, 0.5 * s);
        vec3  ctop  = vec3(ttop, ttop, s + h - ttop);
        vec3  dftop = q - ctop;
        float d2_top = dot(dftop, dftop);
        if( d2_top < best_d2 ) { best_d2 = d2_top; best_df = dftop; }

        d  = sqrt(max(best_d2, 0.0));
        gq = (d > 1e-12) ? (best_df / d)
                         : ((m1 >= m2) ? vec3(INV_SQRT2, INV_SQRT2, 0.0) : vec3(INV_SQRT2, 0.0, INV_SQRT2));
    }

    vec3 g = sw ? vec3(gq.y, gq.x, gq.z) : gq;
    return vec4(d, sgn_p * g);
}

float sdElongatedDodecahedron( in vec3 p, in float s, in float h )
{
    return sdgElongatedDodecahedron(p, s, h).x;
}

#endif // SDF_GEMS_SD_ELONGATED_DODECAHEDRON_GLSL
