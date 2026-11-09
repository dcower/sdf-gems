#ifndef SDF_GEMS_HAS_SD_OCTABOX
#define SDF_GEMS_HAS_SD_OCTABOX

#include "../common/common.glsl"

// ============================================================================
// sdOctaBox / sdgOctaBox — Exact 3D Octahedron-Rounded Box (Crystal-Beveled Octa-Box)
// @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt`)
// ============================================================================

// Helper: sort q = abs(p) - b into x >= y >= z and invert the permutation on a Weyl-chamber vector g.
void _sdgOctaBoxSort( vec3 p, vec3 b, out vec3 q, out vec3 w, out vec3 s0, out bvec3 sw )
{
    s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    q  = abs(p) - b;
    sw.x = (q.x < q.y); if( sw.x ) q.xy = q.yx;
    sw.y = (q.y < q.z); if( sw.y ) q.yz = q.zy;
    sw.z = (q.x < q.y); if( sw.z ) q.xy = q.yx;
    w  = max(q, 0.0);
}

vec3 _sdgOctaBoxUnsort( vec3 g, vec3 s0, bvec3 sw )
{
    if( sw.z ) g.xy = g.yx;
    if( sw.y ) g.yz = g.zy;
    if( sw.x ) g.xy = g.yx;
    return s0 * g;
}

float sdOctaBox( vec3 p, vec3 b, float s )
{
    vec3 q = abs(p) - b;
    vec3 w = max(q, 0.0);
    float sum_w = w.x + w.y + w.z;

    if (sum_w > s)
    {
        float m = sum_w - s;
        float z = min(w.x, min(w.y, w.z));
        float a = max(m - 3.0 * z, 0.0);
        float w_max = max(w.x, max(w.y, w.z));
        float w_mid = sum_w - w_max - z;
        float bb = max(w_max - w_mid - s, 0.0);
        return sqrt((2.0 * m * m + a * a + 3.0 * bb * bb) * (1.0 / 6.0));
    }
    
    float q_max = max(q.x, max(q.y, q.z));
    float q_min = min(q.x, min(q.y, q.z));
    float sum_q = q.x + q.y + q.z;
    
    float d1 = q_max - s;
    float d2 = (sum_q - q_min - s) * 0.7071067811865475;
    float d3 = (sum_q - s) * 0.5773502691896258;

    return max(d1, max(d2, d3));
}

vec4 sdgOctaBox( vec3 p, vec3 b, float s )
{
    vec3  q, w, s0;
    bvec3 sw;
    _sdgOctaBoxSort(p, b, q, w, s0, sw);

    float sum_w = w.x + w.y + w.z;
    
    if (sum_w > s)
    {
        float m  = sum_w - s;
        float a  = max(m - 3.0 * w.z, 0.0);
        float bb = max(w.x - w.y - s, 0.0);
        float d  = sqrt((2.0 * m * m + a * a + 3.0 * bb * bb) * (1.0 / 6.0));
        
        vec3 g = (d > 0.0) 
            ? (vec3(2.0 * m + a + 3.0 * bb, 2.0 * m + a - 3.0 * bb, 2.0 * m - 2.0 * a) / (6.0 * d))
            : vec3(0.5773502691896258);
            
        g *= step(0.0, q);
        return vec4(d, _sdgOctaBoxUnsort(g, s0, sw));
    }
    
    float sum_q = q.x + q.y + q.z;
    float d1 = q.x - s;
    float d2 = (sum_q - q.z - s) * 0.7071067811865475;
    float d3 = (sum_q - s) * 0.5773502691896258;
    
    float d = max(d1, max(d2, d3));
    
    vec3 g = (d1 >= d2 && d1 >= d3) ? vec3(1.0, 0.0, 0.0) :
             (d2 >= d3)             ? vec3(0.7071067811865475, 0.7071067811865475, 0.0) :
                                      vec3(0.5773502691896258);
    
    return vec4(d, _sdgOctaBoxUnsort(g, s0, sw));
}

// ============================================================================
// opExtrudeCone / sdgExtrudeCone — Exact Euclidean Conical / Chamfered Extrusion
// ============================================================================

float opExtrudeCone( float d2d, float z, float h, float r, float c )
{
    vec2 p = vec2(d2d, abs(z) - h);
    vec2 w = max(p, 0.0);
    
    float sq_rc = r * r + c * c;
    if (sq_rc <= 0.0) {
        return max(p.x, p.y);
    }
    
    if (w.x * c + w.y * r > r * c)
    {
        vec2 pa = w - vec2(r, 0.0);
        vec2 ba = vec2(-r, c);
        float h_seg = clamp(dot(pa, ba) / sq_rc, 0.0, 1.0);
        return length(pa - ba * h_seg);
    }
    
    float d1 = p.x - r;
    float d2 = p.y - c;
    float d3 = (p.x * c + p.y * r - r * c) / sqrt(sq_rc);
    return max(d1, max(d2, d3));
}

vec4 sdgExtrudeCone( vec3 dg2d, float z, float h, float r, float c )
{
    vec2 p = vec2(dg2d.x, abs(z) - h);
    vec2 w = max(p, 0.0);
    
    float sq_rc = r * r + c * c;
    if (sq_rc <= 0.0) {
        float d = max(p.x, p.y);
        vec2 g_p = (p.x >= p.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
        float sz0 = (z < 0.0) ? -1.0 : 1.0;
        return vec4(d, g_p.x * dg2d.y, g_p.x * dg2d.z, g_p.y * sz0);
    }
    
    float d;
    vec2 g_p;
    
    if (w.x * c + w.y * r > r * c)
    {
        vec2 pa = w - vec2(r, 0.0);
        vec2 ba = vec2(-r, c);
        float h_seg = clamp(dot(pa, ba) / sq_rc, 0.0, 1.0);
        vec2 v = pa - ba * h_seg;
        d = length(v);
        g_p = (d > 0.0) ? (v / d) : vec2(0.0, 1.0);
        g_p *= step(0.0, p);
    }
    else
    {
        float inv_len = inversesqrt(sq_rc);
        float d1 = p.x - r;
        float d2 = p.y - c;
        float d3 = (p.x * c + p.y * r - r * c) * inv_len;
        
        d = max(d1, max(d2, d3));
        g_p = (d1 >= d2 && d1 >= d3) ? vec2(1.0, 0.0) :
              (d2 >= d3)             ? vec2(0.0, 1.0) :
                                       vec2(c, r) * inv_len;
    }
    
    float sz = (z < 0.0) ? -1.0 : 1.0;
    return vec4(d, g_p.x * dg2d.y, g_p.x * dg2d.z, g_p.y * sz);
}

#endif
