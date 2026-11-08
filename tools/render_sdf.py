"""
Repeatable 3D & 2D SDF rendering utilities for `sdf-gems`.

Provides:
  - render_3d_sdf(sdf_fn, ...)          : Raymarched 3D surface preview with ground-plane isolines
  - render_slice_3d(fn_3d, ...)         : 2D cross-section slice of a 3D SDF along any plane (u, v)
  - render_2d_sdf(fn_2d, ...)           : 2D SDF isoline & zero-set preview
  - render_eikonal_error_map(fn, ...)   : Visualizes Eikonal deviation (||grad d|| - 1)
  - render_error_heatmap(fn, gt, ...)   : log10 |fn - gt| slice heatmap with colorbar (preferred before/after view)
  - combine_panels(images, titles=None) : Horizontally concatenates PIL images with optional titles

Run via UV:
  uv run tools/render_sdf.py
"""

import numpy as np
from PIL import Image, ImageDraw


def smoothstep(e0, e1, x):
    t = np.clip((x - e0) / (e1 - e0 + 1e-12), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def render_3d_sdf(
    sdf_fn,
    res=(360, 270),
    cam_pos=(1.75, 1.35, 1.85),
    look_at=(0.0, -0.05, 0.0),
    max_steps=80,
    plane_y=-0.68,
    label=None,
):
    """
    Raymarches a 3D SDF `sdf_fn(p)` (where `p` has shape `(N, 3)`) and shades
    both the 3D zero-isosurface and a horizontal ground plane at `y = plane_y`
    displaying distance isolines. Returns a `PIL.Image`.
    """
    W, H = res if isinstance(res, tuple) else (res, res)
    ro = np.array(cam_pos, dtype=np.float64)
    ta = np.array(look_at, dtype=np.float64)
    fw = ta - ro
    fw /= np.linalg.norm(fw)
    rt = np.cross(fw, np.array([0.0, 1.0, 0.0]))
    rt /= np.linalg.norm(rt)
    up = np.cross(rt, fw)

    xs = (2.0 * (np.arange(W) + 0.5) / W - 1.0) * (W / H)
    ys = 1.0 - 2.0 * (np.arange(H) + 0.5) / H
    uu, vv = np.meshgrid(xs, ys)
    rd = uu[..., None] * rt + vv[..., None] * up + 1.85 * fw
    rd /= np.linalg.norm(rd, axis=-1, keepdims=True)
    rd_flat = rd.reshape(-1, 3)
    N = len(rd_flat)

    t = np.full(N, 0.1)
    hit_obj = np.zeros(N, dtype=bool)
    for _ in range(max_steps):
        p = ro[None, :] + t[:, None] * rd_flat
        d = sdf_fn(p)
        hit = d < 0.001
        hit_obj |= hit
        t = np.where(hit_obj | (t > 6.0), t, t + d)

    t_plane = np.where(rd_flat[:, 1] < -1e-4, (plane_y - ro[1]) / rd_flat[:, 1], 1e9)
    hit_plane = (t_plane > 0.0) & (t_plane < 8.0) & ((~hit_obj) | (t_plane < t))

    col = np.tile(np.array([0.32, 0.36, 0.42])[None, :], (N, 1))

    if np.any(hit_plane):
        idx = np.where(hit_plane)[0]
        p_pl = ro[None, :] + t_plane[idx, None] * rd_flat[idx]
        d_pl = sdf_fn(p_pl)
        c_pl = np.array([0.85, 0.65, 0.45])[None, :] * (0.85 + 0.15 * np.cos(70.0 * d_pl)[:, None])
        c_pl *= 1.0 - 0.45 * np.exp(-3.0 * np.maximum(d_pl, 0.0))[:, None]
        fade = np.exp(-0.08 * t_plane[idx] * t_plane[idx])[:, None]
        col[idx] = c_pl * fade + np.array([0.32, 0.36, 0.42])[None, :] * (1.0 - fade)

    hit_only_obj = hit_obj & (~hit_plane)
    if np.any(hit_only_obj):
        idx = np.where(hit_only_obj)[0]
        p_ob = ro[None, :] + t[idx, None] * rd_flat[idx]
        eps = 0.0015
        nx = sdf_fn(p_ob + [eps, 0, 0]) - sdf_fn(p_ob - [eps, 0, 0])
        ny = sdf_fn(p_ob + [0, eps, 0]) - sdf_fn(p_ob - [0, eps, 0])
        nz = sdf_fn(p_ob + [0, 0, eps]) - sdf_fn(p_ob - [0, 0, eps])
        n = np.stack([nx, ny, nz], axis=-1)
        n /= np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), 1e-8)

        lig = np.array([0.55, 0.72, 0.42])
        lig /= np.linalg.norm(lig)
        dif = np.clip(n @ lig, 0.0, 1.0)
        amb = np.clip(0.5 + 0.5 * n[:, 1], 0.0, 1.0)
        bac = np.clip(n @ np.array([-lig[0], -0.2, -lig[2]]), 0.0, 1.0)

        albedo = 0.66 + 0.20 * n
        c_ob = albedo * (
            0.62 * dif[:, None] * np.array([1.00, 0.95, 0.90])[None, :]
            + 0.34 * amb[:, None] * np.array([0.75, 0.84, 0.94])[None, :]
            + 0.18 * bac[:, None] * np.array([0.65, 0.62, 0.68])[None, :]
        )
        col[idx] = np.clip(c_ob, 0.0, 1.0)

    img = (np.clip(col.reshape(H, W, 3), 0.0, 1.0) * 255).astype(np.uint8)
    pil_img = Image.fromarray(img)
    if label:
        draw = ImageDraw.Draw(pil_img, "RGBA")
        text_w = len(label) * 6 + 14
        draw.rounded_rectangle([8, 8, 8 + text_w, 26], radius=4, fill=(15, 20, 28, 195))
        draw.text((15, 11), label, fill=(240, 244, 248, 255))
    return pil_img


def render_slice_3d(
    fn_3d,
    axis_u=(1.0, 0.0, 0.0),
    axis_v=(0.0, 1.0, 0.0),
    origin=(0.0, 0.0, 0.0),
    bounds=(-1.15, 1.15),
    res=320,
    freq=110.0,
    label=None,
):
    """
    Renders a 2D planar slice of a 3D SDF `fn_3d(p)` spanned by `(axis_u, axis_v)`
    around `origin`, using Inigo Quilez's orange/blue isoline palette and a crisp
    white zero-level contour. Returns a `PIL.Image`.
    """
    au = np.array(axis_u, dtype=np.float64)
    au /= np.linalg.norm(au)
    av = np.array(axis_v, dtype=np.float64)
    av /= np.linalg.norm(av)
    orig = np.array(origin, dtype=np.float64)

    lin = np.linspace(bounds[0], bounds[1], res)
    uu, vv = np.meshgrid(lin, -lin)
    p = orig[None, None, :] + uu[..., None] * au[None, None, :] + vv[..., None] * av[None, None, :]
    d = fn_3d(p.reshape(-1, 3)).reshape(res, res)

    col = np.where(
        d[..., None] > 0.0,
        np.array([0.90, 0.60, 0.30])[None, None, :],
        np.array([0.40, 0.70, 0.85])[None, None, :],
    )
    col *= 1.0 - np.exp(-6.0 * np.abs(d))[..., None]
    col *= 0.8 + 0.2 * np.cos(freq * d)[..., None]

    px = (bounds[1] - bounds[0]) / res * 1.8
    edge = 1.0 - smoothstep(0.0, px, np.abs(d))[..., None]
    col = col * (1.0 - edge) + edge

    img_rgb = np.clip(col * 255.0, 0, 255).astype(np.uint8)
    pil_img = Image.fromarray(img_rgb)
    if label:
        draw = ImageDraw.Draw(pil_img, "RGBA")
        text_w = len(label) * 6 + 14
        draw.rounded_rectangle([8, 8, 8 + text_w, 26], radius=4, fill=(15, 20, 28, 195))
        draw.text((15, 11), label, fill=(240, 244, 248, 255))
    return pil_img


def render_2d_sdf(fn_2d, bounds=(-1.15, 1.15), res=320, freq=110.0, label=None):
    """
    Renders a 2D SDF `fn_2d(p)` (where `p` has shape `(N, 2)` or `(..., 2)`).
    """
    return render_slice_3d(
        lambda p3: fn_2d(p3[..., :2]),
        axis_u=(1.0, 0.0, 0.0),
        axis_v=(0.0, 1.0, 0.0),
        bounds=bounds,
        res=res,
        freq=freq,
        label=label,
    )


def render_eikonal_error_map(fn_2d, bounds=(-1.15, 1.15), res=320, label=None):
    """
    Renders the signed Eikonal deviation `||grad d|| - 1.0` of a 2D SDF `fn_2d`:
      - Dark emerald: ||grad d|| == 1.0 (exact Euclidean SDF)
      - Cyan/Blue:    ||grad d|| <  1.0 (underestimation / lower bound)
      - Red/Yellow:   ||grad d|| >  1.0 (overestimation / raymarching overstep)
    """
    lin = np.linspace(bounds[0], bounds[1], res)
    xx, yy = np.meshgrid(lin, -lin)
    p = np.stack([xx, yy], axis=-1)
    d = fn_2d(p.reshape(-1, 2)).reshape(res, res)

    eps = (bounds[1] - bounds[0]) / res
    dx = (fn_2d((p + [eps, 0.0]).reshape(-1, 2)) - fn_2d((p - [eps, 0.0]).reshape(-1, 2))).reshape(res, res) / (2.0 * eps)
    dy = (fn_2d((p + [0.0, eps]).reshape(-1, 2)) - fn_2d((p - [0.0, eps]).reshape(-1, 2))).reshape(res, res) / (2.0 * eps)
    err = np.hypot(dx, dy) - 1.0

    col = np.full((res, res, 3), [0.08, 0.10, 0.12], dtype=np.float64)
    col += np.clip(-err / 0.5, 0.0, 1.0)[..., None] * np.array([0.0, 0.55, 0.90])
    col += np.clip(err / 0.3, 0.0, 1.0)[..., None] * np.array([1.0, 0.15, 0.05])
    exact = np.exp(-((err / 0.015) ** 2))
    col = col * (1.0 - 0.3 * exact[..., None]) + exact[..., None] * np.array([0.12, 0.28, 0.18])
    zero_line = np.exp(-((d / 0.012) ** 2))
    col = col * (1.0 - zero_line[..., None]) + zero_line[..., None]

    pil_img = Image.fromarray(np.clip(col * 255.0, 0, 255).astype(np.uint8))
    if label:
        draw = ImageDraw.Draw(pil_img, "RGBA")
        text_w = len(label) * 6 + 14
        draw.rounded_rectangle([8, 8, 8 + text_w, 26], radius=4, fill=(15, 20, 28, 195))
        draw.text((15, 11), label, fill=(240, 244, 248, 255))
    return pil_img


def render_error_heatmap(
    fn_3d,
    gt_3d,
    axis_u=(1.0, 0.0, 0.0),
    axis_v=(0.0, 1.0, 0.0),
    origin=(0.0, 0.0, 0.0),
    bounds=(-1.15, 1.15),
    res=320,
    log_range=(-8.0, -0.5),
    colorbar=True,
    label=None,
):
    """
    Renders log10 |fn - gt| on a planar slice (magma colormap) with the true zero set in cyan.
    Use this for before/after comparisons: unlike an eikonal map, it also exposes SDFs whose
    gradient has unit length but whose value is wrong (e.g. max() of exact distances).
    Returns a `PIL.Image` (with a labelled colorbar strip underneath when `colorbar=True`).
    """
    from matplotlib import cm

    au = np.array(axis_u, dtype=np.float64)
    au /= np.linalg.norm(au)
    av = np.array(axis_v, dtype=np.float64)
    av /= np.linalg.norm(av)
    lin = np.linspace(bounds[0], bounds[1], res)
    uu, vv = np.meshgrid(lin, -lin)
    p = (np.array(origin, dtype=np.float64)[None, None, :] + uu[..., None] * au + vv[..., None] * av).reshape(-1, 3)
    gt = np.asarray(gt_3d(p), dtype=np.float64).reshape(res, res)
    err = np.abs(np.asarray(fn_3d(p), dtype=np.float64).reshape(res, res) - gt)
    lo, hi = log_range
    t = np.clip((np.log10(np.maximum(err, 10.0**lo)) - lo) / (hi - lo), 0.0, 1.0)
    col = cm.magma(t)[..., :3]
    px = (bounds[1] - bounds[0]) / res * 1.6
    z = np.exp(-((gt / px) ** 2))[..., None]
    col = col * (1.0 - z) + z * np.array([0.3, 0.9, 1.0])
    img = Image.fromarray((np.clip(col, 0, 1) * 255).astype(np.uint8))
    draw = ImageDraw.Draw(img, "RGBA")
    txt = (label + "  " if label else "") + f"max {err.max():.1e}"
    draw.rounded_rectangle([8, 8, 8 + len(txt) * 6 + 14, 26], radius=4, fill=(15, 20, 28, 195))
    draw.text((15, 11), txt, fill=(240, 244, 248, 255))
    if not colorbar:
        return img
    bar = Image.new("RGB", (res, 30), (20, 22, 28))
    d = ImageDraw.Draw(bar)
    x0, x1 = 8, res - 8
    for x in range(x0, x1):
        c = cm.magma((x - x0) / (x1 - x0))[:3]
        d.line([(x, 4), (x, 13)], fill=tuple(int(255 * v) for v in c))
    for e in range(int(np.ceil(lo)), int(np.floor(hi)) + 1, 2):
        x = x0 + (e - lo) / (hi - lo) * (x1 - x0)
        d.text((x - 10, 16), f"1e{e}", fill=(230, 230, 230))
    out = Image.new("RGB", (res, res + 30))
    out.paste(img, (0, 0))
    out.paste(bar, (0, res))
    return out


def combine_panels(images):
    """Horizontally concatenates a sequence of `PIL.Image` panels."""
    widths, heights = zip(*(im.size for im in images))
    combined = Image.new("RGB", (sum(widths), max(heights)))
    x_offset = 0
    for im in images:
        combined.paste(im, (x_offset, 0))
        x_offset += im.size[0]
    return combined


if __name__ == "__main__":
    print("sdf-gems rendering toolkit ready. Import functions from `tools.render_sdf`.")
