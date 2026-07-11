"""imgutil.py — trim the white border around a PNG so it embeds tightly in a notebook
(WISP's native waveform plots carry a large white top/side margin). Preserves the figure
content exactly; only crops surrounding white."""
import os
from PIL import Image, ImageChops


def trim(path, pad=12, out=None):
    """Return path to a whitespace-trimmed copy of the PNG (suffix _trim)."""
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    bbox = ImageChops.difference(im, bg).getbbox()
    if bbox is None:
        return path
    l, t, r, b = bbox
    l, t = max(0, l - pad), max(0, t - pad)
    r, b = min(im.size[0], r + pad), min(im.size[1], b + pad)
    out = out or os.path.splitext(path)[0] + "_trim.png"
    im.crop((l, t, r, b)).save(out)
    return out
