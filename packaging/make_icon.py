"""Draw minimol's icon (a small ball-and-stick molecule) as PNG and Windows ICO.

Run:  python packaging/make_icon.py   (needs Pillow; only used when changing the icon)
Writes src/minimol/static/icon.png (window icon) and packaging/windows/minimol.ico.
"""

from pathlib import Path

from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
SIZE = 512


def draw(size: int = SIZE) -> Image.Image:
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    g = ImageDraw.Draw(img)
    s = size / 512
    # Rounded tile.
    g.rounded_rectangle([16 * s, 16 * s, 496 * s, 496 * s], radius=96 * s, fill=(37, 99, 235, 255))
    atoms = {"c1": (190, 300), "c2": (330, 220), "o": (380, 360), "h": (120, 175)}
    for a, b in (("c1", "c2"), ("c2", "o"), ("c1", "h")):
        (x1, y1), (x2, y2) = atoms[a], atoms[b]
        g.line([x1 * s, y1 * s, x2 * s, y2 * s], fill=(255, 255, 255, 255), width=round(34 * s))
    style = {"c1": (52, (70, 70, 70)), "c2": (52, (70, 70, 70)), "o": (48, (232, 49, 42)), "h": (36, (242, 242, 242))}
    for name, (x, y) in atoms.items():
        r, colour = style[name]
        g.ellipse([(x - r - 7) * s, (y - r - 7) * s, (x + r + 7) * s, (y + r + 7) * s], fill=(17, 17, 17, 255))
        g.ellipse([(x - r) * s, (y - r) * s, (x + r) * s, (y + r) * s], fill=colour + (255,))
        g.ellipse([(x - r * 0.45) * s, (y - r * 0.65) * s, (x - r * 0.05) * s, (y - r * 0.25) * s], fill=(255, 255, 255, 90))
    return img


def main() -> None:
    img = draw()
    png = ROOT / "src/minimol/static/icon.png"
    img.resize((256, 256), Image.LANCZOS).save(png)
    ico = ROOT / "packaging/windows/minimol.ico"
    img.save(ico, sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    print(f"wrote {png} and {ico}")


if __name__ == "__main__":
    main()
