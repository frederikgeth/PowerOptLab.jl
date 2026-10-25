# PowerOptLab.jl logo

The mark is a lab flask whose liquid is layered in the four phase colours, with
a small feeder growing out of its neck. The colours are Julia's four colours,
standing for phases a, b and c (red, green, purple) and the neutral (blue).
The style (cream ground, dark-brown outlines, flat fills, Poppins) matches the
[BMOPFTools.jl](https://github.com/frederikgeth/BMOPFTools.jl) logo, so the two
packages read as a family.

All files live in [`docs/src/assets/logo/`](../../docs/src/assets/logo/), each
as SVG, PNG and PDF. Text in the SVGs is converted to outlines, so they render
the same without the font installed.

| File stem | Use |
|---|---|
| `poweroptlab-hero` | README and docs home banner (1600×720; also `@2x` PNG) |
| `poweroptlab-social` | GitHub social preview (1280×640) |
| `poweroptlab-lockup`, `-dark` | Mark + wordmark + tagline, transparent background (also `@2x` PNG) |
| `poweroptlab-lockup-stacked`, `-dark` | Mark above wordmark, for square spaces (also `@2x` PNG) |
| `poweroptlab-mark`, `-dark` | Mark only, full colour, transparent background |
| `poweroptlab-mark-mono`, `-mono-reversed` | One-colour mark for print and single-colour use |
| `poweroptlab-icon`, `-dark` | App-icon tile: avatars, slides, package registries |
| `poweroptlab-mark-small`, `-small-dark`, `poweroptlab-icon-small` | Simplified geometry for 16–48 px; the favicon source |

Wired up elsewhere:

- `docs/src/assets/logo.svg` and `logo-dark.svg`: Documenter's sidebar logo
  (copies of `poweroptlab-mark` and `poweroptlab-mark-dark`).
- `docs/src/assets/favicon.ico`: 16, 32 and 48 px from `poweroptlab-icon-small`,
  registered in `docs/make.jl`.
- GitHub social preview: upload `poweroptlab-social.png` under the repository's
  Settings › General › Social preview. GitHub has no API for this.

## Palette

| Role | Light | Dark |
|---|---|---|
| Ink (outlines, wordmark) | `#26221F` | `#F2EBDD` |
| Paper | `#FBF7EF` | `#26221F` |
| Ground | `#F3EDE2` | `#1A1714` |
| Flask rim, accents | `#D3BC8C` | `#7A6A4C` |
| Phase a / b / c | `#CB3C33` / `#389826` / `#9558B2` | `#D66661` / `#6BAB5B` / `#AA79C1` |
| Neutral | `#4063D8` | `#6682DF` |

## Regenerating

Every file is generated from the geometry in
[`generate_logo_assets.py`](generate_logo_assets.py); edit that, not the
outputs. It needs Python with `fonttools` and `pillow`, the `rsvg-convert`
binary (librsvg), and the Poppins TTF files from
[Google Fonts](https://github.com/google/fonts) (SIL Open Font License):

```bash
python3 -m venv /tmp/logo-venv && /tmp/logo-venv/bin/pip install fonttools pillow
/tmp/logo-venv/bin/python scripts/logo/generate_logo_assets.py --fonts /path/to/fonts
```

The fonts directory must contain `Poppins-Regular.ttf`, `Poppins-Medium.ttf`,
`Poppins-SemiBold.ttf` and `Poppins-Bold.ttf`.
