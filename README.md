# Planar PCB Line Injector Card

**Version 1.0, 21 November 2026** · Himanshu Gautam, Independent Researcher, India · ORCID [0009-0008-5291-496X](https://orcid.org/0009-0008-5291-496X)
Archived version (DOI):** https://doi.org/10.5281/zenodo.24417902

A credit-card-size printed circuit board with no mounted components that works as a simple line injector. Two interleaved planar spirals in the copper form a 1:1 air-core coupled inductor. A signal generator drives one spiral, the DC supply rail passes through the other, and the generator's signal is added in series with the rail. This is the basic job of the line injector used to measure the power-supply ripple rejection (PSRR) of a voltage regulator.

This repository holds everything behind the paper: the paper and its LaTeX source, the Gerber and Altium design files, the measured data, the photographs, the simulation files and the Python scripts that compute every number in the paper.

**Paper:** [`paper/Planar_Line_Injector_Card_v1.0.pdf`](paper/Planar_Line_Injector_Card_v1.0.pdf) (preprint, not peer reviewed)

## Main results

| Quantity | Calculated from the Gerber files | Measured (one card) |
|---|---|---|
| Self-inductance L (each winding) | 17.6 µH | 17.95 µH |
| Mutual inductance M | 16.6 µH | 16.3 µH |
| Coupling coefficient k | 0.94 | 0.91 |
| Winding resistance (DC) | 5.82 Ω (1 oz copper assumed) | 7.31 Ω |
| Inter-winding capacitance | 126–186 pF (159 pF nominal) | 274 pF (probably includes the test leads) |

- Bench test: on a 5 V rail with a 111.6 Ω resistive load, the generator's signal appeared on the DC side from 100 kHz to 5.18 MHz. Below 1 MHz the frequency response agrees with the coupled-inductor model.
- Predicted usable band with a 50 Ω generator: about 0.33 MHz to 5–10 MHz.

## Limitations (please read)

- Only one card was measured, with a handheld LCR meter, and the oscilloscope screens were photographed with a phone.
- The card has **not** been used to measure a real regulator's PSRR.
- The card is **not rated**. The design notes printed on the card ("MAX 30 V · 500 mA", "k ≈ 0.8–0.9", "DCR 5–7 Ω") are not supported by the measurements. Use it only on low-voltage rails (5 V or less) at low current (100 mA or less) from a current-limited supply. At 100 mA the winding drops about 0.73 V.
- There is no isolation rating; the windings have not been hi-pot tested.

The paper lists every open item in its Appendix A, and missing or unverified information is printed in red.

## Repository contents

| Folder | What it is |
|---|---|
| `paper/` | The paper (PDF) and its LaTeX source in `paper/src/` (`main_fallback.tex` builds the PDF; `main.tex` uses the IEEEtran class). |
| `hardware/` | Gerber files (`Planar_Line_Injector_gerbers.zip`, with SHA-256 checksum) and the Altium Designer project in `hardware/altium/`. Read `hardware/README.md` before fabricating. |
| `card_images/` | Renders and layer views of the card artwork. |
| `photos/` | Photographs of the card and the bench, layout-tool screenshots and the oscilloscope screens of the bench test. |
| `data/` | LCR-meter and multimeter readings, oscilloscope screen readings, the reconciled bench table, geometry tables, instrument specifications and all computed values (`values*.json`). |
| `analysis/` | Python scripts that compute every number and figure in the paper. |
| `simulation/` | The initial LTspice schematic, the netlists of the initial and corrected circuit models (LTspice syntax) and their results. |
| `LICENSE`, `LICENSES/`, `LICENSING.md` | Licence texts and which licence applies to which files. |
| `CITATION.cff` | Citation information. |
| `MANIFEST.sha256` | SHA-256 checksums of all files. |

## How the card is connected

| Terminals | Connect to |
|---|---|
| `AC_IN`, `AC_RTN` | Signal generator (the ripple you want to inject). No DC here. |
| `DC_IN`, `DC_RTN` | DC supply (current-limited). |
| `OUT`, `OUT_RTN` | The device under test. The rail passes from `DC_IN` to `OUT` through the second winding. |

The unpopulated "Noise injection" footprint on the right edge connects to `AC_IN` (signal pad) and `AC_RTN` (other pads) and can take a connector.

## Rebuilding the numbers and the paper

Requirements: Python 3 with numpy, scipy, pandas, matplotlib and Pillow; tesseract (used by the geometry script to read the legend); TeX Live with siunitx, circuitikz and booktabs.

Run from the repository folder, in this order:

```
python3 analysis/gerber_geometry.py
python3 analysis/closed_form.py
python3 analysis/inductance_neumann.py      # about 5 min
python3 analysis/capacitance.py             # about 2 min
python3 analysis/theory_bench_model.py
python3 analysis/fig_windings.py
python3 analysis/k_extraction.py
python3 analysis/reconcile_bench.py
python3 analysis/circuit_models.py
python3 analysis/derived_values.py
python3 analysis/headline_values.py
python3 analysis/table_values.py
python3 analysis/build_values.py
python3 analysis/band_values.py
python3 analysis/fabrication_values.py
python3 analysis/bench_shape_fit.py
python3 analysis/figure_sim_vs_bench.py
python3 analysis/figure_scope_screens.py
python3 analysis/build_values.py
```

`build_values.py` writes `data/values.json` and `paper/src/values.tex`. Then compile `paper/src/main_fallback.tex` with pdfLaTeX and BibTeX (for example `latexmk -pdf main_fallback.tex`). Intermediate files go to `work/`. The distributed model (M3) uses the precomputed section inductance matrices in `data/model_inputs/`.

## Contributing measurements

The author does not plan further measurements on this card. Independent measurements are very welcome, including results that disagree. Useful contributions are listed at the end of the paper; the most useful are:

1. LCR measurements with Kelvin clips and a recorded open/short correction;
2. a swept injection measurement with recorded generator settings, a reference channel, and screens and data saved by the oscilloscope (not photographed);
3. a PSRR measurement of a real regulator, compared against a reference injector;
4. measurements of a second card, made from these Gerber files.

Please open an issue or a pull request in this repository with your data and a short description of the setup. Contributed results can be added to later versions with credit to the contributor.

## Licence

- Paper, documentation, photographs, renders and data: [CC BY 4.0](LICENSES/CC-BY-4.0.txt).
- Hardware design files (Gerber and Altium) and scripts: [CERN-OHL-P v2](LICENSE).

See [`LICENSING.md`](LICENSING.md) for details.

## Citation

Please cite the paper if you use this work. Citation details are in [`CITATION.cff`](CITATION.cff) (GitHub shows them under "Cite this repository").
