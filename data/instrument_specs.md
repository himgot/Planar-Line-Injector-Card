# Instrument specifications, compiled 17 Nov 2026

Every figure below was read from the manufacturer's page or datasheet, opened on 17 Nov 2026. Figures are the manufacturer's published values, **not** calibration data for the author's units. Nothing is known about the calibration state of the author's instruments → GAP for all three.

The datasheets were read as extracted text, not as rendered pages. Where the extraction lost symbols, this is noted.

---

## 1. MECO LCR999A LCR meter

**Primary source opened:** MECO catalogue sheet "LCR999A" (PDF, title "105 (LCR 999A).cdr")
- https://www.meco.in/uploads/products/catalog/LCR999A.pdf
- Same sheet, distributor copy: https://pravininstruments.com/download/Meco%20LCR999A%20LCR%20Meter.pdf (title "101.cdr"); identical text.
- Retail pages opened for cross-checks (secondary): https://electrixzone.com/product/meco-lcr-meter-lcr-999a, https://www.toolbuy.com/tp/meco-lcr-meter-lcr-999a, https://envmart.com/63317-meco-lcr-meter-lcr999a.html

| Item | Published figure (verbatim where possible) | Source |
|---|---|---|
| Test frequencies | "100Hz, 120Hz, 1KHz, 10KHz" | catalogue |
| Test signal level | "0.6 Vrms" | catalogue |
| Output (source) impedance | "Constant output impedance : 100 ohm" | catalogue |
| Basic accuracy | "Basic Accuracy 0.1%" / "Basic measurement accuracy 0.1% & speed upto 10 meas / sec" | catalogue |
| Equivalent circuit | "Series, Parallel" ("Equivalent Mode") | catalogue |
| Auto LCR | "Manual / Auto"; "Automatic identification function (Ai)" | catalogue |
| Ranging | "Auto" | catalogue |
| Correction | "Short , Open" | catalogue |
| Test terminals | "2-Terminal, 4-Terminal, 5-Terminal"; accessory "5 Terminal Kelvin Test Leads x 1", "Short Circuit Plate x 1" | catalogue |
| Speed | "10 meas/sec (fast), 5 meas/sec (med), 2 meas/sec (slow)" | catalogue |
| L measuring range | text extracts as "L 4 H ~ 1000H Range for Display 0.001 H ~ 1000.0H" followed by two stray glyphs; the prefix glyph before the first "H" and before "0.001 H" is lost in extraction. Read as **4 µH – 1000 H measuring range, 0.001 µH – 1000.0 H display range**: the lower display limit "0.001 µH" is confirmed with the µ glyph intact on envmart.com ("inductance (0.001 µH to 1000 H)"); the "4 µH" lower measuring limit is an inference from context (the same reading is used for the "4 µH nominal lower range"). | catalogue + envmart (secondary) |
| R / Z measuring range | "R / Z 0.4 ~ 10M Range for Display 0.0001 ~ 10.000M" (Ω glyphs lost) → 0.4 Ω – 10 MΩ measuring, 0.0001 Ω – 10.000 MΩ display | catalogue |
| D, Q display | "D Range for Display 0.0000 ~ 9.999, Resolution : 0.0001"; "Q Range for Display 0.0000 ~ 9999, Resolution : 0.0001" | catalogue |
| Phase display | "Range for Display - 179.9 ~ 179.9 , Resolution : 0.01" (θ glyph lost; unit presumably degrees) | catalogue |
| Interface | "Mini-USB (Virtual Serial Port)", "SCPI compatible" | catalogue |

**GAPS (LCR999A):**
- GAP: no accuracy table versus range, frequency, D/Q or measurement speed was found (the catalogue gives only "Basic Accuracy 0.1%"). The user manual was not found online. Accuracy at 1 kHz and 10 kHz for ~18 µH / D ≈ 6.5–65 is therefore **not specified** by any source opened.
- GAP: whether "4 µH" is the lower L limit at all frequencies or only at 10 kHz is not stated.
- GAP: phase/D accuracy, the auto Ls/Lp selection rule, and the Ai rule for low-Q parts are not stated.
- Not used: a search snippet for the different model "LCR 999" (non-A; 0.3 % basic accuracy, 0.3 Vrms, L 0.1 µH–999.9 H at 1/10 kHz; https://5.imimg.com/data5/SELLER/Doc/2025/3/496822853/EF/BD/RN/939413/meco-lcr999a-digital-lcr-meter.pdf). That page was opened, but it describes "LCR 999", not "LCR999A", so it is **not** applied to the author's meter.

**Generic convention used only as an assumption (different manufacturer):** IET Labs DE-5000 handheld LCR meter manual, https://www.ietlabs.com/pdf/Manuals/DE_5000_im.pdf (opened): "Adjustment to accuracy (Z) based on dissipation (D) reading: D > 0.1: Z * √(1+D2)"; "In inductance mode, D ≤ 0.1: ZL = 2πƒL"; "Phase angle Θ accuracy: ΘZ= ±(180/π)*AZ (deg)". The uncertainty budget applies this convention with the LCR999A's 0.1 % basic accuracy as a **best-case lower bound**. It is not an LCR999A specification.

---

## 2. Keysight DSO1052B oscilloscope

**Primary source opened:** Keysight "DSO1000A/B Series Portable Oscilloscopes Data Sheet", 5989-9368EN (published 26 Apr 2018)
- https://www.keysight.com/content/dam/keysight/en/doc/ungate/data-sheets/5989-9368.pdf

| Item | Published figure | Note |
|---|---|---|
| Bandwidth (−3 dB) | "DSO1052B: DC to 50 MHz" | warranted (footnote 1); footnote 2: "20 MHz (when vertical scale is set to < 5 mV)" |
| Sample rate (B models) | "1 GSa/sec half channel interleaved, 500 MSa/sec all channels" | |
| Memory (B models) | "16 kpts half channel interleaved, 8 kpts all channels" | |
| Vertical resolution | "8 bits" | |
| Vertical sensitivity | "2 mV/div to 10 V/div" | |
| DC gain accuracy | "2 mV/div to 5 mV/div: ± 4.0% full scale (A and B models)"; "10 mV/div to 10 V/div: ± 3.0% full scale (B models only)" | warranted; valid after 30-min warm-up and ±10 °C from firmware-calibration temperature |
| Input impedance (B models) | "1 MΩ ± 2% in parallel with 15 pF ± 3 pF" | |
| Selectable BW limit | "20 MHz" | |
| Dynamic range | "+6 divisions from center screen" (sic) | |
| Max input | "CAT I 300 Vrms, 400 Vpk; transient overvoltage 1.6 kVpk" | |
| Probes, standard | "10:1 passive probe for each input channel"; "N2862B 150 MHz 10:1 passive probe (standard with 50, 60, 70, 100 MHz models)" | |
| Probe options | "10070D 20 MHz 1:1 passive probe"; "N2863B 300 MHz 10:1"; "10076B 250 MHz, 100:1, 4 kV"; "N2791A 25 MHz, 700V differential probe" | |
| Save formats | waveforms "CSV saved"; images "8-bit BMP, 24-bit BMP, PNG saved" | relevant: CSV capture was possible but none is archived |

**Graticule and probe-factor handling.** Source: Keysight "DSO1000 Oscilloscope Educator's Training Resources", manual part 54136-97000 (May 2012), https://www.keysight.com/dz/en/assets/9018-03712/training-materials/9018-03712.pdf (opened):
- "Our X-axis (horizontal) consists of 12 major divisions (if the menu display is off)", with "5 minor divisions (0.2 divisions each)"; "Our Y-axis (vertical) consists of 8 major divisions vertically".
- "when using an Keysight 1000 Series oscilloscope you must either manually enter the probe attenuation factor (10:1) or press the [Default Setup] front panel key to set the probe attenuation to 10:1". This supports the possibility of a ×1/×10 menu-vs-probe mismatch. The actual state remains UNKNOWN.

**Probe catalogue (secondary, Keysight probe resource centre):** https://prc.keysight.com/Content/Browse%20by%20Scope/1000_series.html lists N2862B as "Passive | BNC | 150 MHz | 10:1 | 10 MΩ | 300V" and the 1:1 probes "10700D … 20 MHz | 1:1 | 1 MΩ, 70 pF" (the model number there differs from the data sheet's "10070D").

**GAPS (DSO1052B):** the probe actually used, its ×1/×10 switch position and the channel probe-factor menu setting are UNKNOWN. N2862B input capacitance is not verified (the product page was bot-blocked; only a search snippet was seen, which is not used). No AC (frequency-dependent) vertical accuracy is specified: the 3 % is a DC full-scale figure.

---

## 3. Scientific SM5072 function generator (5 MHz)

**Primary source opened:** Scientific Mes-Technik "5 MHz General Purpose Function Generator SM5072 – Technical Specifications" (manufacturer-branded data sheet hosted by a distributor)
- https://5.imimg.com/data5/SELLER/Doc/2025/3/496629217/JP/KX/NI/1605384/scientific-sm5072-5-mhz-general-purpose-function-generator.pdf
- The manufacturer product page https://www.scientificindia.com/product-details/5-mhz-general-purpose-function-generator returned a bot check (not readable).

| Item | Published figure |
|---|---|
| Frequency range | "0.5 Hz to 5 MHz" |
| Waveforms | "Sine, Square, Triangle, DC" (+ TTL trigger output) |
| Output voltage | "10 Vpp into 50 Ω, max. 20 Vpp open circuit" ("All waveform characteristics are with 50 Ω termination") |
| Attenuation | "max. 60 dB, 2 steps : 20 dB ± 0.2 dB each"; "Variable : 0 to 20 dB" |
| Amplitude flatness (sine/triangle) | "0.5 Hz - 0.5 MHz : max. 0.2 dB"; "0.5 MHz - 5 MHz : max. 0.5 dB" |
| Output impedance | "switchable 50 Ω / 600 Ω" |
| Sine distortion | "0.5 Hz - 100 kHz : max. 0.5 %"; "0.1 MHz-0.5 MHz : max. 1.5%"; "0.5 MHz - 5 MHz : max. 3%" |
| Frequency display accuracy | "5 Hz - 5 MHz : ± (5 × 10^-5 + 1 D)" |
| DC offset | "max. ± 2.5 V into 50 Ω, max. ± 5 V open circuit" |
| Square rise time | "typ. < 30 ns"; overshoot "< 5%" |

(Ω glyphs were blank in the extracted text and are restored here from context, e.g. "50 [Ω] termination".)

**GAPS (SM5072):** the states of the two −20 dB buttons and of the 50 Ω/600 Ω switch during the photo sessions are UNKNOWN. The flatness figure applies into a 50 Ω termination, not into the low-impedance AC winding (7.31 Ω dc in series with jωL, L ≈ 18 µH), so the delivered amplitude versus frequency into the card is not given by this spec.
