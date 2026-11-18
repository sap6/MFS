# Austin RCS Benchmark Suite, problem set III-A (PEC almonds): reference data

Unmodified copies of the files in `IIIA-Reference Data.zip` from
<https://github.com/UTAustinCEMGroup/AustinCEMBenchmarks> (commit 300010c,
folder `Austin-RCS-Benchmarks/Problem III-Almonds/Problem Set IIIA-PEC Almonds`).
Copyright 2023 The University of Texas at Austin; distributed under the terms
in `LICENSE-Austin-RCS.txt` and `LICENSE-CC-BY-SA-4.0.txt` (both copied from
that repository). These files are not covered by this repository's own licence.

Columns: frequency [Hz], theta [deg], phi [deg], monostatic RCS [dBsm].

- `meas.m1_rcs.III.A.s1.fx{1..4}.{H,V}.txt`: measured, 9.936 in almond, at 3.5,
  5.125, 7, 10.25 GHz, HH and VV, every 0.25 degrees.
- `meas.m2_rcs...`: measured on the almond twice as large at 1.75, 2.575, 3.5,
  5.125 GHz, shifted by -10 log10(4) dB (so comparable with the m1 files).
- `ref_rcs...`: the suite's reference simulation (ARCHIE-AIM, finest mesh,
  about 0.6 mm average edge length), every 0.5 degrees.

Source of the measurements: J. T. Kelley, D. A. Chamulak, C. C. Courtney,
A. E. Yilmaz, "Rye Canyon radar cross-section measurements of benchmark almond
targets", IEEE Antennas Propag. Mag. 62(1), 96-106, 2020,
doi:10.1109/MAP.2019.2955702.

Used by `../fig_almond_benchmark.py`.
