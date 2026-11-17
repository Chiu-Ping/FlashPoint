# FlashPointHexShape

Geant4 optical simulation of a regular hexagonal BC-408 plastic scintillator
with configurable side-coupled 6 mm SiPM readout. The requested baseline has
a 400 mm nominal edge length, 10 mm thickness, one SiPM at each of the six
vertices, and one SiPM at the centre of each of the six edges (12 total).
Geometry and SiPM dimensions are selected in the Geant4 macro without
recompiling.

The primary observables are per-SiPM detected-photoelectron arrival time
(`dN/dt`), detected PE count, first/mean time, and integrated avalanche
charge. These provide the amplitude and timing features needed to develop an
incident-position reconstruction.

## Default detector and SiPM

- BC-408 regular hexagon: 400 mm edge length and 10 mm thickness.
- Flat-topped orientation: opposite vertices lie on the X axis.
- Six fixed vertex sensors on short contact faces.
- One centred sensor on each straight edge by default (12 sensors total).
- onsemi MicroFJ-60035-TSV active area: 6.07 x 6.07 mm2.
- Package: 6.13 x 6.13 x 0.37 mm3.
- 22,292 microcells, 35 um pitch, 75% fill factor.
- Nominal +6 V values used by the model: 50% peak PDE at 420 nm and gain
  6.3e6.
- The wavelength-dependent PDE is an approximate digitization of the +6 V
  curve in the supplied MICROJ-SERIES/D datasheet.

The datasheet also gives a 180-250 ps anode-cathode rise time and a 3.0 ns
fast-output FWHM for the 6 mm sensor. These describe the electrical output;
they are documented here but are not incorrectly applied as photon-arrival
time smearing. Add a front-end/shaper model when voltage-pulse amplitude or
electronics-level timing is required.

## Geometry controls

Set all geometry commands before `/run/initialize`:

```text
/detector/edgeLength 400 mm
/detector/thickness 10 mm

/detector/sipmLength 6.07 mm
/detector/sipmHeight 6.07 mm
/detector/sipmPackageLength 6.13 mm
/detector/sipmPackageHeight 6.13 mm
/detector/sipmThickness 0.37 mm
/detector/sipmGain 6.3e6

/detector/edgeSiPMsPerEdge 1
```

`edgeSiPMsPerEdge` may be 0, 1, 2, 3, 4, and so on. The total is

```text
number of SiPMs = 6 + 6 * edgeSiPMsPerEdge
```

For example, 0 gives 6 sensors, 1 gives 12, 2 gives 18, and 4 gives 30.
For `N` sensors on an edge, their centres are placed at the interior fractions
`1/(N+1), 2/(N+1), ..., N/(N+1)` of the straight section between the two
vertex contact faces. Thus the requested `N=1` sensor is at the edge centre.
Construction stops with an error if a requested count would make packages
overlap each other or a vertex contact face.

The nominal `edgeLength` is the side length of the regular hexagon before its
six small vertex contact faces are made. With a 400 mm edge, the bounding XY
span is 800 mm by 692.820 mm. Each vertex is trimmed by
`sipmLength/sqrt(3)` along both adjacent edges so its contact face exactly
matches the configured active SiPM length.

### SiPM IDs

The six vertex sensors are numbered counter-clockwise from +X:

```text
SiPM 0: +X vertex (0 degrees)
SiPM 1: upper-right vertex (60 degrees)
SiPM 2: upper-left vertex (120 degrees)
SiPM 3: -X vertex (180 degrees)
SiPM 4: lower-left vertex (240 degrees)
SiPM 5: lower-right vertex (300 degrees)
```

For `N = edgeSiPMsPerEdge`, the remaining IDs are grouped by edge in the
same counter-clockwise order. Edge 0 joins vertices 0 and 1, edge 1 joins
vertices 1 and 2, and so on through edge 5, which joins vertices 5 and 0:

```text
edge k IDs: 6+k*N ... 6+(k+1)*N-1
```

With the requested `N=1` baseline, SiPMs 6-11 are the centres of edges 0-5.

## Preview before running

From the `run` directory, start the application without a macro:

```bash
./plasticMuon
```

The application loads `macros/vis.mac`. Edit the geometry block at the top of
that file, then restart the viewer. It intentionally stops after drawing the
geometry and does not run events. After inspection, type commands such as:

```text
/source/position 0 0
/run/beamOn 1
```

For a fast non-graphical placement/overlap check (initialization only, with no
events):

```bash
./plasticMuon macros/check_geometry.mac
```

That macro initializes and exits without running a particle.

## Build and run

Load the local Geant4 environment, then:

```bash
cd /Volumes/ExtremeSSD/Software/FlashPointHexShape
cmake -S . -B build
cmake --build build -j8
cd run
./plasticMuon macros/run.mac
```

The main output is `sipm_timing.root`. To change the incident XY point, edit:

```text
/source/position 0 0
```

For a single file with uniformly randomized truth positions, use:

```bash
./plasticMuon macros/position_scan.mac
```

This enables `/source/randomXY true` after initialization and writes
`sipm_position_scan.root`. It is the recommended input for the position-
reconstruction notebook.

## Timing and amplitude outputs

For every configured sensor, the ROOT file contains:

- `sipmN_time`: detected optical-photon arrival time, 0-200 ns with 50 ps
  histogram bins.
- `sipmN_path`: detected-photon optical path length.
- `response_sipmN_xy`: detected-PE response versus true incident XY.

The `hits` tree stores one row per detected photon:

```text
event, sipm, t_ns, wavelength_nm, hit_x_mm, hit_y_mm,
source_x_mm, source_y_mm, creator, path_mm
```

The `events` tree stores one row per muon. For every SiPM N it includes:

```text
n_reach_sipmN
n_detected_sipmN
n_scint_detected_sipmN
n_cherenkov_detected_sipmN
first_time_sipmN_ns
mean_time_sipmN_ns
mean_path_sipmN_mm
charge_sipmN_pc
```

`n_detected_sipmN` is the photoelectron-count amplitude observable.
`charge_sipmN_pc = N_PE * gain * elementary_charge` is the corresponding
integrated sensor charge before front-end shaping. A voltage peak requires a
readout impedance and pulse-shaping model.

Plot normalized `dN/dt` for every available sensor with:

```bash
root -l -b -q 'analysis/plot_dndt.C("sipm_timing.root")'
```

This automatically discovers the runtime SiPM count and produces
`dNdt_all_SiPMs.pdf` and `dNdt_overlay.pdf`. It normalizes each bin as

```text
dN_PE/dt = bin content / (generated events * bin width)
```

The 50 ps bin size is sampling, not a claim of 50 ps detector resolution.

## Position-reconstruction observables

Use a calibration scan over known `(x,y)` positions and reserve a separate
test sample. Useful features are:

- all per-SiPM `N_PE` or `charge_pc` values;
- first-photon times and robust early-photon times (for example the 5th PE);
- time differences between opposing edges and between opposite vertices;
- normalized charge fractions `Q_i / sum(Q)`;
- log charge ratios for opposing sensors;
- total charge and hit multiplicity.

Report reconstruction bias and the central-68% or RMS68 error separately for
X, Y, and radial distance. Do not quote training-sample resolution.

Generate a uniform random-position calibration sample and open the supplied
step-by-step notebook with:

```bash
cd run
./plasticMuon macros/position_scan.mac
cd ../analysis
jupyter lab FlashPoint_spatiotemporal_position_reconstruction.ipynb
```

The notebook automatically discovers the SiPM count and compares calibrated
charge centre-of-gravity, amplitude asymmetry, amplitude ridge regression,
first-time ridge regression, combined amplitude-plus-time regression, and
weighted nearest neighbours.  Change `ROOT_FILES` in its configuration cell
to combine several statistically independent ROOT files.

## Fair comparison with an RPC

Use the same active area, muon angular/energy distribution, trigger region,
position grid, event count, dead regions, threshold definition, and reference
clock. Also state whether the comparison holds readout-channel count,
segmentation, bandwidth, and material budget fixed; equal area alone is not a
fully controlled comparison.

Recommended primary metrics:

| Metric | Plastic/SiPM definition | RPC comparison |
|---|---|---|
| Detection efficiency | Fraction with `sum(N_PE)` above a threshold fixed from the noise study | Same trigger and false-hit probability |
| Time resolution | Width of reconstructed time minus true crossing time after position correction | Same clock and estimator convention |
| X/Y position resolution | Test-sample residual width plus bias | Same illuminated grid and fiducial cuts |
| Response uniformity | P5/P50/P95 of efficiency, charge, and timing over XY | Same map binning |
| Edge/corner degradation | Metric change in edge bands and corner regions | Same excluded/included boundary widths |
| Multiplicity | Number of sensors/strips over threshold per muon | Compare clustering and ambiguity |
| Noise occupancy | False hits per area and time window | Same threshold and gate |
| Rate capability | Efficiency and timing versus incident flux | Same flux and run duration |
| Dead time | Lost-event fraction versus rate | Include recovery/readout effects |
| Material budget | Radiation/interaction length and mass per m2 | Include supports and electronics |
| System cost/power | Channels, power, cooling, and cost per m2 | Compare complete detector systems |

A practical optimization should scan `edgeSiPMsPerEdge = 0, 1, 2, 3, 4` and
compare Pareto curves: efficiency versus channel count, timing resolution
versus channel count, and position resolution versus channel count. A useful
project-level gate is to require the scintillator design to remain within the
agreed tolerance of the measured RPC baseline for efficiency, time, and
position resolution simultaneously, rather than matching only the best one.

## Model limitations to address before an RPC claim

The current optical model does not yet include SiPM dark noise, crosstalk,
afterpulsing, microcell saturation/recovery, temperature dependence, coupling
layer losses, mechanical gaps, surface roughness variations, electronics
noise, discriminator threshold, or time walk. The datasheet gives 50-150
kHz/mm2 dark count, 8-25% crosstalk, 0.75-5% afterpulsing, and a 50 ns
microcell recharge constant for the 60035, depending on overvoltage. Add these
effects before making a quantitative claim that the panel matches an RPC.
