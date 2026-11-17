# FlashPoint: Position Reconstruction in Monolithic Plastic Scintillators

FlashPoint is a detector research and development project focused on
high-resolution particle position reconstruction in large-area
monolithic plastic scintillators with distributed silicon
photomultiplier (SiPM) readout.

This repository combines Geant4 optical simulations and machine-learning
reconstruction to optimize SiPM configurations and investigate detector
performance under different readout constraints.

## Key Features

- Geant4-based optical photon generation, transport, and detection.
- Configurable square and hexagonal scintillator geometries.
- Flexible corner, edge, and face-mounted SiPM configurations.
- Muon position reconstruction using light patterns and early photon timing.
- SiPM layout and readout-channel optimization.
- Evaluation of spatial resolution, uniformity, and timing sensitivity.

## Repository Structure

| Directory | Description |
|---|---|
| [square/](square/) | Square detector simulation and reconstruction |
| [hexagonal/](hexagonal/) | Hexagonal detector simulation and reconstruction |

## Software Requirements

![Geant4](https://img.shields.io/badge/Geant4-11.4.0-blue)
![CMake](https://img.shields.io/badge/CMake-3.31.0-green)
![Python](https://img.shields.io/badge/Python-3.13.0-blue)

Development environment:

- Geant4 11.4.0
- CMake 3.31.0
- Apple Clang 21.0.0
- Python 3.13.0
- Python packages: NumPy, pandas, SciPy, scikit-learn, Matplotlib, and Jupyter.

Geometry-specific compilation and execution instructions are provided
in the corresponding directories.

## Associated Publication

**FlashPoint: Joint Optimization of SiPM Layout and Early Photon Timing for Muon Localization in Large-Area Monolithic Plastic Scintillators**

Manuscript in preparation.

## Contact

Qiuping Shen (shenqp@ihep.ac.cn)  
Xilei Sun (sunxl@ihep.ac.cn)

Institute of High Energy Physics, Chinese Academy of Sciences.

