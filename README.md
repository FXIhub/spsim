[![Build Status](https://travis-ci.org/FXIhub/spsim.svg?branch=master)](https://travis-ci.org/FXIhub/spsim)

# SPSIM: X-ray Diffraction Pattern Simulation

SPSIM is a program designed to realistically simulate X-ray diffraction patterns from Protein Data Bank (PDB) files. It models various aspects of the diffraction experiment, including detector characteristics, beam properties, and noise effects, to generate simulated diffraction data.

## Features

*   **PDB File Input**: Reads atomic coordinates and element information from PDB files.
*   **Chemical Formula Input**: Can generate random molecular structures based on a chemical formula.
*   **Diffraction Simulation**: Calculates diffraction patterns based on theoretical models.
*   **Detector Modeling**: Simulates various detector properties such as distance, pixel size, quantum efficiency, readout noise, and dark current.
*   **Noise Simulation**: Incorporates Poisson noise and readout noise for realistic output.
*   **Real-Space and Reciprocal-Space Analysis**: Provides outputs for both real-space electron density and reciprocal-space diffraction patterns.
*   **Multiple Output Formats**: Generates output in VTK and CXI formats for visualization and analysis.
*   **Parallel Processing**: Supports MPI for distributed computation.
*   **Vectorization and CUDA Support**: Optimizes calculations for performance.

## Installation

For detailed installation instructions, please refer to the User Manual: [doc/UserManual.pdf](doc/UserManual.pdf).

A general outline of the build process is as follows:

1.  Extract the source code.
2.  Inside the source code directory create a `build` directory.
3.  Inside the `build` directory run `ccmake ..` followed by `make` and optionally `make install`.

This will produce a binary executable named `spsim`.

You can enable and disable multiple options inside ccmake.

## Usage

### Configuration (`spsim.conf`)

The program reads its configuration from a file named `spsim.conf`. All values are expected to be in SI units. Key parameters include:

*   `number of dimensions`: Currently only `2` is supported.
*   `input type`: Specifies the input structure type (e.g., `"pdb"`).
*   `pdb filename`: Path to the PDB file.
*   Detector parameters: `detector_distance`, `detector_width`, `detector_height`, `detector_pixel_width`, etc.
*   Experiment parameters: `experiment_wavelength`, `experiment_exposure_time`, `experiment_beam_intensity`, etc.

**Example `spsim.conf`:**

```
number_of_dimensions = 2;
input_type = "pdb";
pdb_filename = "DNA_triangle.pdb";
# 5cm detector distance
detector_distance = 5.000000e-02;
# 26.8mm x 26mm CCD
detector_width = 2.680000e-02;
detector_height = 2.600000e-02;
#20um square pixels
detector_pixel_width = 2.000000e-05;
detector_pixel_height = 2.000000e-05;
# 15% quantum efficiency
detector_quantum_efficiency = 1.500000e-01;
# 3.6 eV
detector_electron_hole_production_energy = 5.800000e-19;
detector_readout_noise = 1.000000e+01;
detector_dark_current = 1.000000e-01;
detector_linear_full_well = 2.000000e+05;
# 4x4 binning
detector_binning = 4;
detector_maximum_value = 65535.0;
#13 nm
experiment_wavelength = 1.300000e-08;
# 100fs
experiment_exposure_time = 1.000000e-13;
# 1e15 photons, 100um^2 area
experiment_beam_intensity = 1.0e25;
```

### Running the Simulation

After configuring `spsim.conf`, run the simulation from the command line:

```bash
./spsim
```

### Output Files

The simulation generates several output files, typically in VTK format, which can be visualized using tools like VisIt:

*   `spsim.confout`: Contains the configuration read by the program.
*   `scattering_factor.vtk`: Scattering factors of the input molecule on the detector.
*   `thomson correction.vtk`: Thomson correction factor for each pixel.
*   `solid angle.vtk`: Solid angle for each pixel.
*   `photon count.vtk`: Number of photons detected by each pixel.
*   `electrons per pixel.vtk`: Electrons generated on each pixel.
*   `real output.vtk`: Detector output including noise effects.
*   `noiseless output.vtk`: Detector output without noise.

## Theory

SPSIM simulates diffraction patterns using principles of crystallography and scattering theory. The process involves:

1.  **Reciprocal Space Coordinates**: Calculating reciprocal space coordinates (h, k, l) for detector pixels using the Ewald sphere construction.
2.  **Molecular Scattering Factors**: Computing atomic scattering factors using CCP4's `atomsf.lib` and then the overall molecular scattering factor.
3.  **Thomson Correction**: Applying corrections based on beam polarization and scattering geometry.
4.  **Photon Count**: Calculating the number of photons detected per pixel, considering quantum efficiency and Poisson noise.
5.  **Electron Generation**: Determining the number of electrons generated based on detected photons.
6.  **Detector Output**: Simulating the final detector output, including readout noise and other effects.

For a comprehensive understanding of the underlying theory, please consult the User Manual: [doc/UserManual.pdf](doc/UserManual.pdf).

## Python Interface

The project includes Python bindings (`spsim_pybackend.py`) that allow for programmatic interaction with the simulation core.

## License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

