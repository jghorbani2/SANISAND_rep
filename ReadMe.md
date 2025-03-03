# SANISAND Model Implementation

## Overview
This is a C++ implementation of the SANISAND (Simple ANIsotropic SAND) constitutive model for simulating sand behavior under various loading conditions.

## Configuration

### Environment Variables (.env file)
The program uses a `.env` file for configuration. Create a `.env` file in the project root with these settings:

```
# File paths for input and output
INPUT_FILE_PATH=input_data.txt    # Path to input parameters file
OUTPUT_FILE_PATH=q_vs_epsilon2.csv # Path for simulation results
```

### Input File Format
The input file should contain model parameters in the following format:

```
# Material Properties (all required)
Property: G0 125.0      # Shear modulus constant
Property: K0 150.0      # Bulk modulus constant
Property: Mc 1.25       # Critical state slope in compression
Property: Lambda 0.025  # Compression index
Property: N_c 0.2       # Reference void ratio
Property: alpha_c 1.5   # Dilatancy parameter
Property: n_b 2.0       # Bounding surface parameter
Property: ch 0.968      # Hardening parameter
Property: n_d 1.1       # Dilatancy parameter
Property: h0 7.05       # Hardening modulus
Property: A0 0.704      # Plastic potential parameter
Property: Me 1.25       # Critical state slope in extension
Property: cz 600.0      # Critical state parameter
Property: zmax 4.0      # Maximum fabric tensor
Property: m_iso 0.05    # Isotropic hardening parameter
Property: Patm 101.3    # Atmospheric pressure
Property: P_min 0.5     # Minimum mean stress
Property: STOL 1e-5     # Convergence tolerance
Property: FTOL 1e-5     # Yield tolerance
Property: LTOL 1e-5     # Load/unload detection tolerance

# State Variables (all required)
StateVariable: VoidRatio 0.7
StateVariable: StressXX -100.0
StateVariable: StressYY -100.0
StateVariable: StressZZ -100.0
StateVariable: StressZY 0.0
StateVariable: StressZX 0.0
StateVariable: StressXY 0.0

# Analysis Configuration
AnalysisType: 2        # Analysis type (2, 3, or 4)
NumSteps: 100          # Number of simulation steps
EpsilonYY: -1e-4      # Strain increment
QUP: 50.0             # Upper q limit (required for type 4)
QD: -50.0             # Lower q limit (required for type 4)
```

### Analysis Types
1. Type 2: Modified strain calculation with stress-dependent moduli
   - Uses elastic moduli (G, K) dependent on mean stress
   - Computes strain increments using consistency parameter

2. Type 3: Alternative strain calculation
   - Simpler strain increment calculation
   - Uses plastic multiplier and hardening parameter

3. Type 4: Analysis with q-value limits
   - Controls deviatoric stress path
   - Requires valid QUP and QD values (QUP must be > QD)

### Output Format
Results are written to a CSV file with columns:
- Step: Simulation step number
- Epsilon_2: Accumulated strain value
- p: Mean effective stress
- q: Deviatoric stress
- specificVoid: Current void ratio

## Building and Running

### Prerequisites
- C++ compiler with C++17 support
- Eigen library (version 3.4.0 or later)
- macOS development tools

### Directory Structure
```
project_root/
├── .env                    # Configuration file
├── .vscode/               # VSCode configuration (if you are using vs code, note that these are tailored for MAC OS so chnages are needed for other systems/configurations)
│   ├── c_cpp_properties.json
│   ├── launch.json
│   └── tasks.json
├── input_data.txt         # Input parameters
├── main.cpp               # Main program
├── SANISANDModelUMAT.hpp # Model header
├── SANISANDModelUMAT.cpp # Model implementation
└── UMATBase.hpp          # Base class header
```

### Build Instructions
1. Ensure Eigen is installed:
   ```bash
   brew install eigen
   ```

2. Configure VSCode:
   - Update include paths in `c_cpp_properties.json`
   - Set correct SDK paths
   - Configure debugger in `launch.json`

3. Build the project:
   ```bash
   make
   ```

### Running
1. Create and configure `.env` file
2. Prepare input file with model parameters
3. Run the program:
   ```bash
   ./main
   ```

## Error Handling
The program includes comprehensive error checking:
- Validates all required input parameters
- Checks file permissions and existence
- Validates analysis parameters:
  - NumSteps must be > 0
  - AnalysisType must be 2, 3, or 4
  - EpsilonYY must be nonzero
  - For Type 4: QUP must be > QD, both nonzero

## Debugging
- Use VSCode's built-in debugger
- Check working directory settings in `launch.json`
- Enable debug output for detailed information
- Monitor file operations and parameter validation

## Known Issues
- File permission issues may occur on macOS
- Working directory must be set correctly
- All material properties must be provided

## Analysis Types
1. Type 1: Standard analysis
2. Type 2: Modified strain calculation with stress-dependent moduli
3. Type 3: Alternative strain calculation
4. Type 4: Analysis with q-value limits (QUP and QD)

## Troubleshooting
If you encounter file access issues:
1. Verify file existence and permissions
2. Check working directory
3. Ensure paths in .env are correct
4. Run the program from the project directory

## Overview

This repository contains a C++ implementation of the SANISAND (Simple ANIsotropic SAND) constitutive model for soil mechanics. SANISAND is an advanced critical state-compatible elastoplastic model for sand that can capture important behaviors such as:

- Stress-dilatancy
- Critical state
- Fabric evolution
- Stress-induced anisotropy

This implementation is designed for integration with finite element analysis software and can be used for geotechnical simulations.

## Features

- Full implementation of the SANISAND model with kinematic hardening
- Explicit integration scheme with adaptive substepping for numerical stability
- Support for fabric evolution tensor
- Support for different types of stress paths 

## Requirements

- C++ compiler with C++11 support
- Eigen library (version 3.3 or newer)
- CMake (version 3.10 or newer) for building

## Building the Project

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

The main interface for the SANISAND model is through the `SANISANDModel` class, which implements the `UMATBase` interface. The model can be used in two main ways:

1. As a standalone library for single element tests
2. As a UMAT (User MATerial) subroutine for integration with FEM software

### Example: Running a Simple Test

```cpp
#include "SANISANDModelUMAT.hpp"

int main() {
    // Define material properties
    std::unordered_map<std::string, double> properties;
    properties["G0"] = 125.0;     // Shear modulus constant
    properties["K0"] = 150.0;     // Bulk modulus constant
    properties["Mc"] = 1.25;      // Critical state slope in compression
    properties["Lambda"] = 0.37;  // Compression index
    // ... other properties ...
    
    // Create SANISAND model instance
    SANISANDModel model(properties);
    
    // Setup input and output data structures for stress calculations
    UMATBase::InputData inputData;
    UMATBase::OutputData outputData;
    
    // Define strain increments and initial state
    // ... [initialization code] ...
    
    // Calculate stress increment
    model.calculateStressIncrement(inputData, outputData);
    
    // Process results
    // ... [results processing code] ...
    
    return 0;
}
```

## Input Data Format

The model requires input data in a specific format, which includes:
- Material properties
- State variables (stress, void ratio)
- Custom state variables (fabric tensor, kinematic hardening variables)

See the `input_data.txt` file for an example of the required input format.

## Model Parameters

| Parameter | Description |
|-----------|-------------|
| G0 | Shear modulus constant |
| K0 | Bulk modulus constant |
| Mc | Critical state slope in compression |
| Lambda | Compression index |
| N_c | Reference void ratio |
| alpha_c | Dilatancy parameter |
| n_b | Bounding surface parameter |
| ch | Hardening parameter |
| n_d | Dilatancy parameter |
| h0 | Hardening modulus |
| A0 | Plastic potential parameter |
| Me | Critical state slope in extension |
| cz | Critical state parameter |
| zmax | Maximum fabric tensor value |
| m_iso | Isotropic hardening parameter  (may be set to zero)|
| Patm | Atmospheric reference pressure |
| P_min | Minimum pressure (numerical) |
| STOL | Convergence tolerance |

## Analysis Types

The implementation supports different analysis types:
1. Standard strain-controlled test
2. Triaxial compression test with constant lateral stress
3. Constant volume test
4. Stress-controlled cyclic test

## Citation

If you use this implementation in your research, please cite:

@software{SANISANDcpp2025,
  author       = {Javad Ghorbani},
  title        = {{SANISAND Model C++ Implementation: A Critical 
                  State Constitutive Model for Sand}},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/jghorbani2/SANISAND_rep},
  version      = {1.0.0}
}

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

This implementation is based on the theoretical work on SANISAND model by:
- Dafalias, Y.F. and Manzari, M.T., 2004. Simple plasticity sand model accounting for fabric change effects. Journal of Engineering mechanics, 130(6), pp.622-634.

The CSL used in the model based on the following paper:
- Ghorbani, J. and Airey, D.W., 2021. Modelling stress-induced anisotropy in multi-phase granular soils. Computational Mechanics, 67(2), pp.497-521.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
