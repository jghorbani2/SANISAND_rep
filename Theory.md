# SANISAND Model: Technical Documentation

## 1. Theoretical Framework

### 1.1 Overview of SANISAND

The SANISAND (Simple ANIsotropic SAND) model is an elastoplastic constitutive model for sands based on the critical state soil mechanics framework. The model incorporates several key features:

- Bounding surface plasticity formulation
- State-dependent dilatancy
- Kinematic hardening to account for fabric anisotropy
- Fabric tensor evolution to capture cyclic mobility

The model is capable of simulating the mechanical behavior of sands under various loading conditions, including monotonic and cyclic loading paths.

### 1.2 Elasticity

The elastic behavior is defined using pressure-dependent bulk and shear moduli:

- Shear modulus: $G = G_0 \cdot p_{atm} \cdot (p/p_{atm})^{0.5} \cdot (2.97 - e)^2 / (1 + e)$
- Bulk modulus: $K = K_0 \cdot p_{atm} \cdot (p/p_{atm})^{2/3} \cdot (e+1) / e$

Where:
- $G_0$ and $K_0$ are dimensionless elastic constants
- $p$ is the mean effective stress
- $p_{atm}$ is the atmospheric pressure (reference pressure)
- $e$ is the void ratio

### 1.3 Critical State Line

The critical state line (CSL) is defined in $e-\ln p$ space as:

$e_c = N_c - \lambda \cdot \ln(p + \alpha_c)$

Where:
- $e_c$ is the critical state void ratio
- $N_c$ is the reference void ratio
- $\lambda$ is the slope of the CSL
- $\alpha_c$ is a parameter that controls the shape of the CSL

### 1.4 Yield/Loading Surface

The yield surface is defined in the stress ratio space:

$\sqrt{(\mathbf{s}/p' - \boldsymbol{\alpha}):(\mathbf{s}/p' - \boldsymbol{\alpha})} - m = 0$

Where:
- $\mathbf{s}$ is the deviatoric stress tensor
- $p'$ is the mean effective stress
- $\boldsymbol{\alpha}$ is the backstress (kinematic hardening) tensor
- $m$ is a small constant defining the size of the yield surface (may be assummed zero)

### 1.5 Flow Rule

The plastic strain increment is determined by a non-associative flow rule:

$d\boldsymbol{\varepsilon}^p = \langle d\lambda \rangle \cdot \mathbf{R}$

Where:
- $d\lambda$ is the plastic multiplier
- $\mathbf{R}$ is the direction of plastic flow, given by:
  $\mathbf{R} = \mathbf{n} + \frac{1}{3}D \cdot \mathbf{I}$
- $\mathbf{n}$ is the normalized deviatoric direction
- $D$ is the dilatancy coefficient
- $\mathbf{I}$ is the second-order identity tensor

### 1.6 Dilatancy

The dilatancy coefficient is defined as:

$D = A_0 \cdot (1 + \langle \mathbf{z}:\mathbf{n} \rangle) \cdot (M_d \cdot exp({n_d \cdot \psi}) - \boldsymbol{\alpha}:\mathbf{n})$

Where:
- $A_0$ is a model parameter
- $\mathbf{z}$ is the fabric tensor
- $M_d$ is the critical state stress ratio modified for the current Lode angle
- $n_d$ is a model parameter
- $\psi = e - e_c$ is the state parameter
- $\boldsymbol{\alpha}$ is the backstress tensor

### 1.7 Hardening Law

The backstress (kinematic hardening) evolution is governed by:

$d\boldsymbol{\alpha} = \langle d\lambda \rangle \cdot h \cdot (M_b \cdot \mathbf{n} - \boldsymbol{\alpha})$

Where:
- $h$ is the plastic hardening modulus
- $M_b$ is the bounding stress ratio, defined as $M_b = M \cdot exp(<{-n_b \cdot \psi}>)$
- $n_b$ is a model parameter

### 1.8 Fabric Evolution

The fabric tensor evolves to capture sand behavior under cyclic loading:

$d\mathbf{z} = -c_z \cdot \max(-d\varepsilon^p_v, 0) \cdot \left(\sqrt{\frac{2}{3}} \cdot z_{max} \cdot \mathbf{n} + \mathbf{z}\right)$

Where:
- $c_z$ is a model parameter
- $z_{max}$ is a parameter associated with the fabric tensor
- $D$ is the dilatancy coefficient

## 2. Numerical Implementation

### 2.1 Stress Integration Algorithm

The implementation uses an explicit integration scheme with adaptive substepping to ensure numerical stability. 
  

### 2.2 Stress Invariants and Lode Angle Handling

The implementation includes special handling of stress invariants and Lode angle to ensure proper model behavior at different stress paths:

- Mean stress: $p = \frac{1}{3}(\sigma_1 + \sigma_2 + \sigma_3)$
- Deviatoric stress: $\mathbf{s} = \boldsymbol{\sigma} - p \cdot \mathbf{I}$
- Second invariant: $J_2 = \frac{1}{2} \mathbf{s}:\mathbf{s}$
- Third invariant: $J_3 = \det(\mathbf{s})$

The critical state stress ratio $M$ is adjusted based on the Lode angle to handle different stress paths:

$M(\theta) = M_c \cdot \left(\frac{2c^4}{1+c^4-(1-c^4)\cos(3\theta)}\right)^{1/4}$

Where $c = M_e/M_c$ is the ratio of extension to compression critical state stress ratios.

### 2.3 Special Numerical Considerations

Several numerical techniques are implemented to ensure robust performance:

1. **Substepping with error control**: Automatic adjustment of integration step size based on error estimates
2. **Stress scaling for low pressures**: Implementation of P_min to avoid numerical issues at very low stresses
3. **Fabric tensor initialization**: Special initialization of fabric tensor to ensure proper cyclic behavior
4. **Limit checks**: Prevention of unphysical values through various limit checks

## 3. Code Structure

### 3.1 Class Hierarchy

The implementation follows an object-oriented approach with the following classes:

1. **UMATBase**: Abstract base class defining the interface for material models !! stay tuned for updates !!
2. **SANISANDModel**: Concrete implementation of SANISAND model
3. **StressInvariants**: Helper class for stress invariant calculations

### 3.2 Key Methods

The key methods in the implementation include:

- `computeStressInvariants`: Calculates stress invariants and related quantities
- `GradientOfYieldToSigma`: Computes the gradient of the yield function
- `GradientOfPlasticToSigma`: Computes the gradient of the plastic potential
- `PlasticModulus`: Calculates the plastic modulus for consistency condition
- `ComputeDlambdaAndUpdateStressInc`: Implements the return mapping algorithm
- `calculateStressIncrement`: Main method for stress integration
- `computeModuli`: Calculates the tangent modulus for FEM implementation

### 3.3 Data Structures

The implementation uses several specialized data structures:

- **InputData**: Contains strain increments and state variables
- **OutputData**: Contains stress increments and updated state variables
- **StressInvariants**: Contains stress invariants and related quantities

## 4. Analysis Types

The implementation supports different types of analyses through the `AnalysisType` parameter:

1. **Type 1**: Undrained monotonic triaxial
2. **Type 2**: Drained monotonic triaxial
3. **Type 3**: Constant-P monotonic 
4. **Type 4**: Stress-controlled cyclic test

Each analysis type applies specific constraints to the strain increments to simulate different laboratory tests.

## 5. Validation and Verification

The model has been validated against Toyuora Sand Data.

## 6. References

1. Dafalias, Y.F. and Manzari, M.T. (2004). "Simple Plasticity Sand Model Accounting for Fabric Change Effects." Journal of Engineering Mechanics, ASCE, 130(6), 622-634.

2. Dafalias, Y.F. and Manzari, M.T. (2018). "Simple Plasticity Sand Model Accounting for Fabric Change Effects—Part II: Calibration and Validation." Journal of Engineering Mechanics, ASCE, 144(7), 04018048.

3. Taiebat, M. and Dafalias, Y.F. (2008). "SANISAND: Simple Anisotropic Sand Plasticity Model." International Journal for Numerical and Analytical Methods in Geomechanics, 32(8), 915-948.

4. Ghorbani, J. and Airey, D.W. (2021). "Modelling stress-induced anisotropy in multi-phase granular soils." Computational Mechanics, 67(2), pp.497-521.
