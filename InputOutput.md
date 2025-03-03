# Manual for Soil Behavior Simulations using SANISAND

## 1. Overview

This manual provides a guide for running different types of soil behavior simulations using the SANISAND plasticity framework. The available simulations include:

- **Cyclic Controlled Stress Test (AnalysisType: 4)**
- **Drained Monotonic Test (AnalysisType: 2)**
- **Constant Mean Stress (p-constant) Monotonic Test (AnalysisType: 3)**
- **Undrained Test (AnalysisType: 1)**

Each of these tests uses specific stress paths and material parameters, which are defined in the input files.

**Default Notation:** In this implementation, **stress and strain in compression are defined as negative**.

---

## 2. Input File Structure

The input file (commonly named `input_data.txt`) is divided into several sections. Below is an explanation of each section with sample values:

### 2.1 Properties

These are material and model parameters required for the simulation:

```
Property: G0 125.0
Property: K0 150.0
Property: Mc 1.25
Property: Lambda 0.37
Property: N_c 18.7
Property: alpha_c 3370
Property: n_b 1.25
Property: ch 0.968
Property: n_d 2.3
Property: h0 12.0
Property: A0 0.4
Property: Me 0.89
Property: cz 600
Property: zmax 4
Property: m_iso 0.01
Property: Patm 100.0
Property: P_min 1.0e-1
Property: FTOL 1.0e-5
Property: STOL 1.0e-5
Property: LTOL 1.0e-6
```

**Note:** The parameter `FTOL` plays a **dummy role** in this version. It was originally intended to control the yield surface tolerance, but in this implementation, the yield surface is assumed, meaning the stress is always on the yield surface and there is no elasticity.

### 2.2 State Variables

Initial stress conditions and void ratio:

```
StateVariable: VoidRatio 0.735
StateVariable: StressXX -100.0
StateVariable: StressYY -100.0
StateVariable: StressZZ -100.0
StateVariable: StressZY 0.0
StateVariable: StressZX 0.0
StateVariable: StressXY 0.0
```

### 2.3 Custom Variables

Custom variables allow users to define additional model-specific parameters:

```
CustomVariable: AlphaXX 0.0
CustomVariable: AlphaYY 0.0
CustomVariable: AlphaZZ 0.0
CustomVariable: AlphaZY 0.0
CustomVariable: AlphaZX 0.0
CustomVariable: AlphaXY 0.0
CustomVariable: AlphaInitialXX 0.0
CustomVariable: AlphaInitialYY 0.0
CustomVariable: AlphaInitialZZ 0.0
CustomVariable: AlphaInitialZY 0.0
CustomVariable: AlphaInitialZX 0.0
CustomVariable: AlphaInitialXY 0.0
CustomVariable: ZXX 0.0
CustomVariable: ZYY 0.0
CustomVariable: ZZZ 0.0
CustomVariable: ZZY 0.0
CustomVariable: ZZX 0.0
CustomVariable: ZXY 0.0
CustomVariable: eps_p_q 0.0
```

### **Explanation of Custom Variables**
- **AlphaXX, AlphaYY, AlphaZZ, AlphaZY, AlphaZX, AlphaXY**: These correspond to the **kinematic hardening tensor** (also referred to as **backstress**, $\mathbf{\alpha}$).
- **AlphaInitialXX, AlphaInitialYY, AlphaInitialZZ, AlphaInitialZY, AlphaInitialZX, AlphaInitialXY**: These represent the **initial kinematic hardening tensor** (  $\mathbf{\alpha_{init}}$).
- **ZXX, ZYY, ZZZ, ZZY, ZZX, ZXY**: These correspond to the **fabric tensor** ($\mathbf{z}$), which accounts for anisotropy in the soil structure.

### 2.4 Test Control Parameters

```
AnalysisType: X
NumSteps: 2500
EpsilonYY: -1.0e-4
```

- **AnalysisType: 1** → Undrained Test
- **AnalysisType: 2** → Drained Test
- **AnalysisType: 3** → Constant Mean Stress (p-constant) Monotonic Test
- **AnalysisType: 4** → Cyclic Controlled Stress Test

### 2.5 Cyclic Stress Control (For AnalysisType 4)

```
QUP: 50
QD: -50
```

---

## 3. State Variables

State variables are predefined in **structures.hpp** and track key aspects of the soil state.

```cpp
enum StateVariable {
    StressXX, StressYY, StressZZ, StressZY, StressZX, StressXY,
    DStrainXX, DStrainYY, DStrainZZ, DStrainZY, DStrainZX, DStrainXY,
    PoreWaterPressure, PoreAirPressure, InitialPoreWaterPressure, InitialPoreAirPressure,
    VoidRatio, DegreeOfSaturation,
    PermW_XX, PermW_YY, PermW_ZZ, PermW_ZY, PermW_ZX, PermW_XY,
    PermA_XX, PermA_YY, PermA_ZZ, PermA_ZY, PermA_ZX, PermA_XY,
    Xi, Damping,
    NumVariables
};
```

---

## 4. Custom Variables

Custom variables allow users to define additional state variables specific to a model, such as:
- **Hardening variables** (e.g., backstress components for kinematic hardening).
- **Fabric variables** (e.g., directional anisotropy factors in granular materials).

They offer flexibility to extend the model beyond the predefined state variables.

---

## 5. Troubleshooting

- **File Errors:** Check that `input_data.txt` exists.
- **Parameter Errors:** Ensure that `NumSteps > 0`, `AnalysisType` is 1, 2, 3, or 4, `EpsilonYY` is nonzero, and `QUP > QD` for cyclic tests.
- **Compilation Issues:** Ensure required headers (`SANISANDModelUMAT.hpp`, `UMATBase.hpp`) are available.

---

This manual provides a detailed guide for setting up and running different types of soil behavior simulations using the SANISAND plasticity framework. If you need further clarification, refer to the source code and examples provided in the repository.

