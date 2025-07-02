# Nonlinear HSDT Implementation Guide for cane Multiphysics

## Overview

This document describes the implementation of **Nonlinear Higher-order Shear Deformation Theory (HSDT)** for isogeometric shell analysis in the cane Multiphysics framework. The implementation extends the linear HSDT formulation to include **geometric nonlinearity** through large deformation analysis.

## Table of Contents

1. [Theory and Formulation](#theory-and-formulation)
2. [Implementation Architecture](#implementation-architecture)
3. [Key Functions](#key-functions)
4. [Usage Examples](#usage-examples)
5. [Performance Considerations](#performance-considerations)
6. [Validation and Testing](#validation-and-testing)

## Theory and Formulation

### HSDT vs Kirchhoff-Love Theory

| **Aspect** | **Kirchhoff-Love** | **HSDT** |
|------------|-------------------|----------|
| **DOF per Control Point** | 3 [u, v, w] | 5 [u, v, w, θx, θy] |
| **Transverse Shear** | Neglected (κ = 0) | Included with correction |
| **Shell Thickness** | Thin shells only | Thin and thick shells |
| **Shear Strains** | γξz = γηz = 0 | γξz ≠ 0, γηz ≠ 0 |
| **Rotation DOFs** | Derived from displacements | Independent variables |

### Strain Formulations

#### Membrane Strains (Green-Lagrange)
For nonlinear analysis, we use the Green-Lagrange strain tensor:
```
εξξ = ½(gξξ - Gξξ)
εηη = ½(gηη - Gηη)  
γξη = gξη - Gξη
```
where `g` and `G` are current and reference metric tensors.

#### Bending Strains (Curvature Changes)
```
κξξ = -∂θx/∂ξ
κηη = -∂θy/∂η
κξη = -(∂θx/∂η + ∂θy/∂ξ)
```

#### Shear Strains (Transverse Shear)
```
γξz = ∂w/∂ξ + θx
γηz = ∂w/∂η + θy
```

### Material Matrices

#### Membrane Material Matrix
```
Dm = (Et)/(1-ν²) * [1  ν  0]
                   [ν  1  0]
                   [0  0  (1-ν)/2]
```

#### Bending Material Matrix
```
Db = (Et³)/(12(1-ν²)) * [1  ν  0]
                        [ν  1  0]
                        [0  0  (1-ν)/2]
```

#### Shear Material Matrix
```
Ds = κGt * [1  0]    where κ = 5/6 (shear correction factor)
           [0  1]          G = E/(2(1+ν))
```

### Nonlinear Formulation

#### Tangent Stiffness Matrix
The tangent stiffness includes both material and geometric contributions:
```
KT = KMaterial + KGeometric
```

where:
- `KMaterial`: Material stiffness in current configuration
- `KGeometric`: Geometric stiffness due to stress-displacement coupling

#### Residual Vector
```
R = Finternal - Fexternal
```

where:
- `Finternal = BTσ`: Internal forces from current stress state
- `Fexternal`: Applied external loads

## Implementation Architecture

### File Structure

```
isogeometricThinStructureAnalysis/
├── solutionMatricesAndVectors/
│   ├── computeTangentStiffMtxResVctIGAHSDTShellNLinear.m
│   ├── computeCurrentBaseVectorsHSDT.m
│   ├── computeCurrentStrainsRotationsHSDT.m
│   ├── computeElementTangentStiffnessHSDT.m
│   └── computeElementInternalForcesHSDT.m
├── solvers/
│   └── solve_IGAHSDTShellNLinear.m
└── main/
    └── main_scordelisLoRoof_HSDT_Nonlinear.m
```

### Algorithm Flow

```
Newton-Raphson Iteration:
1. Compute current geometry (deformed configuration)
2. Calculate current strains and rotations
3. Compute stresses from constitutive relations
4. Assemble tangent stiffness matrix
5. Compute residual force vector
6. Solve linearized system: KT Δu = -R
7. Update displacement field: u = u + Δu
8. Check convergence: ||R|| < tolerance
9. If not converged, return to step 1
```

## Key Functions

### 1. Main Nonlinear Solver
**`solve_IGAHSDTShellNLinear.m`**
- Implements Newton-Raphson iteration with load stepping
- **Automatic center point tracking** for load-displacement curves
- Manages convergence criteria and iteration history
- Outputs: displacement field, convergence history, **center displacement tracking**

### 2. Tangent Stiffness and Residual Computation
**`computeTangentStiffMtxResVctIGAHSDTShellNLinear.m`**
- Computes global tangent stiffness matrix and residual vector
- Loops over all elements and Gauss points
- Assembles contributions from membrane, bending, and shear

### 3. Current Configuration Geometry
**`computeCurrentBaseVectorsHSDT.m`**
- Updates base vectors in deformed configuration
- Accounts for displacement gradients: g₁ = G₁ + ∇u
- Essential for geometric nonlinearity

### 4. Current Strains and Rotations
**`computeCurrentStrainsRotationsHSDT.m`**
- Computes Green-Lagrange membrane strains
- Calculates bending strains from rotation gradients
- Evaluates shear strains including displacement and rotation contributions

### 5. Element-Level Computations
**`computeElementTangentStiffnessHSDT.m`** & **`computeElementInternalForcesHSDT.m`**
- Element tangent stiffness: material + geometric contributions
- Element internal forces from virtual work principle

### 6. Load-Displacement Tracking System
**`plotLoadDisplacementHistory.m`**
- **Automatic center point identification** (no manual setup required)
- Real-time displacement tracking during Newton-Raphson iterations
- Comprehensive convergence analysis with visual indicators
- Linear vs nonlinear comparison capabilities
- Stiffness evolution analysis and trend detection

## Usage Examples

### Basic Nonlinear HSDT Analysis with Center Tracking

```matlab
%% Setup nonlinear analysis parameters
propNLinearAnalysis.method = 'newtonRaphson';
propNLinearAnalysis.noLoadSteps = 10;
propNLinearAnalysis.eps = 1e-6;
propNLinearAnalysis.maxIter = 20;

%% Solve nonlinear system with automatic center tracking
[dHat, CPHistory, resHistory, isConverged, BSplinePatch, minElSize, ...
 centerDisplacementHistory, loadHistory] = solve_IGAHSDTShellNLinear...
    (BSplinePatch, propNLinearAnalysis, @solve_LinearSystemMatlabBackslashSolver, ...
     'undefined', graph, 'outputEnabled');

%% Extract results
u = dHat(1:5:end);      % u-displacements
v = dHat(2:5:end);      % v-displacements  
w = dHat(3:5:end);      % w-displacements
theta_x = dHat(4:5:end); % θx-rotations
theta_y = dHat(5:5:end); % θy-rotations

%% Center point analysis (automatically provided)
fprintf('Center displacement: %.3f mm\n', abs(centerDisplacementHistory(end)) * 1000);
fprintf('Load steps converged: %d/%d\n', sum(isConverged), length(isConverged));

%% Generate comprehensive load-displacement analysis
plotLoadDisplacementHistory(centerDisplacementHistory, loadHistory, ...
                           isConverged, BSplinePatch, propNLinearAnalysis, ...
                           true, linearCenterDisplacement);
```

### Scordelis-Lo Roof Analysis

```matlab
% Run the complete analysis
run('main_scordelisLoRoof_HSDT_Nonlinear.m');

% Results automatically include:
% - Linear vs nonlinear comparison
% - Convergence history plots
% - Load-displacement curves
% - Comprehensive result summaries
```

### Material Parameter Setup

```matlab
% Material constants for HSDT
parameters.E = 4.32e8;          % Young's modulus [Pa]
parameters.nue = 0.0;           % Poisson ratio
parameters.t = 0.25;            % Shell thickness [m]
parameters.shearCorrection = 5/6; % Shear correction factor
```

### Boundary Conditions (5 DOF per Control Point)

```matlab
% Example: Fix displacement and rotation at a corner
xiSup = [0 0]; etaSup = [0 0];
for dirSupp = [1 2 3 4 5]  % u, v, w, θx, θy
    homDOFs = findDofs5D_HSDT(homDOFs, xiSup, etaSup, dirSupp, CP);
end
```

## Performance Considerations

### Computational Complexity

| **Aspect** | **Linear HSDT** | **Nonlinear HSDT** |
|------------|----------------|-------------------|
| **DOF System Size** | 5N | 5N |
| **Stiffness Assembly** | Once | Every iteration |
| **Memory Usage** | O(N²) | O(N²) |
| **Solution Time** | O(N) | O(kN) where k = iterations |

### Optimization Strategies

1. **Load Stepping**: Use incremental loading for better convergence
2. **Adaptive Time Stepping**: Adjust load increments based on convergence rate
3. **Preconditioned Solvers**: Use appropriate linear solvers for large systems
4. **Element-Level Optimization**: Vectorize element loop computations

### Memory Requirements

For a mesh with `N` control points:
- **DOF vector**: `5N` elements
- **Stiffness matrix**: `(5N)²` elements (sparse storage recommended)
- **History storage**: `5N × noLoadSteps` for displacement history

## Validation and Testing

### Test Cases

1. **Scordelis-Lo Roof**: Standard benchmark for shell analysis
2. **Cylindrical Shell under Pressure**: Validates shear deformation effects
3. **Spherical Cap under Point Load**: Tests geometric nonlinearity
4. **Cantilever Plate**: Validates bending-membrane coupling

### Verification Metrics

1. **Convergence Rate**: Quadratic convergence for Newton-Raphson
2. **Load-Displacement Curves**: Should show nonlinear stiffening/softening
3. **Energy Conservation**: Check virtual work consistency
4. **Mesh Independence**: Results should converge with refinement

### Expected Results for Scordelis-Lo Roof

- **Linear vs Nonlinear Differences**: Typically 5-15% for moderate loads
- **Convergence**: 3-6 iterations per load step for well-conditioned problems
- **Maximum Displacement**: Should match published benchmarks within 2%

## Troubleshooting Common Issues

### Convergence Problems

1. **Reduce Load Step Size**: Use more increments for better convergence
2. **Check Boundary Conditions**: Ensure no rigid body modes
3. **Material Parameters**: Verify shear correction factor is reasonable
4. **Mesh Quality**: Poor element aspect ratios can cause problems

### Accuracy Issues

1. **Integration Order**: Ensure sufficient Gauss points for nonlinear terms
2. **Element Type**: Higher-order elements often perform better
3. **Geometric Representation**: Accurate NURBS representation is crucial

### Performance Issues

1. **Sparse Matrix Storage**: Use sparse matrices for large problems
2. **Parallel Assembly**: Vectorize element loop operations
3. **Linear Solver Choice**: Consider iterative solvers for very large systems

## Conclusion

The nonlinear HSDT implementation provides a robust framework for large deformation shell analysis with:

- ✅ **Complete 5 DOF formulation** with independent rotations
- ✅ **Geometric nonlinearity** through Green-Lagrange strains
- ✅ **Newton-Raphson solver** with load stepping
- ✅ **Automatic load-displacement tracking** for center point analysis
- ✅ **Real-time convergence monitoring** with comprehensive visualization
- ✅ **Linear vs nonlinear comparison** capabilities
- ✅ **Stiffness evolution analysis** and nonlinearity assessment
- ✅ **Validated implementation** against standard benchmarks
- ✅ **Comprehensive documentation** and examples

### Key Load-Displacement Features:
- **Zero Setup Required**: Automatically identifies center control point
- **Real-Time Tracking**: Monitors displacement during each Newton-Raphson iteration
- **Comprehensive Visualization**: Multi-panel plots showing load-displacement curves, convergence history, and stiffness evolution
- **Intelligent Analysis**: Automatic detection of stiffness softening/hardening trends
- **Comparison Tools**: Direct comparison with linear analysis results
- **Export Capabilities**: All tracking data saved for post-processing

This implementation significantly extends the capabilities of cane Multiphysics for thick shell analysis and large deformation problems while providing immediate visual feedback on nonlinear behavior and maintaining computational efficiency and accuracy.