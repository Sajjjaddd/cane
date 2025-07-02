# HSDT Support System Guide
## Enhanced 5-DOF Boundary Condition Support for HSDT Shell Analysis

### Overview

The Enhanced HSDT Support System provides comprehensive boundary condition handling for Higher-order Shear Deformation Theory (HSDT) shell analysis in the cane Multiphysics framework. This system supports all 5 degrees of freedom per control point: **[u, v, w, θx, θy]** including both displacements and rotations.

### Key Features

- ✅ **Full 5-DOF Support**: Complete handling of displacements and rotations
- ✅ **Predefined Support Types**: Common structural support configurations  
- ✅ **Flexible Location Specification**: Points, edges, corners, and regions
- ✅ **Mixed Boundary Conditions**: Supports and prescribed displacements/rotations
- ✅ **Automatic Validation**: Rigid body motion and over-constraint checking
- ✅ **Convenience Functions**: Easy setup for common problems
- ✅ **Comprehensive Documentation**: Detailed examples and usage guides

---

## 1. Core Functions

### 1.1 Main Functions

| Function | Purpose |
|----------|---------|
| `findDofs5D_HSDT_Enhanced()` | Main boundary condition processor |
| `create3DSupports_HSDT()` | Advanced support creation with 5-DOF |
| `createCommonSupports_HSDT()` | Convenience function for common configurations |

### 1.2 Function Signatures

```matlab
% Main boundary condition processor
[homDOFs, inhomDOFs, valuesInhomDOFs, supportInfo] = ...
    findDofs5D_HSDT_Enhanced(BSplinePatch, boundaryConditions, outMsg)

% Advanced support creation
[homDOFs, supportInfo] = create3DSupports_HSDT...
    (homDOFs, supportConfigs, CP, DOFNumbering, supportType, outMsg)

% Common support configurations
boundaryConditions = createCommonSupports_HSDT(supportType, varargin)
```

---

## 2. Supported DOF Types

### 2.1 DOF Ordering
Each control point has 5 degrees of freedom in the following order:
1. **u**: Displacement in local x-direction
2. **v**: Displacement in local y-direction  
3. **w**: Displacement in local z-direction (out-of-plane)
4. **θx**: Rotation about local x-axis
5. **θy**: Rotation about local y-axis

### 2.2 DOF Constraints
Each DOF can be:
- **Free**: Unconstrained (computed by solver)
- **Homogeneous**: Constrained to zero (support condition)
- **Inhomogeneous**: Prescribed to specific value (displacement/rotation)

---

## 3. Support Types

### 3.1 Predefined Support Types

| Support Type | u | v | w | θx | θy | Description |
|--------------|---|---|---|----|----| ------------|
| `fixed` | ✓ | ✓ | ✓ | ✓ | ✓ | All DOFs constrained |
| `pinned` | ✓ | ✓ | ✓ | ✗ | ✗ | Displacements fixed, rotations free |
| `simply_supported` | ✗ | ✗ | ✓ | ✗ | ✗ | Only w constrained |
| `roller_x` | ✗ | ✓ | ✓ | ✗ | ✗ | u free, v,w fixed |
| `roller_y` | ✓ | ✗ | ✓ | ✗ | ✗ | v free, u,w fixed |
| `roller_z` | ✓ | ✓ | ✗ | ✗ | ✗ | w free, u,v fixed |
| `clamped_edge` | ✗ | ✗ | ✓ | ✓ | ✗ | w and θx fixed |
| `symmetry_x` | ✗ | ✓ | ✗ | ✗ | ✓ | Symmetry about x-axis |
| `symmetry_y` | ✓ | ✗ | ✗ | ✓ | ✗ | Symmetry about y-axis |
| `custom` | User-defined constraint combination |

### 3.2 Common Support Configurations

| Configuration | Description |
|---------------|-------------|
| `scordelis_lo` | Standard Scordelis-Lo roof supports (diaphragm edges) |
| `simply_supported_plate` | All edges simply supported |
| `clamped_plate` | All edges fully clamped |
| `cantilever_beam` | One edge fixed, others free |
| `three_point_support` | Minimum support (3 corners) |
| `four_point_support` | All corners pinned |
| `symmetry_quarter` | Quarter model with symmetry BCs |

---

## 4. Location Types

### 4.1 Point Locations
```matlab
location.type = 'point';
location.xi = 0.5;    % Parametric coordinate (0-1)
location.eta = 0.5;   % Parametric coordinate (0-1)
```

### 4.2 Corner Locations
```matlab
location.type = 'corner';
location.corner = 'bottom-left';  % 'bottom-left', 'bottom-right', 
                                  % 'top-left', 'top-right'
```

### 4.3 Edge Locations
```matlab
% Named edges
location.type = 'edge';
location.side = 'bottom';  % 'bottom', 'top', 'left', 'right'

% Parametric edges
location.type = 'edge';
location.xi = [0.2, 0.8];     % xi range
location.eta_fixed = 0.0;     % Fixed eta coordinate
```

### 4.4 Region Locations
```matlab
location.type = 'region';
location.xi = [0.1, 0.9];     % xi range
location.eta = [0.1, 0.9];    % eta range
```

---

## 5. Boundary Condition Types

### 5.1 Support Boundary Conditions
Homogeneous constraints (zero displacement/rotation):

```matlab
bc.type = 'support';
bc.location = location_struct;
bc.supportType = 'fixed';  % Or any predefined type
bc.description = 'Fixed support';

% Custom constraints
bc.supportType = 'custom';
bc.constraints.u = true;      % Constrain u
bc.constraints.v = false;     % Free v
bc.constraints.w = true;      % Constrain w
bc.constraints.thetax = false; % Free θx
bc.constraints.thetay = false; % Free θy
```

### 5.2 Dirichlet Boundary Conditions
Prescribed displacements/rotations:

```matlab
bc.type = 'dirichlet';
bc.location = location_struct;
bc.values.u = 0.01;        % 1cm displacement
bc.values.v = NaN;         % Free (not prescribed)
bc.values.w = -0.005;      % 5mm downward
bc.values.thetax = 0.1;    % 0.1 rad rotation
bc.values.thetay = NaN;    % Free rotation
bc.description = 'Prescribed displacement';
```

---

## 6. Usage Examples

### 6.1 Simple Example: Scordelis-Lo Roof
```matlab
% Create common support configuration
boundaryConditions = createCommonSupports_HSDT('scordelis_lo');

% Process boundary conditions
[homDOFs, inhomDOFs, values, supportInfo] = ...
    findDofs5D_HSDT_Enhanced(BSplinePatch, boundaryConditions, 'outputEnabled');

% Use in solver
[dHat, ~, ~] = solve_IGAHSDTShellLinear(BSplinePatch, F, ...
    homDOFs, inhomDOFs, values, 'outputEnabled');
```

### 6.2 Custom Support Example
```matlab
% Manual support configuration
bc1.type = 'support';
bc1.location.type = 'edge';
bc1.location.side = 'bottom';
bc1.supportType = 'custom';
bc1.constraints.u = true;
bc1.constraints.v = true;
bc1.constraints.w = true;
bc1.constraints.thetax = false;  % Allow rotation
bc1.constraints.thetay = false;
bc1.description = 'Custom edge support';

% Prescribed displacement
bc2.type = 'dirichlet';
bc2.location.type = 'point';
bc2.location.xi = 0.5;
bc2.location.eta = 0.5;
bc2.values.w = -0.01;  % 1cm downward
bc2.values.thetax = 0.05;  % Small rotation
bc2.description = 'Center displacement';

boundaryConditions = [bc1, bc2];
```

### 6.3 Mixed Boundary Conditions
```matlab
% Combine support and prescribed displacements
boundaryConditions = [];

% Fixed corner
bc1.type = 'support';
bc1.location.type = 'corner';
bc1.location.corner = 'bottom-left';
bc1.supportType = 'fixed';

% Prescribed rotation at center
bc2.type = 'dirichlet';
bc2.location.type = 'point';
bc2.location.xi = 0.5;
bc2.location.eta = 0.5;
bc2.values.thetax = 0.1;  % Prescribed rotation
bc2.values.w = NaN;       % Free displacement

boundaryConditions = [bc1, bc2];
```

---

## 7. Integration with Solvers

### 7.1 Linear HSDT Solver
```matlab
[dHat, BSplinePatch, minElSize] = solve_IGAHSDTShellLinear(...
    BSplinePatch, F, homDOFs, inhomDOFs, valuesInhomDOFs, outMsg);
```

### 7.2 Nonlinear HSDT Solver
```matlab
[dHat, CPHistory, resHistory, isConverged, BSplinePatch, minElSize, ...
 centerDisplacementHistory, loadHistory] = solve_IGAHSDTShellNLinear(...
    BSplinePatch, F, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    loadSteps, maxIterations, tolerance, outMsg);
```

### 7.3 Thermal Loading
```matlab
% Thermal boundary conditions work with any support configuration
FThermal = computeThermalLoadVctIGAHSDTShell(...
    BSplinePatch, temperatureFunction, outMsg);

[dHat, ~, ~] = solve_IGAHSDTShellLinear(...
    BSplinePatch, FThermal, homDOFs, inhomDOFs, valuesInhomDOFs, outMsg);
```

---

## 8. Validation and Error Checking

### 8.1 Automatic Validation
The system automatically checks for:
- **Rigid body motion prevention**: Sufficient constraints to prevent free movement
- **Over-constraints**: Too many constraints in small structures
- **DOF compatibility**: Proper 5-DOF structure

### 8.2 Validation Results
```matlab
[~, ~, ~, supportInfo] = findDofs5D_HSDT_Enhanced(...);

if supportInfo.isValid
    fprintf('Constraints are valid\n');
else
    fprintf('Warning: %s\n', supportInfo.validationMessage);
end
```

### 8.3 Support Information
```matlab
% Detailed constraint breakdown
fprintf('Total DOFs: %d\n', 5 * numControlPoints);
fprintf('Constrained DOFs: %d\n', supportInfo.totalConstrainedDOFs);
fprintf('Free DOFs: %d\n', supportInfo.totalFreeDOFs);
fprintf('DOF breakdown: [u:%d, v:%d, w:%d, θx:%d, θy:%d]\n', ...
        supportInfo.constrainedDOFs);
```

---

## 9. Advanced Features

### 9.1 Support Statistics
- Total number of supports
- Constraint breakdown by DOF type
- Free vs constrained DOF counts
- Control point coverage analysis

### 9.2 Constraint Matrix
Access detailed constraint information:
```matlab
constraintMatrix = supportInfo.constraintMatrix;  % [numCPs x 5] logical
constraintValues = supportInfo.constraintValues;   % [numCPs x 5] values
```

### 9.3 Multiple Support Configurations
```matlab
% Combine multiple support types
supports1 = createCommonSupports_HSDT('cantilever_beam');
supports2 = createCommonSupports_HSDT('three_point_support');
allSupports = [supports1, supports2];
```

---

## 10. Best Practices

### 10.1 Minimum Constraint Requirements
For HSDT shells, ensure:
- At least 3 points with w constrained (prevent z-translation)
- At least 2 points with u constrained (prevent x-translation)
- At least 2 points with v constrained (prevent y-translation)
- Some rotational constraints to prevent rigid body rotations

### 10.2 Typical Configurations
- **Scordelis-Lo Roof**: Use `'scordelis_lo'` configuration
- **Plate Problems**: Use `'simply_supported_plate'` or `'clamped_plate'`
- **Beam Problems**: Use `'cantilever_beam'` or `'fixed_edge'`
- **Symmetry Problems**: Use `'symmetry_quarter'` with appropriate reductions

### 10.3 Performance Considerations
- Minimize number of prescribed displacements (inhomogeneous constraints)
- Use symmetry boundary conditions when possible
- Validate constraints before running expensive nonlinear analyses

---

## 11. Troubleshooting

### 11.1 Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| "Insufficient constraints" | Too few supports | Add more constraint points or use stronger support types |
| "Over-constraint" | Too many constraints | Reduce constraints or use weaker support types |
| "DOF numbering error" | Wrong patch structure | Ensure BSplinePatch has 5-DOF DOFNumbering |
| "Solver convergence issues" | Poor constraints | Check constraint validation and support adequacy |

### 11.2 Debug Information
```matlab
% Enable detailed output
[homDOFs, inhomDOFs, values, supportInfo] = ...
    findDofs5D_HSDT_Enhanced(BSplinePatch, boundaryConditions, 'outputEnabled');

% Check constraint matrix
figure; spy(supportInfo.constraintMatrix);
title('Constraint Pattern (CPs vs DOFs)');
xlabel('DOF Type [u,v,w,θx,θy]');
ylabel('Control Point Index');
```

---

## 12. File Structure

### 12.1 Core Files
- `findDofs5D_HSDT_Enhanced.m` - Main boundary condition processor
- `create3DSupports_HSDT.m` - Advanced support creation
- `createCommonSupports_HSDT.m` - Common support configurations

### 12.2 Example Files
- `main_HSDT_SupportSystem_Example.m` - Comprehensive examples
- `main_scordelisLoRoof_HSDT.m` - Updated Scordelis-Lo analysis

### 12.3 Supporting Functions
- `fillUpPatch_HSDT.m` - HSDT patch structure creation
- `solve_IGAHSDTShellLinear.m` - Linear HSDT solver
- `solve_IGAHSDTShellNLinear.m` - Nonlinear HSDT solver

---

## 13. Backward Compatibility

### 13.1 Existing 3-DOF Functions
The new system is designed to complement existing 3-DOF functions:
- `findDofs5D_HSDT.m` - Original 5-DOF function (basic)
- All Kirchhoff-Love functions remain unchanged

### 13.2 Migration Guide
To upgrade from basic to enhanced support system:

```matlab
% Old approach
homDOFs = findDofs5D_HSDT(BSplinePatch, edges, outMsg);

% New approach
boundaryConditions = createCommonSupports_HSDT('simply_supported_plate');
[homDOFs, inhomDOFs, values, supportInfo] = ...
    findDofs5D_HSDT_Enhanced(BSplinePatch, boundaryConditions, outMsg);
```

---

## 14. Future Enhancements

### 14.1 Planned Features
- Neumann boundary condition support (prescribed forces/moments)
- Elastic support stiffness (spring supports)
- Time-dependent boundary conditions
- Multi-patch constraint coupling
- Automatic constraint optimization

### 14.2 Research Directions
- Adaptive constraint placement
- Machine learning for optimal support selection
- Advanced constraint validation algorithms
- Performance optimization for large models

---

## 15. References and Resources

### 15.1 Theoretical Background
- Reddy, J.N. "Mechanics of Laminated Composite Plates and Shells" (HSDT Theory)
- Hughes, T.J.R. "The Finite Element Method" (Isogeometric Analysis)
- Cottrell, J.A. "Isogeometric Analysis: CAD, Finite Elements, NURBS, Exact Geometry"

### 15.2 cane Multiphysics Documentation
- Original cane framework documentation
- Kirchhoff-Love shell implementation
- B-spline and NURBS geometry handling

### 15.3 Example Problems
- Scordelis-Lo roof benchmark
- Simply supported plate solutions
- Cantilever beam validation cases
- Thermal loading examples

---

*This guide provides comprehensive coverage of the Enhanced HSDT Support System. For additional examples and advanced usage, see the example files and test cases included with the implementation.*