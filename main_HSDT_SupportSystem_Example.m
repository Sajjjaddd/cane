%% Comprehensive Example: HSDT 5-DOF Support System
%
% This example demonstrates the enhanced support system for HSDT shell 
% analysis with full 5 DOF support including rotations. Shows various
% support configurations and boundary condition handling.
%
% Author: [Your Name] based on cane Multiphysics framework
% Date: 2024
%
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
%% Initialization

clear; close all; clc;

% Add paths
addpath(genpath('isogeometricThinStructureAnalysis'));
addpath(genpath('BSplineBasisFunctions'));
addpath(genpath('CAGDKernel'));

%% Example 1: Scordelis-Lo Roof with Enhanced Support System

fprintf('=================================================================\n');
fprintf('EXAMPLE 1: Scordelis-Lo Roof with Enhanced HSDT Support System\n');
fprintf('=================================================================\n\n');

%% 1.1 Geometry and Material Properties

% Geometric properties
R = 25.0;           % Radius [m]
L = 50.0;           % Length [m]
theta = 40*pi/180;  % Opening angle [rad]
t = 0.25;           % Thickness [m]

% Material properties
E = 4.32e8;         % Young's modulus [Pa]
nu = 0.0;           % Poisson's ratio
rho = 1.0;          % Density [kg/m³]

% Additional HSDT parameters
shearCorrection = 5/6;  % Shear correction factor

%% 1.2 Create NURBS Surface (Scordelis-Lo roof geometry)

% Control points for cylindrical shell
numCPs_xi = 9;      % Along length
numCPs_eta = 9;     % Along circumferential direction

% Create control point mesh
CP = zeros(numCPs_xi, numCPs_eta, 4);
for i = 1:numCPs_xi
    xi = (i-1)/(numCPs_xi-1);
    x = xi * L;
    for j = 1:numCPs_eta
        eta = (j-1)/(numCPs_eta-1);
        angle = (eta - 0.5) * theta;  % Center around 0
        y = R * sin(angle);
        z = R * cos(angle);
        CP(i,j,:) = [x, y, z, 1.0];
    end
end

% Create knot vectors
p = 2; q = 2;  % Polynomial degrees
Xi = [zeros(1,p+1), linspace(0,1,numCPs_xi-p-1), ones(1,p+1)];
Eta = [zeros(1,q+1), linspace(0,1,numCPs_eta-q-1), ones(1,q+1)];

%% 1.3 Create HSDT B-Spline Patch with 5 DOF

% Create patch structure
BSplinePatch = fillUpPatch_HSDT(CP, p, q, Xi, Eta, [], []);

% Add material properties
BSplinePatch.material.E = E;
BSplinePatch.material.nu = nu;
BSplinePatch.material.rho = rho;
BSplinePatch.material.t = t;
BSplinePatch.material.shearCorrection = shearCorrection;

fprintf('B-Spline patch created:\n');
fprintf('  Control points: %d x %d = %d\n', numCPs_xi, numCPs_eta, numCPs_xi*numCPs_eta);
fprintf('  Total DOFs: %d (5 per control point)\n', 5*numCPs_xi*numCPs_eta);
fprintf('  Polynomial degrees: p=%d, q=%d\n\n', p, q);

%% 1.4 Method A: Using Convenience Function

fprintf('--- Method A: Using createCommonSupports_HSDT convenience function ---\n');

% Create Scordelis-Lo supports using convenience function
boundaryConditions_A = createCommonSupports_HSDT('scordelis_lo');

% Process boundary conditions
[homDOFs_A, inhomDOFs_A, valuesInhomDOFs_A, supportInfo_A] = ...
    findDofs5D_HSDT_Enhanced(BSplinePatch, boundaryConditions_A, 'outputEnabled');

%% 1.5 Method B: Manual Support Configuration

fprintf('\n--- Method B: Manual support configuration ---\n');

% Create boundary conditions manually for comparison
clear boundaryConditions_B;

% Left edge diaphragm support
bc1.type = 'support';
bc1.location.type = 'edge';
bc1.location.side = 'left';
bc1.supportType = 'custom';
bc1.constraints.u = false;      % Free longitudinal movement
bc1.constraints.v = true;       % Constrain transverse movement
bc1.constraints.w = true;       % Constrain vertical movement
bc1.constraints.thetax = false; % Free rotation about x
bc1.constraints.thetay = false; % Free rotation about y
bc1.description = 'Manual left edge diaphragm';

% Right edge diaphragm support
bc2.type = 'support';
bc2.location.type = 'edge';
bc2.location.side = 'right';
bc2.supportType = 'custom';
bc2.constraints.u = false;      % Free longitudinal movement
bc2.constraints.v = true;       % Constrain transverse movement
bc2.constraints.w = true;       % Constrain vertical movement
bc2.constraints.thetax = false; % Free rotation about x
bc2.constraints.thetay = false; % Free rotation about y
bc2.description = 'Manual right edge diaphragm';

boundaryConditions_B = [bc1, bc2];

% Process boundary conditions
[homDOFs_B, inhomDOFs_B, valuesInhomDOFs_B, supportInfo_B] = ...
    findDofs5D_HSDT_Enhanced(BSplinePatch, boundaryConditions_B, 'outputEnabled');

%% 1.6 Verify Both Methods Give Same Result

fprintf('\n--- Verification: Method A vs Method B ---\n');
if isequal(sort(homDOFs_A), sort(homDOFs_B))
    fprintf('✓ Both methods produce identical constraints\n');
else
    fprintf('✗ Methods produce different constraints!\n');
end

%% Example 2: Various Support Types Demonstration

fprintf('\n\n=================================================================\n');
fprintf('EXAMPLE 2: Various Support Types for Rectangular Plate\n');
fprintf('=================================================================\n\n');

%% 2.1 Create Simple Rectangular Plate Geometry

% Create flat rectangular plate
numCPs_xi = 5;
numCPs_eta = 5;
Length = 10.0;
Width = 8.0;

CP_rect = zeros(numCPs_xi, numCPs_eta, 4);
for i = 1:numCPs_xi
    xi = (i-1)/(numCPs_xi-1);
    x = xi * Length;
    for j = 1:numCPs_eta
        eta = (j-1)/(numCPs_eta-1);
        y = eta * Width;
        z = 0.0;
        CP_rect(i,j,:) = [x, y, z, 1.0];
    end
end

% Create patch
Xi_rect = [0 0 0 0.5 1 1 1];
Eta_rect = [0 0 0 0.5 1 1 1];
BSplinePatch_rect = fillUpPatch_HSDT(CP_rect, 2, 2, Xi_rect, Eta_rect, [], []);
BSplinePatch_rect.material = BSplinePatch.material;  % Same material

%% 2.2 Test Different Support Types

supportTypes = {'simply_supported_plate', 'clamped_plate', 'cantilever_beam', ...
                'three_point_support', 'four_point_support', 'symmetry_quarter'};

for iType = 1:length(supportTypes)
    supportType = supportTypes{iType};
    
    fprintf('--- Testing support type: %s ---\n', supportType);
    
    % Create boundary conditions
    boundaryConditions = createCommonSupports_HSDT(supportType);
    
    % Process boundary conditions
    [homDOFs, inhomDOFs, valuesInhomDOFs, supportInfo] = ...
        findDofs5D_HSDT_Enhanced(BSplinePatch_rect, boundaryConditions, 'outputDisabled');
    
    % Display summary
    fprintf('  Total boundary conditions: %d\n', supportInfo.totalBCs);
    fprintf('  Constrained DOFs: %d\n', supportInfo.totalConstrainedDOFs);
    fprintf('  Free DOFs: %d\n', supportInfo.totalFreeDOFs);
    fprintf('  DOF breakdown: [u:%d, v:%d, w:%d, θx:%d, θy:%d]\n', ...
            supportInfo.constrainedDOFs);
    fprintf('  Validation: %s\n', supportInfo.validationMessage);
    fprintf('\n');
end

%% Example 3: Custom Support Configuration with Prescribed Displacements

fprintf('\n=================================================================\n');
fprintf('EXAMPLE 3: Custom Supports with Prescribed Displacements\n');
fprintf('=================================================================\n\n');

%% 3.1 Create Mixed Boundary Conditions

clear boundaryConditions_mixed;

% Support 1: Fixed corner
bc1.type = 'support';
bc1.location.type = 'corner';
bc1.location.corner = 'bottom-left';
bc1.supportType = 'fixed';
bc1.description = 'Fixed corner support';

% Support 2: Pinned edge
bc2.type = 'support';
bc2.location.type = 'edge';
bc2.location.side = 'bottom';
bc2.supportType = 'pinned';
bc2.description = 'Pinned bottom edge';

% Prescribed displacement: Top edge with 5mm downward displacement
bc3.type = 'dirichlet';
bc3.location.type = 'edge';
bc3.location.side = 'top';
bc3.values.u = NaN;        % Free u-displacement
bc3.values.v = NaN;        % Free v-displacement
bc3.values.w = -0.005;     % 5mm downward displacement
bc3.values.thetax = 0.0;   % No x-rotation
bc3.values.thetay = NaN;   % Free y-rotation
bc3.description = 'Prescribed displacement top edge';

% Prescribed rotation: Right edge with 0.1 rad rotation about x
bc4.type = 'dirichlet';
bc4.location.type = 'point';
bc4.location.xi = 1.0;     % Right edge
bc4.location.eta = 0.5;    % Middle point
bc4.values.u = NaN;        % Free u-displacement
bc4.values.v = NaN;        % Free v-displacement
bc4.values.w = NaN;        % Free w-displacement
bc4.values.thetax = 0.1;   % 0.1 rad rotation about x
bc4.values.thetay = NaN;   % Free y-rotation
bc4.description = 'Prescribed rotation point';

boundaryConditions_mixed = [bc1, bc2, bc3, bc4];

%% 3.2 Process Mixed Boundary Conditions

fprintf('--- Processing mixed boundary conditions ---\n');

[homDOFs_mixed, inhomDOFs_mixed, valuesInhomDOFs_mixed, supportInfo_mixed] = ...
    findDofs5D_HSDT_Enhanced(BSplinePatch_rect, boundaryConditions_mixed, 'outputEnabled');

%% 3.3 Display Detailed Results

fprintf('\n--- Detailed constraint analysis ---\n');
fprintf('Homogeneous constraints (zero displacement/rotation): %d DOFs\n', length(homDOFs_mixed));
fprintf('Inhomogeneous constraints (prescribed values): %d DOFs\n', length(inhomDOFs_mixed));

if ~isempty(inhomDOFs_mixed)
    fprintf('\nPrescribed values:\n');
    for i = 1:length(inhomDOFs_mixed)
        fprintf('  DOF %d: %g\n', inhomDOFs_mixed(i), valuesInhomDOFs_mixed(i));
    end
end

%% Example 4: Support Validation and Error Checking

fprintf('\n\n=================================================================\n');
fprintf('EXAMPLE 4: Support Validation and Error Checking\n');
fprintf('=================================================================\n\n');

%% 4.1 Test Insufficient Constraints (Should Fail Validation)

fprintf('--- Testing insufficient constraints (should warn) ---\n');

% Create insufficient support (only one point fixed)
bc_insufficient.type = 'support';
bc_insufficient.location.type = 'corner';
bc_insufficient.location.corner = 'bottom-left';
bc_insufficient.supportType = 'custom';
bc_insufficient.constraints.u = true;
bc_insufficient.constraints.v = false;
bc_insufficient.constraints.w = false;
bc_insufficient.constraints.thetax = false;
bc_insufficient.constraints.thetay = false;
bc_insufficient.description = 'Insufficient constraint test';

[~, ~, ~, supportInfo_insufficient] = ...
    findDofs5D_HSDT_Enhanced(BSplinePatch_rect, bc_insufficient, 'outputEnabled');

%% 4.2 Display Support System Summary

fprintf('\n\n=================================================================\n');
fprintf('SUPPORT SYSTEM SUMMARY\n');
fprintf('=================================================================\n\n');

fprintf('Enhanced HSDT Support System Features:\n\n');

fprintf('✓ 5 DOF Support: [u, v, w, θx, θy] per control point\n');
fprintf('✓ Predefined Support Types:\n');
fprintf('  - fixed: All DOFs constrained\n');
fprintf('  - pinned: Displacements fixed, rotations free\n');
fprintf('  - simply_supported: w fixed, others free\n');
fprintf('  - roller_x/y/z: One displacement free\n');
fprintf('  - symmetry_x/y: Symmetry boundary conditions\n');
fprintf('  - clamped_edge: w and one rotation fixed\n');
fprintf('  - custom: User-defined constraints\n\n');

fprintf('✓ Location Types:\n');
fprintf('  - point: Single control point\n');
fprintf('  - corner: Specific corners (bottom-left, etc.)\n');
fprintf('  - edge: Named edges (bottom, top, left, right)\n');
fprintf('  - region: Rectangular regions\n\n');

fprintf('✓ Boundary Condition Types:\n');
fprintf('  - support: Homogeneous constraints (zero displacement/rotation)\n');
fprintf('  - dirichlet: Prescribed non-zero values\n');
fprintf('  - neumann: Prescribed forces/moments (future)\n\n');

fprintf('✓ Convenience Functions:\n');
fprintf('  - createCommonSupports_HSDT(): Predefined configurations\n');
fprintf('  - create3DSupports_HSDT(): Advanced support creation\n');
fprintf('  - findDofs5D_HSDT_Enhanced(): Comprehensive BC processing\n\n');

fprintf('✓ Validation Features:\n');
fprintf('  - Rigid body motion prevention checks\n');
fprintf('  - Over-constraint detection\n');
fprintf('  - Detailed constraint breakdown\n');
fprintf('  - Support statistics and summaries\n\n');

fprintf('✓ Integration with HSDT Solvers:\n');
fprintf('  - Linear: solve_IGAHSDTShellLinear()\n');
fprintf('  - Nonlinear: solve_IGAHSDTShellNLinear()\n');
fprintf('  - Thermal: Enhanced thermal loading support\n\n');

fprintf('Usage in analysis scripts:\n');
fprintf('  1. Create geometry and material properties\n');
fprintf('  2. Define boundary conditions using convenience functions or manual setup\n');
fprintf('  3. Process with findDofs5D_HSDT_Enhanced()\n');
fprintf('  4. Use resulting DOF arrays in HSDT solvers\n\n');

fprintf('=================================================================\n');
fprintf('EXAMPLES COMPLETED SUCCESSFULLY\n');
fprintf('=================================================================\n');

%% Function Definitions for Support System Usage

function showSupportSystemUsage()
    fprintf('\n=== Quick Usage Guide ===\n\n');
    
    fprintf('1. Common Support Types:\n');
    fprintf('   boundaryConditions = createCommonSupports_HSDT(''scordelis_lo'');\n');
    fprintf('   boundaryConditions = createCommonSupports_HSDT(''simply_supported_plate'');\n\n');
    
    fprintf('2. Custom Support Configuration:\n');
    fprintf('   bc.type = ''support'';\n');
    fprintf('   bc.location.type = ''edge'';\n');
    fprintf('   bc.location.side = ''bottom'';\n');
    fprintf('   bc.supportType = ''fixed'';\n\n');
    
    fprintf('3. Prescribed Displacements:\n');
    fprintf('   bc.type = ''dirichlet'';\n');
    fprintf('   bc.location.type = ''point'';\n');
    fprintf('   bc.location.xi = 0.5; bc.location.eta = 0.5;\n');
    fprintf('   bc.values.w = -0.01;  %% 1cm downward\n\n');
    
    fprintf('4. Process Boundary Conditions:\n');
    fprintf('   [homDOFs, inhomDOFs, values, info] = ...\n');
    fprintf('       findDofs5D_HSDT_Enhanced(BSplinePatch, boundaryConditions, ''outputEnabled'');\n\n');
end

% Call usage guide
showSupportSystemUsage();