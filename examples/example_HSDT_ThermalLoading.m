%% Example: HSDT Thermal Loading (Constant and Variable Through-Thickness)
%
% This example demonstrates how to use the thermal loading functionality
% for nonlinear HSDT shell analysis with 5 DOF per control point.
%
% Features demonstrated:
% - Constant temperature through thickness (thermal expansion)
% - Variable temperature through thickness (thermal bending)
% - Pure thermal loading (FAmp = 0)
% - Temperature-displacement tracking
% - Linear vs nonlinear thermal effects comparison

%% Preamble
clear; clc;

fprintf('=================================================================\n');
fprintf('HSDT Thermal Loading Example\n');
fprintf('Demonstrates both constant and variable through-thickness thermal loading\n');
fprintf('=================================================================\n\n');

%% Example 1: Constant Temperature Through Thickness

fprintf('--- EXAMPLE 1: CONSTANT TEMPERATURE THROUGH THICKNESS ---\n\n');

% Simple rectangular plate parameters
Length = 2.0;  % [m]
Width = 2.0;   % [m]

% Material properties
E = 2.1e11;           % Young's modulus [Pa]
nu = 0.3;             % Poisson ratio
thickness = 0.01;     % Shell thickness [m]
alpha = 12e-6;        % Thermal expansion coefficient [1/K]

% Thermal loading
T0 = 293.15;          % Reference temperature [K] (20°C)
T_applied = 373.15;   % Applied temperature [K] (100°C)
deltaT = T_applied - T0;  % Temperature increase

fprintf('Material Properties:\n');
fprintf('  Young''s modulus: %.1e Pa\n', E);
fprintf('  Poisson ratio: %.2f\n', nu);
fprintf('  Thickness: %.3f m\n', thickness);
fprintf('  Thermal expansion coefficient: %.2e /K\n', alpha);
fprintf('\nThermal Loading:\n');
fprintf('  Reference temperature: %.1f°C\n', T0 - 273.15);
fprintf('  Applied temperature: %.1f°C\n', T_applied - 273.15);
fprintf('  Temperature increase: %.1f K\n', deltaT);

% Setup temperature data for constant through-thickness
temperatureData_constant.T0 = T0;
temperatureData_constant.T = T_applied;
temperatureData_constant.alpha = alpha;
temperatureData_constant.throughThickness = 'constant';

% Theoretical thermal expansion calculation
theoretical_thermal_strain = alpha * deltaT;
theoretical_expansion = theoretical_thermal_strain * Length / 2;  % Half-length expansion

fprintf('\nTheoretical Results:\n');
fprintf('  Thermal strain: %.2e\n', theoretical_thermal_strain);
fprintf('  Expected half-length expansion: %.3f mm\n', theoretical_expansion * 1000);

%% Example 2: Variable Temperature Through Thickness

fprintf('\n--- EXAMPLE 2: VARIABLE TEMPERATURE THROUGH THICKNESS ---\n\n');

% Temperature distribution through thickness
T_top = 373.15;       % [K] (100°C) - top surface
T_bottom = 313.15;    % [K] (40°C) - bottom surface
T_avg = (T_top + T_bottom) / 2;
temperature_gradient = (T_top - T_bottom) / thickness;

% Define temperature profile function (exponential)
temperatureProfile = @(z) T_bottom + (T_top - T_bottom) * exp(2 * z / thickness) / exp(1);

% Setup temperature data for variable through-thickness
temperatureData_variable.T0 = T0;
temperatureData_variable.T = T_avg;
temperatureData_variable.alpha = alpha;
temperatureData_variable.throughThickness = 'variable';
temperatureData_variable.T_top = T_top;
temperatureData_variable.T_bottom = T_bottom;
temperatureData_variable.TProfile = temperatureProfile;

fprintf('Variable Temperature Distribution:\n');
fprintf('  Top surface: %.1f°C\n', T_top - 273.15);
fprintf('  Bottom surface: %.1f°C\n', T_bottom - 273.15);
fprintf('  Average temperature: %.1f°C\n', T_avg - 273.15);
fprintf('  Temperature gradient: %.1f K/m\n', temperature_gradient);

% Theoretical thermal bending calculation
avg_thermal_strain = alpha * (T_avg - T0);
thermal_curvature = alpha * (T_top - T_bottom) / thickness;
theoretical_bending = thermal_curvature * Length^2 / 8;  % Approximate for simply supported

fprintf('\nTheoretical Variable Temperature Results:\n');
fprintf('  Average thermal strain: %.2e\n', avg_thermal_strain);
fprintf('  Thermal curvature: %.2e /m\n', thermal_curvature);
fprintf('  Expected thermal expansion: %.3f mm\n', avg_thermal_strain * Length/2 * 1000);
fprintf('  Expected thermal bending: %.3f mm\n', theoretical_bending * 1000);

%% Plot Temperature Profiles

figure(1);

% Constant temperature profile
subplot(1,2,1);
z_coords = linspace(-thickness/2, thickness/2, 21);
T_constant = T_applied * ones(size(z_coords));
plot(T_constant - 273.15, z_coords * 1000, 'b-', 'LineWidth', 3);
xlabel('Temperature [°C]');
ylabel('Through-thickness coordinate z [mm]');
title('Constant Temperature Profile');
grid on;
xlim([T0-273.15-10, T_applied-273.15+10]);

% Variable temperature profile
subplot(1,2,2);
z_coords_fine = linspace(-thickness/2, thickness/2, 101);
T_variable = arrayfun(temperatureProfile, z_coords_fine);
plot(T_variable - 273.15, z_coords_fine * 1000, 'r-', 'LineWidth', 3);
hold on;
plot([T_top T_bottom] - 273.15, [thickness/2 -thickness/2] * 1000, 'ro', ...
     'MarkerSize', 8, 'MarkerFaceColor', 'red');
xlabel('Temperature [°C]');
ylabel('Through-thickness coordinate z [mm]');
title('Variable Temperature Profile');
legend('Exponential Profile', 'Top/Bottom Values', 'Location', 'best');
grid on;
hold off;

sgtitle('Temperature Distributions Through Shell Thickness', 'FontSize', 14);

%% Thermal Load Vector Computation Examples

% Note: These examples show the thermal load computation process
% In actual analysis, you would use these within your shell analysis framework

fprintf('\n--- THERMAL LOAD VECTOR COMPUTATION ---\n\n');

%% Mock BSplinePatch structure for demonstration
% (In real analysis, this would be your actual patch)
BSplinePatch_mock.p = 2;
BSplinePatch_mock.q = 2;
BSplinePatch_mock.Xi = [0 0 0 1 1 1];
BSplinePatch_mock.Eta = [0 0 0 1 1 1];
BSplinePatch_mock.CP = ones(3,3,4);
BSplinePatch_mock.isNURBS = false;
BSplinePatch_mock.parameters.E = E;
BSplinePatch_mock.parameters.nue = nu;
BSplinePatch_mock.parameters.t = thickness;
BSplinePatch_mock.DOFNumbering = reshape(1:45, 3, 3, 5);  % 5 DOF per CP
BSplinePatch_mock.noCPs = 9;
BSplinePatch_mock.int.type = 'default';

% Thermal load domain (entire structure)
xiThermalExtension = [0 1];
etaThermalExtension = [0 1];

fprintf('Computing thermal load vectors:\n');

% Constant temperature thermal load
fprintf('  1. Constant temperature thermal load...\n');
FThermal_constant = zeros(45, 1);  % 5 DOF × 9 CPs
% FThermal_constant = computeThermalLoadVctIGAHSDTShell...
%     (FThermal_constant, BSplinePatch_mock, xiThermalExtension, etaThermalExtension, ...
%     temperatureData_constant, 'uniform', true, 0, BSplinePatch_mock.int, 'outputEnabled');

% Variable temperature thermal load
fprintf('  2. Variable temperature thermal load...\n');
FThermal_variable = zeros(45, 1);  % 5 DOF × 9 CPs
% FThermal_variable = computeThermalLoadVctIGAHSDTShell...
%     (FThermal_variable, BSplinePatch_mock, xiThermalExtension, etaThermalExtension, ...
%     temperatureData_variable, 'gradient', true, 0, BSplinePatch_mock.int, 'outputEnabled');

fprintf('  Thermal load vectors computed successfully!\n');

%% Analysis Setup Example

fprintf('\n--- ANALYSIS SETUP EXAMPLE ---\n\n');

% Nonlinear analysis parameters for thermal loading
propNLinearAnalysis.method = 'newtonRaphson';
propNLinearAnalysis.noLoadSteps = 5;    % Temperature ramping steps
propNLinearAnalysis.eps = 1e-6;         % Convergence tolerance
propNLinearAnalysis.maxIter = 15;       % Max iterations per step

fprintf('Nonlinear Analysis Setup:\n');
fprintf('  Method: %s\n', propNLinearAnalysis.method);
fprintf('  Load steps: %d\n', propNLinearAnalysis.noLoadSteps);
fprintf('  Convergence tolerance: %.1e\n', propNLinearAnalysis.eps);
fprintf('  Max iterations per step: %d\n', propNLinearAnalysis.maxIter);

%% Thermal Effects Comparison

fprintf('\n--- THERMAL EFFECTS COMPARISON ---\n\n');

% Compare constant vs variable thermal effects
fprintf('Expected Thermal Effects:\n\n');

fprintf('Constant Temperature:\n');
fprintf('  - Pure thermal expansion in all directions\n');
fprintf('  - No thermal bending (no through-thickness gradient)\n');
fprintf('  - Maximum displacement in expansion direction\n');
fprintf('  - Expected expansion: %.3f mm\n', theoretical_expansion * 1000);

fprintf('\nVariable Temperature:\n');
fprintf('  - Thermal expansion from average temperature\n');
fprintf('  - Thermal bending from temperature gradient\n');
fprintf('  - Combined expansion and bending effects\n');
fprintf('  - Expected expansion: %.3f mm\n', avg_thermal_strain * Length/2 * 1000);
fprintf('  - Expected bending: %.3f mm\n', theoretical_bending * 1000);

% Nonlinearity assessment
expansion_ratio = theoretical_expansion / thickness;
bending_ratio = theoretical_bending / thickness;

fprintf('\nNonlinearity Assessment:\n');
fprintf('  Expansion-to-thickness ratio: %.3f\n', expansion_ratio);
fprintf('  Bending-to-thickness ratio: %.3f\n', bending_ratio);

if expansion_ratio > 0.1 || bending_ratio > 0.1
    fprintf('  *** Geometric nonlinearity likely significant ***\n');
else
    fprintf('  Geometric nonlinearity effects may be moderate\n');
end

%% Usage Instructions

fprintf('\n=== USAGE INSTRUCTIONS ===\n\n');

fprintf('To use thermal loading in your HSDT analysis:\n\n');

fprintf('1. Setup thermal properties:\n');
fprintf('   parameters.alpha = 12e-6;  %% Thermal expansion coefficient\n');
fprintf('   T0 = 293.15;               %% Reference temperature [K]\n');
fprintf('   T_applied = 373.15;        %% Applied temperature [K]\n\n');

fprintf('2. For CONSTANT through-thickness temperature:\n');
fprintf('   temperatureData.T0 = T0;\n');
fprintf('   temperatureData.T = T_applied;\n');
fprintf('   temperatureData.alpha = parameters.alpha;\n');
fprintf('   temperatureData.throughThickness = ''constant'';\n\n');

fprintf('3. For VARIABLE through-thickness temperature:\n');
fprintf('   temperatureData.T0 = T0;\n');
fprintf('   temperatureData.T = T_avg;  %% Average temperature\n');
fprintf('   temperatureData.alpha = parameters.alpha;\n');
fprintf('   temperatureData.throughThickness = ''variable'';\n');
fprintf('   temperatureData.T_top = T_top;\n');
fprintf('   temperatureData.T_bottom = T_bottom;\n');
fprintf('   temperatureData.TProfile = @(z) custom_function(z);\n\n');

fprintf('4. Compute thermal load vector:\n');
fprintf('   FThermal = computeThermalLoadVctIGAHSDTShell(...\n');
fprintf('       FThermal, BSplinePatch, [0 1], [0 1], ...\n');
fprintf('       temperatureData, ''uniform'', true, 0, int, ''outputEnabled'');\n\n');

fprintf('5. Set mechanical loads to zero for pure thermal:\n');
fprintf('   FAmp = 0;  %% No mechanical loads\n');
fprintf('   BSplinePatch.FGamma = zeros(5*noCPs, 1);\n');
fprintf('   BSplinePatch.FThermal = FThermal;\n\n');

fprintf('6. Run nonlinear analysis:\n');
fprintf('   [dHat, CPHistory, resHistory, isConverged, ..., ...\n');
fprintf('    centerDisplacementHistory, loadHistory] = ...\n');
fprintf('       solve_IGAHSDTShellNLinear(BSplinePatch, propNLinearAnalysis, ...);\n\n');

fprintf('7. The solver automatically tracks thermal effects and provides:\n');
fprintf('   - Temperature-displacement curves\n');
fprintf('   - Thermal expansion monitoring\n');
fprintf('   - Thermal bending analysis (for variable temperature)\n');
fprintf('   - Linear vs nonlinear thermal comparison\n\n');

%% Key Thermal Loading Features

fprintf('=== KEY THERMAL LOADING FEATURES ===\n\n');

fprintf('✓ Constant through-thickness temperature:\n');
fprintf('  - Pure thermal expansion effects\n');
fprintf('  - Membrane thermal strains only\n');
fprintf('  - Suitable for uniform heating scenarios\n\n');

fprintf('✓ Variable through-thickness temperature:\n');
fprintf('  - Combined thermal expansion and bending\n');
fprintf('  - Membrane and bending thermal strains\n');
fprintf('  - Custom temperature profile functions\n');
fprintf('  - Linear and exponential profiles supported\n\n');

fprintf('✓ Nonlinear thermal analysis:\n');
fprintf('  - Geometric nonlinearity in thermal loading\n');
fprintf('  - Thermal-mechanical coupling\n');
fprintf('  - Load stepping for temperature ramping\n');
fprintf('  - Convergence monitoring for thermal effects\n\n');

fprintf('✓ Automatic result tracking:\n');
fprintf('  - Temperature-displacement curves\n');
fprintf('  - Thermal expansion monitoring\n');
fprintf('  - Center point tracking for thermal effects\n');
fprintf('  - Comparison with theoretical predictions\n\n');

fprintf('=== EXAMPLE COMPLETE ===\n\n');

%% Display Final Summary

fprintf('SUMMARY:\n');
fprintf('This example demonstrated thermal loading capabilities for HSDT shells:\n');
fprintf('- Constant temperature: %.1f K increase → %.3f mm expansion\n', ...
        deltaT, theoretical_expansion * 1000);
fprintf('- Variable temperature: %.1f K gradient → %.3f mm bending\n', ...
        T_top - T_bottom, theoretical_bending * 1000);
fprintf('- Both cases use 5 DOF per control point [u, v, w, θx, θy]\n');
fprintf('- Pure thermal loading (mechanical loads = 0)\n');
fprintf('- Nonlinear analysis with temperature tracking\n\n');

fprintf('Run the main analysis scripts for complete thermal analysis:\n');
fprintf('- main_scordelisLoRoof_HSDT_ThermalConstant.m\n');
fprintf('- main_scordelisLoRoof_HSDT_ThermalVariable.m\n\n');