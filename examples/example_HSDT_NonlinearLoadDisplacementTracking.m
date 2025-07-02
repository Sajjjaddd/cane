%% Example: Nonlinear HSDT Load-Displacement Tracking
%
% This example demonstrates how to use the enhanced nonlinear HSDT solver
% with automatic load-displacement curve generation for the center point
% of any shell structure.
%
% Features demonstrated:
% - Automatic center point identification
% - Load-displacement tracking during iterations
% - Comprehensive convergence analysis
% - Linear vs nonlinear comparison
% - Stiffness evolution analysis

%% Preamble
clear; clc;

%% Example Parameters (modify for your specific problem)

% Example: Simple plate under pressure
Length = 10;     % [m]
Width = 10;      % [m]  
pressure = -1000; % [N/m²] - distributed load

% Material properties
E = 2.1e11;      % Young's modulus [Pa]
nu = 0.3;        % Poisson ratio
thickness = 0.01; % Shell thickness [m]

%% Setup Analysis (simplified - adapt to your geometry)

% For this example, we'll assume you have already set up:
% - BSplinePatch structure with HSDT formulation
% - Boundary conditions with homDOFs
% - Load vector with NBC structure

% Mock patch structure for demonstration
BSplinePatch.noCPs = 25; % 5x5 control points
BSplinePatch.CP = ones(5,5,4); % Simplified control points
BSplinePatch.DOFNumbering = reshape(1:125, 5, 5, 5); % 5 DOF per CP

% Material parameters for HSDT
BSplinePatch.parameters.E = E;
BSplinePatch.parameters.nue = nu;
BSplinePatch.parameters.t = thickness;
BSplinePatch.parameters.shearCorrection = 5/6;

% Load configuration
BSplinePatch.NBC.loadAmplitude = {pressure};

%% Nonlinear Analysis Parameters

propNLinearAnalysis.method = 'newtonRaphson';
propNLinearAnalysis.noLoadSteps = 8;    % Number of load steps
propNLinearAnalysis.eps = 1e-6;         % Convergence tolerance
propNLinearAnalysis.maxIter = 15;       % Max iterations per load step

%% Linear Analysis for Comparison (optional)

fprintf('=== Running Linear HSDT Analysis for Comparison ===\n');
% [dHatLinear, ~, ~] = solve_IGAHSDTShellLinear(BSplinePatch, @solve_LinearSystemMatlabBackslashSolver, 'outputEnabled');

% For demonstration, create mock linear result
dHatLinear = zeros(125, 1);
dHatLinear(3:5:end) = -0.001 * randn(25, 1); % Mock w-displacements

%% Nonlinear Analysis with Load-Displacement Tracking

fprintf('\n=== Running Nonlinear HSDT Analysis with Center Tracking ===\n');

% The enhanced solver automatically tracks the center point
% [dHat, CPHistory, resHistory, isConverged, BSplinePatch_updated, minElSize, ...
%  centerDisplacementHistory, loadHistory] = solve_IGAHSDTShellNLinear...
%     (BSplinePatch, propNLinearAnalysis, @solve_LinearSystemMatlabBackslashSolver, ...
%      'undefined', struct('index', 1), 'outputEnabled');

% For demonstration, create mock results
dHat = zeros(125, 1);
CPHistory = zeros(125, 8);
resHistory = [1e-3, 1e-6; 2e-3, 1e-7; 1e-3, 1e-8; zeros(12, 2)]'; % Mock convergence
isConverged = true(8, 1);

% Mock center displacement tracking (what the solver now provides automatically)
centerDisplacementHistory = [-0.001, -0.002, -0.003, -0.0045, -0.006, -0.0078, -0.0095, -0.012]';
loadHistory = linspace(0.125, 1.0, 8)';

%% Extract Center Point Information

% The solver automatically identifies the center control point
numCPs_xi = 5;
numCPs_eta = 5;
centerCP_xi = ceil(numCPs_xi / 2);  % = 3
centerCP_eta = ceil(numCPs_eta / 2); % = 3
centerCP_index = (centerCP_eta - 1) * numCPs_xi + centerCP_xi; % = 13

fprintf('\n=== Center Point Analysis ===\n');
fprintf('Identified center control point: CP(%d,%d) = Index %d\n', ...
        centerCP_xi, centerCP_eta, centerCP_index);
fprintf('Linear center displacement: %.6f mm\n', abs(dHatLinear(3 + 5*(centerCP_index-1))) * 1000);
fprintf('Nonlinear final displacement: %.6f mm\n', abs(centerDisplacementHistory(end)) * 1000);

%% Automatic Load-Displacement Curve Generation

fprintf('\n=== Generating Load-Displacement Curves ===\n');

% The solver automatically generates load-displacement curves
% This happens automatically in the solver, but you can also call it manually:

plotLoadDisplacementHistory(centerDisplacementHistory, loadHistory, ...
                           isConverged, BSplinePatch, propNLinearAnalysis, ...
                           true, dHatLinear(3 + 5*(centerCP_index-1)));

%% Analysis of Results

% Convergence analysis
convergenceRate = sum(isConverged) / length(isConverged) * 100;
fprintf('\nConvergence Rate: %.1f%% (%d/%d load steps)\n', ...
        convergenceRate, sum(isConverged), length(isConverged));

% Nonlinearity assessment
linearDisp = abs(dHatLinear(3 + 5*(centerCP_index-1)));
nonlinearDisp = abs(centerDisplacementHistory(end));
nonlinearityFactor = nonlinearDisp / linearDisp;

fprintf('\nNonlinearity Assessment:\n');
fprintf('  Linear displacement: %.6f mm\n', linearDisp * 1000);
fprintf('  Nonlinear displacement: %.6f mm\n', nonlinearDisp * 1000);
fprintf('  Nonlinearity factor: %.3f\n', nonlinearityFactor);

if nonlinearityFactor > 1.1
    fprintf('  *** Significant geometric nonlinearity (stiffness softening) ***\n');
elseif nonlinearityFactor < 0.9
    fprintf('  *** Unusual stiffness hardening detected ***\n');
else
    fprintf('  Mild nonlinear effects\n');
end

% Load-displacement curve characteristics
actualLoads = loadHistory * pressure;
maxLoad = max(abs(actualLoads));
maxDisplacement = max(abs(centerDisplacementHistory)) * 1000;

fprintf('\nLoad-Displacement Curve Characteristics:\n');
fprintf('  Maximum load: %.1f N/m²\n', maxLoad);
fprintf('  Maximum displacement: %.6f mm\n', maxDisplacement);
fprintf('  Load steps: %d\n', length(loadHistory));

%% Usage Instructions

fprintf('\n=== Usage Instructions ===\n');
fprintf('To use this functionality in your own analysis:\n\n');
fprintf('1. Set up your BSplinePatch with HSDT formulation (5 DOF per CP)\n');
fprintf('2. Configure nonlinear analysis parameters:\n');
fprintf('   propNLinearAnalysis.method = ''newtonRaphson'';\n');
fprintf('   propNLinearAnalysis.noLoadSteps = 10;\n');
fprintf('   propNLinearAnalysis.eps = 1e-6;\n');
fprintf('   propNLinearAnalysis.maxIter = 20;\n\n');
fprintf('3. Call the enhanced nonlinear solver:\n');
fprintf('   [dHat, CPHistory, resHistory, isConverged, BSplinePatch, minElSize, ...\n');
fprintf('    centerDisplacementHistory, loadHistory] = solve_IGAHSDTShellNLinear(...);\n\n');
fprintf('4. The solver automatically:\n');
fprintf('   - Identifies the center control point\n');
fprintf('   - Tracks displacement at each converged load step\n');
fprintf('   - Generates load-displacement curves\n');
fprintf('   - Provides convergence analysis\n');
fprintf('   - Compares with linear solution if available\n\n');
fprintf('5. Use plotLoadDisplacementHistory() for additional analysis\n\n');

fprintf('=== Example Complete ===\n');

%% Additional Notes

% The implementation provides:
% 1. Automatic center point identification based on control point grid
% 2. Real-time displacement tracking during Newton-Raphson iterations
% 3. Comprehensive convergence monitoring
% 4. Automatic plot generation with:
%    - Main load-displacement curve
%    - Convergence history
%    - Stiffness evolution analysis
%    - Linear vs nonlinear comparison
% 5. Detailed console output with analysis summaries

% Key benefits:
% - No manual setup required for tracking
% - Works with any HSDT shell geometry
% - Provides immediate visual feedback
% - Enables nonlinearity assessment
% - Facilitates parametric studies