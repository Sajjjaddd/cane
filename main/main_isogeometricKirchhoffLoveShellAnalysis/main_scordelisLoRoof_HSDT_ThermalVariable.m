%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Script documentation
% 
% Task : The Scordelis-Lo-Roof benchmark problem under pure thermal loading
%        with VARIABLE (nonlinear) temperature through thickness using 
%        nonlinear HSDT (Higher-order Shear Deformation Theory) with 
%        5 DOF per control point
%        
%        Mechanical loads (FAmp) are set to 0 - pure thermal loading
%        Demonstrates thermal bending effects and geometric nonlinearity
%
% Date : [Current Date]
%
%% Preamble
clear;
clc;

%% Includes 

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Include linear equation system solvers
addpath('../../equationSystemSolvers/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the Isogeometric HSDT shell formulation
addpath('../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../isogeometricThinStructureAnalysis/loads/',...
        '../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../isogeometricThinStructureAnalysis/solvers/',...
        '../../isogeometricThinStructureAnalysis/metrics/',...
        '../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../isogeometricThinStructureAnalysis/output/');

%% NURBS parameters

% Global variables
Length = 50;
Radius = 25;

% Polynomial degrees
p = 1;
q = 2;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [-Length/2 -Length/2 -Length/2
             Length/2  Length/2  Length/2];
         
% y-coordinates
CP(:,:,2) = [-Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)
             -Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)];
         
% z-coordinates
CP(:,:,3) = [Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)
             Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)];
       
% Weights
weight = cos(2*pi/9);
CP(:,:,4) = [1 weight 1
             1 weight 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = 0;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
for i= 1:nxi
    for j=1:neta
        if CP(i,j,4)~=1
            isNURBS = 1;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% Material constants

% Young's modulus
parameters.E = 4.32e8;

% Poisson ratio  
parameters.nue = .0;

% Thickness of the shell
parameters.t = .25;

% Density of the shell (used only for dynamics)
parameters.rho = 7850;

% Shear correction factor for HSDT
parameters.shearCorrection = 5/6;

% Thermal properties
parameters.alpha = 12e-6;  % Thermal expansion coefficient [1/K] for steel

%% Thermal loading parameters - VARIABLE THROUGH THICKNESS

% Reference temperature
T0 = 293.15;  % [K] (20°C)

% Temperature distribution through thickness
% Top surface: Hot (heated)
% Bottom surface: Cool (ambient)
T_top = 373.15;     % [K] (100°C) - top surface
T_bottom = 313.15;  % [K] (40°C) - bottom surface
T_avg = (T_top + T_bottom) / 2;  % Average temperature

% Define nonlinear temperature profile function T(z)
% Example: Exponential decay from top to bottom
% T(z) = T_bottom + (T_top - T_bottom) * exp(2*z/t) for z ∈ [-t/2, t/2]
temperatureProfile = @(z) T_bottom + (T_top - T_bottom) * ...
                         exp(2 * z / parameters.t) / exp(1);

% Temperature data structure for variable through-thickness temperature
temperatureData.T0 = T0;
temperatureData.T = T_avg;  % Use average for spatial variation
temperatureData.alpha = parameters.alpha;
temperatureData.throughThickness = 'variable';
temperatureData.T_top = T_top;
temperatureData.T_bottom = T_bottom;
temperatureData.TProfile = temperatureProfile;

% Display temperature profile information
z_plot = linspace(-parameters.t/2, parameters.t/2, 101);
T_plot = arrayfun(temperatureProfile, z_plot);

fprintf('=== Variable Through-Thickness Thermal Loading ===\n');
fprintf('Reference temperature: %.1f°C\n', T0 - 273.15);
fprintf('Top surface temperature: %.1f°C\n', T_top - 273.15);
fprintf('Bottom surface temperature: %.1f°C\n', T_bottom - 273.15);
fprintf('Average temperature: %.1f°C\n', T_avg - 273.15);
fprintf('Temperature gradient: %.1f K/m\n', (T_top - T_bottom) / parameters.t);
fprintf('Thermal expansion coefficient: %.2e /K\n', parameters.alpha);
fprintf('Shell thickness: %.3f m\n', parameters.t);
fprintf('Through-thickness distribution: %s (exponential)\n', temperatureData.throughThickness);

% Plot temperature profile
figure(1);
plot(T_plot - 273.15, z_plot * 1000, 'r-', 'LineWidth', 2);
xlabel('Temperature [°C]');
ylabel('Through-thickness coordinate z [mm]');
title('Nonlinear Temperature Profile Through Shell Thickness');
grid on;
hold on;
plot([T_top T_bottom] - 273.15, [parameters.t/2 -parameters.t/2] * 1000, 'bo', ...
     'MarkerSize', 8, 'MarkerFaceColor', 'blue');
legend('Exponential Profile', 'Top/Bottom Values', 'Location', 'best');
hold off;

%% GUI

% Analysis type
analysis.type = 'isogeometricHSDTShellAnalysis';

% Integration scheme
int.type = 'default';
if strcmp(int.type,'user')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.xetaNGPForLoad = 6;
end

% Equation system solvers
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% On the graphics
graph.index = 2;  % Start from 2 since figure 1 is used for temperature profile
graph.postprocConfig = 'referenceCurrent';
graph.resultant = 'displacement';
graph.component = 'z';

%% Refinement

% Degree by which to elevate
tp = 1; 
tq = 1; 
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
scaling = 1; 
edgeRatio = ceil(Length/Radius/(sin(4*pi/9)));
refXi = edgeRatio*scaling;
refEta = ceil(4/3)*scaling;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions for HSDT (5 DOF per CP)

% supports (Dirichlet boundary conditions)
homDOFs = [];

% For thermal loading with through-thickness gradient, thermal bending occurs
% We need careful boundary conditions to allow thermal expansion and bending

% back and front curved edges: allow thermal expansion and bending
xiSup = [0 0];   etaSup = [0 1];    
for dirSupp = [2]  % constrain v displacement only (allow thermal expansion and bending)
    homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 1];    
for dirSupp = [2]  % constrain v displacement only
    homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Fix minimal constraints to prevent rigid body motion while allowing thermal effects
% One corner: fix u and w (allow thermal expansion at other points)
xiSup = [0 0];   etaSup = [0 0];   
for dirSupp = [1 3]  % constrain u and w at one corner
    homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Fix one rotation to prevent rigid body rotation
xiSup = [0 0];   etaSup = [0 0];   dirSupp = 5;  % constrain θy at corner
homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak Dirichlet boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables.No = 0;

% Mechanical loads (SET TO ZERO for pure thermal loading)
FAmp = 0;  % No mechanical loads

NBC.noCnd = 0;  % No mechanical loads

%% Create the B-Spline patch array for HSDT
BSplinePatch = fillUpPatch_HSDT...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);

%% Compute thermal load vector
fprintf('\n=== Computing Variable Through-Thickness Thermal Load Vector ===\n');

% Initialize thermal load vector
FThermal = zeros(5*BSplinePatch.noCPs,1);

% Thermal load domain (entire structure)
xiThermalExtension = [0 1];
etaThermalExtension = [0 1];
thermalLoadType = 'gradient';

% Compute thermal load vector with variable through-thickness temperature
FThermal = computeThermalLoadVctIGAHSDTShell...
    (FThermal, BSplinePatch, xiThermalExtension, etaThermalExtension, ...
    temperatureData, thermalLoadType, true, 0, BSplinePatch.int, 'outputEnabled');

% Store thermal load vector in patch
BSplinePatch.FThermal = FThermal;

% No mechanical loads (pure thermal)
BSplinePatch.FGamma = zeros(5*BSplinePatch.noCPs,1);

%% Plot reference configuration
figure(graph.index)
plot_referenceConfigurationIGAThinStructure(p,q,Xi,Eta,CP,isNURBS,homDOFs,FThermal,'outputEnabled');
title('Reference configuration for thermal HSDT shell (5 DOF per CP, Variable T)');
graph.index = graph.index + 1;

%% Nonlinear analysis parameters

% Nonlinear analysis method
propNLinearAnalysis.method = 'newtonRaphson';

% Number of load steps for thermal loading (temperature ramping)
propNLinearAnalysis.noLoadSteps = 6;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-6;

% Assign the maximum number of iterations per load step
propNLinearAnalysis.maxIter = 20;

% On the load conservativeness
propNLinearAnalysis.conservativeLoad = true;

%% Linear analysis for comparison
fprintf('\n=== Starting Linear HSDT Variable Thermal Analysis (for comparison) ===\n');
[dHatLinear_HSDT, F_Linear, minElArea_Linear] = solve_IGAHSDTShellLinear...
    (BSplinePatch, solve_LinearSystem, 'outputEnabled');

%% Nonlinear thermal analysis with variable through-thickness temperature
fprintf('\n=== Starting Nonlinear HSDT Variable Thermal Analysis ===\n');
plot_IGANonlinear = 'undefined';  % No plotting during iterations for simplicity

[dHatNonlinear_HSDT, CPHistory, resHistory, isConverged, BSplinePatch_updated, minElSize_NL, ...
 centerDisplacementHistory, loadHistory] = ...
    solve_IGAHSDTShellNLinear...
    (BSplinePatch, propNLinearAnalysis, solve_LinearSystem, ...
    plot_IGANonlinear, graph, 'outputEnabled');

%% Results comparison and analysis

% Extract displacement components
numCPs = nxi * neta;

% Linear results
u_linear = dHatLinear_HSDT(1:5:end);
v_linear = dHatLinear_HSDT(2:5:end);
w_linear = dHatLinear_HSDT(3:5:end);
theta_x_linear = dHatLinear_HSDT(4:5:end);
theta_y_linear = dHatLinear_HSDT(5:5:end);

% Nonlinear results
u_nonlinear = dHatNonlinear_HSDT(1:5:end);
v_nonlinear = dHatNonlinear_HSDT(2:5:end);
w_nonlinear = dHatNonlinear_HSDT(3:5:end);
theta_x_nonlinear = dHatNonlinear_HSDT(4:5:end);
theta_y_nonlinear = dHatNonlinear_HSDT(5:5:end);

%% Display comprehensive thermal analysis results
fprintf('\n=== HSDT Variable Thermal Analysis: Linear vs Nonlinear ===\n');
fprintf('Loading type: Pure thermal with through-thickness gradient\n');
fprintf('Top surface temperature: %.1f°C (%.1f K increase)\n', T_top - 273.15, T_top - T0);
fprintf('Bottom surface temperature: %.1f°C (%.1f K increase)\n', T_bottom - 273.15, T_bottom - T0);
fprintf('Temperature gradient: %.1f K/m\n', (T_top - T_bottom) / parameters.t);
fprintf('Through-thickness distribution: %s (exponential)\n', temperatureData.throughThickness);
fprintf('Total DOFs: %d (5 per control point)\n', length(dHatNonlinear_HSDT));
fprintf('Control points: %d x %d = %d\n', nxi, neta, numCPs);
fprintf('Converged load steps: %d/%d\n', sum(isConverged), propNLinearAnalysis.noLoadSteps);

fprintf('\n--- Linear HSDT Variable Thermal Results ---\n');
fprintf('Max |u|: %.6e (thermal expansion)\n', max(abs(u_linear)));
fprintf('Max |v|: %.6e\n', max(abs(v_linear)));
fprintf('Max |w|: %.6e (thermal bending)\n', max(abs(w_linear)));
fprintf('Max |θx|: %.6e rad (%.3f deg)\n', max(abs(theta_x_linear)), rad2deg(max(abs(theta_x_linear))));
fprintf('Max |θy|: %.6e rad (%.3f deg)\n', max(abs(theta_y_linear)), rad2deg(max(abs(theta_y_linear))));

fprintf('\n--- Nonlinear HSDT Variable Thermal Results ---\n');
fprintf('Max |u|: %.6e (thermal expansion)\n', max(abs(u_nonlinear)));
fprintf('Max |v|: %.6e\n', max(abs(v_nonlinear)));
fprintf('Max |w|: %.6e (thermal bending)\n', max(abs(w_nonlinear)));
fprintf('Max |θx|: %.6e rad (%.3f deg)\n', max(abs(theta_x_nonlinear)), rad2deg(max(abs(theta_x_nonlinear))));
fprintf('Max |θy|: %.6e rad (%.3f deg)\n', max(abs(theta_y_nonlinear)), rad2deg(max(abs(theta_y_nonlinear))));

% Percentage differences
u_diff_percent = 100 * (max(abs(u_nonlinear)) - max(abs(u_linear))) / max(abs(u_linear));
w_diff_percent = 100 * (max(abs(w_nonlinear)) - max(abs(w_linear))) / max(abs(w_linear));
theta_diff_percent = 100 * (max(abs(theta_x_nonlinear)) - max(abs(theta_x_linear))) / max(abs(theta_x_linear));

fprintf('\n--- Variable Thermal Nonlinearity Effects ---\n');
fprintf('Thermal expansion difference (u): %.1f%%\n', u_diff_percent);
fprintf('Thermal bending difference (w): %.1f%%\n', w_diff_percent);
fprintf('Rotation difference (θx): %.1f%%\n', theta_diff_percent);

% Center point analysis
centerCP_xi = ceil(nxi / 2);
centerCP_eta = ceil(neta / 2);
centerCP_index = (centerCP_eta - 1) * nxi + centerCP_xi;

fprintf('\n--- Center Point Variable Thermal Analysis ---\n');
fprintf('Center control point: CP(%d,%d) = Index %d\n', centerCP_xi, centerCP_eta, centerCP_index);
fprintf('Linear center u: %.6e (thermal expansion)\n', u_linear(centerCP_index));
fprintf('Linear center w: %.6e (thermal bending)\n', w_linear(centerCP_index));
fprintf('Nonlinear center u: %.6e (thermal expansion)\n', u_nonlinear(centerCP_index));
fprintf('Nonlinear center w: %.6e (thermal bending)\n', w_nonlinear(centerCP_index));
fprintf('Final center thermal displacement: u = %.3f mm, w = %.3f mm\n', ...
        u_nonlinear(centerCP_index) * 1000, abs(centerDisplacementHistory(end)) * 1000);

% Theoretical thermal effects estimation
avg_thermal_strain = parameters.alpha * (T_avg - T0);
thermal_curvature = parameters.alpha * (T_top - T_bottom) / parameters.t;
theoretical_expansion = avg_thermal_strain * Length/2;
theoretical_bending = thermal_curvature * (Length/2)^2 / 8; % Approximate for simply supported

fprintf('\n--- Variable Thermal Effects Verification ---\n');
fprintf('Average thermal strain: %.6e\n', avg_thermal_strain);
fprintf('Thermal curvature: %.6e /m\n', thermal_curvature);
fprintf('Theoretical expansion: %.6e m (%.3f mm)\n', ...
        theoretical_expansion, theoretical_expansion * 1000);
fprintf('Theoretical bending deflection: %.6e m (%.3f mm)\n', ...
        theoretical_bending, theoretical_bending * 1000);
fprintf('Actual expansion: %.6e m (%.3f mm)\n', ...
        max(abs(u_nonlinear)), max(abs(u_nonlinear)) * 1000);
fprintf('Actual bending deflection: %.6e m (%.3f mm)\n', ...
        max(abs(w_nonlinear)), max(abs(w_nonlinear)) * 1000);

if abs(w_diff_percent) > 5
    fprintf('*** Significant thermal-geometric coupling in bending (>5%% difference) ***\n');
elseif abs(w_diff_percent) > 2
    fprintf('Moderate thermal-geometric coupling in bending (2-5%% difference)\n');
else
    fprintf('Weak thermal-geometric coupling in bending (<2%% difference)\n');
end

%% Plot convergence history
figure(graph.index)
semilogy(1:size(resHistory,1), resHistory, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Residual Norm');
title('Newton-Raphson Convergence History (Variable Thermal Loading)');
grid on;
legend(arrayfun(@(x) sprintf('Temp Step %d', x), 1:propNLinearAnalysis.noLoadSteps, 'UniformOutput', false));
graph.index = graph.index + 1;

%% Variable thermal load-displacement analysis
fprintf('\n=== Generating Variable Thermal Load-Displacement Analysis ===\n');

% Temperature gradient ramping
temperatureGradients = linspace(0, T_top - T_bottom, propNLinearAnalysis.noLoadSteps + 1);
temperatureGradients = temperatureGradients(2:end);  % Remove zero gradient

figure(graph.index)
subplot(2,2,1);
plot(temperatureGradients, abs(centerDisplacementHistory) * 1000, 'o-', ...
     'LineWidth', 2.5, 'MarkerSize', 8, 'Color', 'red');
xlabel('Temperature Gradient [K]');
ylabel('Center Point |w| Displacement [mm]');
title('Thermal Gradient vs Bending Displacement');
grid on;

subplot(2,2,2);
% Plot thermal expansion in u-direction
u_center_history = zeros(length(centerDisplacementHistory), 1);
for i = 1:length(centerDisplacementHistory)
    u_step = CPHistory(1:5:end, i);
    u_center_history(i) = u_step(centerCP_index);
end

plot(temperatureGradients, u_center_history * 1000, 's-', ...
     'LineWidth', 2.5, 'MarkerSize', 8, 'Color', 'blue');
xlabel('Temperature Gradient [K]');
ylabel('Center Point u Displacement [mm]');
title('Thermal Gradient vs Expansion');
grid on;

subplot(2,2,3);
% Plot rotation history
theta_x_center_history = zeros(length(centerDisplacementHistory), 1);
for i = 1:length(centerDisplacementHistory)
    theta_step = CPHistory(4:5:end, i);
    theta_x_center_history(i) = theta_step(centerCP_index);
end

plot(temperatureGradients, rad2deg(abs(theta_x_center_history)), 'd-', ...
     'LineWidth', 2.5, 'MarkerSize', 8, 'Color', 'green');
xlabel('Temperature Gradient [K]');
ylabel('Center Point |θx| [deg]');
title('Thermal Gradient vs Rotation');
grid on;

subplot(2,2,4);
% Plot thermal bending curvature
thermal_curvatures = parameters.alpha * temperatureGradients / parameters.t;
plot(thermal_curvatures, abs(centerDisplacementHistory) * 1000, 'v-', ...
     'LineWidth', 2.5, 'MarkerSize', 8, 'Color', 'magenta');
xlabel('Thermal Curvature [1/m]');
ylabel('Center Point |w| Displacement [mm]');
title('Thermal Curvature vs Bending');
grid on;

sgtitle('Variable Through-Thickness Thermal Analysis (HSDT)', 'FontSize', 14, 'FontWeight', 'bold');
graph.index = graph.index + 1;

%% Save comprehensive variable thermal results
save('scordelisLoRoof_HSDT_ThermalVariable_Results.mat', ...
     'dHatLinear_HSDT', 'dHatNonlinear_HSDT', 'CPHistory', 'resHistory', ...
     'isConverged', 'u_linear', 'v_linear', 'w_linear', 'theta_x_linear', 'theta_y_linear', ...
     'u_nonlinear', 'v_nonlinear', 'w_nonlinear', 'theta_x_nonlinear', 'theta_y_nonlinear', ...
     'centerDisplacementHistory', 'loadHistory', 'centerCP_xi', 'centerCP_eta', ...
     'temperatureData', 'T0', 'T_top', 'T_bottom', 'T_avg', 'temperatureProfile', ...
     'u_center_history', 'theta_x_center_history', 'temperatureGradients', 'thermal_curvatures', ...
     'z_plot', 'T_plot', 'propNLinearAnalysis', 'parameters', 'CP', 'Xi', 'Eta', 'p', 'q');

fprintf('\n=== Variable Thermal Analysis Completed Successfully ===\n');
fprintf('Results saved to: scordelisLoRoof_HSDT_ThermalVariable_Results.mat\n');
fprintf('Pure thermal loading with variable through-thickness temperature distribution\n');
fprintf('Temperature gradient: %.1f K, Top: %.1f°C, Bottom: %.1f°C\n', ...
        T_top - T_bottom, T_top - 273.15, T_bottom - 273.15);
fprintf('Thermal expansion coefficient: %.2e /K\n', parameters.alpha);
fprintf('Maximum thermal expansion: %.3f mm\n', max(abs(u_nonlinear)) * 1000);
fprintf('Maximum thermal bending: %.3f mm\n', max(abs(w_nonlinear)) * 1000);
fprintf('Center displacement tracking completed for variable thermal loading\n');
fprintf('Load steps converged: %d/%d (%.1f%%)\n', sum(isConverged), length(isConverged), ...
        100*sum(isConverged)/length(isConverged));
fprintf('\n');

%% end