%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Script documentation
% 
% Task : The Scordelis-Lo-Roof benchmark problem is modelled and
%        solved using nonlinear HSDT (Higher-order Shear Deformation Theory)
%        with 5 DOF per control point: [u, v, w, θx, θy]
%        
%        This script demonstrates geometric nonlinearity effects
%        and compares with linear HSDT results
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
graph.index = 1;
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

% back and front curved edges are a rigid diaphragm
xiSup = [0 0];   etaSup = [0 1];    
for dirSupp = [2 3]  % constrain v and w displacements
    homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 1];    
for dirSupp = [2 3]  % constrain v and w displacements
    homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Fix the back left corner to avoid rigid body motions
xiSup = [0 0];   etaSup = [0 0];   dirSupp = 1;  % constrain u displacement
homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);

% For nonlinear analysis, might need to constrain some rotations to avoid rigid body rotation
% (Uncomment if needed for stability)
% xiSup = [0 0];   etaSup = [0 0];   
% for dirSupp = [4 5]  % constrain θx and θy rotations at corner
%     homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);
% end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak Dirichlet boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables.No = 0;

% load (Neuman boundary conditions)
% For nonlinear analysis, we can start with higher loads to see geometric effects
FAmp = - 9e1 * 5;  % 5x higher load to see nonlinear effects

NBC.noCnd = 1;
xib = [0 1];   etab = [0 1];   dirForce = 'z';
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.isFollower(1,1) = false;
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAHSDTShell';
NBC.isConservative(1,1) = true;
NBC.isTimeDependent(1,1) = false;

%% Create the B-Spline patch array for HSDT
BSplinePatch = fillUpPatch_HSDT...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);

%% Compute the load vectors for each patch
FGamma = zeros(5*BSplinePatch.noCPs,1);
for counterNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    FGamma = funcHandle...
        (FGamma,BSplinePatch,NBC.xiLoadExtension{counterNBC},...
        NBC.etaLoadExtension{counterNBC},NBC.loadAmplitude{counterNBC},...
        NBC.loadDirection{counterNBC},NBC.isConservative(counterNBC,1),...
        0,BSplinePatch.int,'outputEnabled');
end

% Store load vector in patch
BSplinePatch.FGamma = FGamma;

%% Plot reference configuration
figure(graph.index)
plot_referenceConfigurationIGAThinStructure(p,q,Xi,Eta,CP,isNURBS,homDOFs,FGamma,'outputEnabled');
title('Reference configuration for nonlinear HSDT shell (5 DOF per CP)');
graph.index = graph.index + 1;

%% Nonlinear analysis parameters

% Nonlinear analysis method
propNLinearAnalysis.method = 'newtonRaphson';

% Number of load steps for the nonlinear analysis
propNLinearAnalysis.noLoadSteps = 10;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-6;

% Assign the maximum number of iterations per load step
propNLinearAnalysis.maxIter = 20;

% On the load conservativeness
propNLinearAnalysis.conservativeLoad = true;

%% Linear analysis for comparison
fprintf('=== Starting Linear HSDT Analysis (for comparison) ===\n');
[dHatLinear_HSDT, F_Linear, minElArea_Linear] = solve_IGAHSDTShellLinear...
    (BSplinePatch, solve_LinearSystem, 'outputEnabled');

%% Nonlinear analysis
fprintf('\n=== Starting Nonlinear HSDT Analysis ===\n');
plot_IGANonlinear = 'undefined';  % No plotting during iterations for simplicity

[dHatNonlinear_HSDT, CPHistory, resHistory, isConverged, BSplinePatch_updated, minElSize_NL] = ...
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

%% Display comprehensive results
fprintf('\n=== HSDT Analysis Comparison: Linear vs Nonlinear ===\n');
fprintf('Load factor: %.1f (%.0f%% of linear load)\n', FAmp/(-9e1), 100*FAmp/(-9e1));
fprintf('Total DOFs: %d (5 per control point)\n', length(dHatNonlinear_HSDT));
fprintf('Control points: %d x %d = %d\n', nxi, neta, numCPs);
fprintf('Converged load steps: %d/%d\n', sum(isConverged), propNLinearAnalysis.noLoadSteps);

fprintf('\n--- Linear HSDT Results ---\n');
fprintf('Max |w|: %.6e\n', max(abs(w_linear)));
fprintf('Max |θx|: %.6e rad (%.3f deg)\n', max(abs(theta_x_linear)), rad2deg(max(abs(theta_x_linear))));
fprintf('Max |θy|: %.6e rad (%.3f deg)\n', max(abs(theta_y_linear)), rad2deg(max(abs(theta_y_linear))));

fprintf('\n--- Nonlinear HSDT Results ---\n');
fprintf('Max |w|: %.6e\n', max(abs(w_nonlinear)));
fprintf('Max |θx|: %.6e rad (%.3f deg)\n', max(abs(theta_x_nonlinear)), rad2deg(max(abs(theta_x_nonlinear))));
fprintf('Max |θy|: %.6e rad (%.3f deg)\n', max(abs(theta_y_nonlinear)), rad2deg(max(abs(theta_y_nonlinear))));

% Percentage differences
w_diff_percent = 100 * (max(abs(w_nonlinear)) - max(abs(w_linear))) / max(abs(w_linear));
theta_x_diff_percent = 100 * (max(abs(theta_x_nonlinear)) - max(abs(theta_x_linear))) / max(abs(theta_x_linear));
theta_y_diff_percent = 100 * (max(abs(theta_y_nonlinear)) - max(abs(theta_y_linear))) / max(abs(theta_y_linear));

fprintf('\n--- Nonlinearity Effects ---\n');
fprintf('Displacement difference (w): %.1f%%\n', w_diff_percent);
fprintf('Rotation difference (θx): %.1f%%\n', theta_x_diff_percent);
fprintf('Rotation difference (θy): %.1f%%\n', theta_y_diff_percent);

if abs(w_diff_percent) > 5
    fprintf('*** Significant geometric nonlinearity detected (>5%% difference) ***\n');
else
    fprintf('Geometric nonlinearity effects are moderate (<5%% difference)\n');
end

%% Plot convergence history
figure(graph.index)
semilogy(1:size(resHistory,1), resHistory, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Residual Norm');
title('Newton-Raphson Convergence History');
grid on;
legend(arrayfun(@(x) sprintf('Load Step %d', x), 1:propNLinearAnalysis.noLoadSteps, 'UniformOutput', false));
graph.index = graph.index + 1;

%% Plot load-displacement curve
figure(graph.index)
loadFactors = linspace(1/propNLinearAnalysis.noLoadSteps, 1, propNLinearAnalysis.noLoadSteps);
maxDisplacements = zeros(propNLinearAnalysis.noLoadSteps, 1);
for i = 1:propNLinearAnalysis.noLoadSteps
    w_step = CPHistory(3:5:end, i);
    maxDisplacements(i) = max(abs(w_step));
end

plot(maxDisplacements * 1000, loadFactors * FAmp, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Maximum |w| Displacement [mm]');
ylabel('Applied Load [N/m²]');
title('Load-Displacement Curve (Nonlinear HSDT)');
grid on;
graph.index = graph.index + 1;

%% Save comprehensive results
save('scordelisLoRoof_HSDT_Nonlinear_Results.mat', ...
     'dHatLinear_HSDT', 'dHatNonlinear_HSDT', 'CPHistory', 'resHistory', ...
     'isConverged', 'u_linear', 'v_linear', 'w_linear', 'theta_x_linear', 'theta_y_linear', ...
     'u_nonlinear', 'v_nonlinear', 'w_nonlinear', 'theta_x_nonlinear', 'theta_y_nonlinear', ...
     'propNLinearAnalysis', 'parameters', 'CP', 'Xi', 'Eta', 'p', 'q', 'FAmp');

fprintf('\n=== Analysis Completed Successfully ===\n');
fprintf('Results saved to: scordelisLoRoof_HSDT_Nonlinear_Results.mat\n');
fprintf('Linear and nonlinear HSDT solutions computed and compared\n\n');

%% end