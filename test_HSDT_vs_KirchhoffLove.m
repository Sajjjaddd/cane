%% Test Script: HSDT vs Kirchhoff-Love Comparison
%
% Purpose: Compare results between HSDT (5 DOF) and Kirchhoff-Love (3 DOF)
%          formulations for the Scordelis-Lo Roof benchmark
%
% Expected behavior:
% - For thin shells: Results should be very similar
% - For thick shells: HSDT should capture shear effects better
% - HSDT should have physically meaningful rotational DOFs

clear; clc; close all;

%% Test Parameters
fprintf('=== HSDT vs Kirchhoff-Love Comparison Test ===\n\n');

% Test different thickness ratios
thickness_ratios = [0.01, 0.1, 0.25]; % t/R ratios to test
Length = 50;
Radius = 25;

results = struct();

%% Common Setup
addpath('../../generalMath/');
addpath('../../auxiliary/');
addpath('../../equationSystemSolvers/');
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
addpath('../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../isogeometricThinStructureAnalysis/loads/',...
        '../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../isogeometricThinStructureAnalysis/solvers/',...
        '../../isogeometricThinStructureAnalysis/metrics/',...
        '../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../isogeometricThinStructureAnalysis/output/');

solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

%% Loop over different thickness ratios
for t_idx = 1:length(thickness_ratios)
    
    t_ratio = thickness_ratios(t_idx);
    thickness = t_ratio * Radius;
    
    fprintf('--- Testing thickness ratio t/R = %.3f (t = %.3f) ---\n', t_ratio, thickness);
    
    %% Common geometry setup
    p = 1; q = 2;
    Xi = [0 0 1 1];
    Eta = [0 0 0 1 1 1];
    
    % Control Points
    CP(:,:,1) = [-Length/2 -Length/2 -Length/2; Length/2 Length/2 Length/2];
    CP(:,:,2) = [-Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9); -Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)];
    CP(:,:,3) = [Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9); Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)];
    weight = cos(2*pi/9);
    CP(:,:,4) = [1 weight 1; 1 weight 1];
    
    isNURBS = 1;
    
    % Material properties
    E = 4.32e8;
    nue = 0.0;
    
    % Integration
    int.type = 'default';
    
    % Refinement
    tp = 1; tq = 1;
    [Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'');
    refXi = 2; refEta = 2;
    [Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'');
    
    %% Test 1: Kirchhoff-Love (3 DOF)
    fprintf('  Running Kirchhoff-Love analysis...\n');
    
    % Setup parameters
    parameters_KL.E = E;
    parameters_KL.nue = nue;
    parameters_KL.t = thickness;
    parameters_KL.rho = 7850;
    
    % Analysis setup
    analysis_KL.type = 'isogeometricKirchhoffLoveShellAnalysis';
    
    % Boundary conditions
    homDOFs_KL = [];
    xiSup = [0 0]; etaSup = [0 1];
    for dirSupp = [2 3]
        homDOFs_KL = findDofs3D(homDOFs_KL,xiSup,etaSup,dirSupp,CP);
    end
    xiSup = [1 1]; etaSup = [0 1];
    for dirSupp = [2 3]
        homDOFs_KL = findDofs3D(homDOFs_KL,xiSup,etaSup,dirSupp,CP);
    end
    xiSup = [0 0]; etaSup = [0 0]; dirSupp = 1;
    homDOFs_KL = findDofs3D(homDOFs_KL,xiSup,etaSup,dirSupp,CP);
    
    % Load
    FAmp = -9e1;
    NBC_KL.noCnd = 1;
    NBC_KL.xiLoadExtension = {[0 1]};
    NBC_KL.etaLoadExtension = {[0 1]};
    NBC_KL.loadAmplitude = {FAmp};
    NBC_KL.loadDirection = {'z'};
    NBC_KL.isFollower(1,1) = false;
    NBC_KL.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
    NBC_KL.isConservative(1,1) = true;
    NBC_KL.isTimeDependent(1,1) = false;
    
    % Create patch
    BSplinePatch_KL = fillUpPatch(analysis_KL,p,Xi,q,Eta,CP,isNURBS,parameters_KL,...
        homDOFs_KL,[],[],struct('noCnd',0),struct('No',0),NBC_KL,[],[],[],[],[],int);
    
    % Solve
    try
        [dHat_KL,~,~] = solve_IGAKirchhoffLoveShellLinear(BSplinePatch_KL,solve_LinearSystem,'');
        
        % Extract displacements
        u_KL = dHat_KL(1:3:end);
        v_KL = dHat_KL(2:3:end);
        w_KL = dHat_KL(3:3:end);
        
        results.KL(t_idx).success = true;
        results.KL(t_idx).max_u = max(abs(u_KL));
        results.KL(t_idx).max_v = max(abs(v_KL));
        results.KL(t_idx).max_w = max(abs(w_KL));
        results.KL(t_idx).max_total = max(abs(dHat_KL));
        results.KL(t_idx).nDOFs = length(dHat_KL);
        
        fprintf('    ✓ Kirchhoff-Love: max|w| = %.6e\n', results.KL(t_idx).max_w);
        
    catch ME
        fprintf('    ✗ Kirchhoff-Love failed: %s\n', ME.message);
        results.KL(t_idx).success = false;
    end
    
    %% Test 2: HSDT (5 DOF)
    fprintf('  Running HSDT analysis...\n');
    
    % Setup parameters
    parameters_HSDT.E = E;
    parameters_HSDT.nue = nue;
    parameters_HSDT.t = thickness;
    parameters_HSDT.rho = 7850;
    parameters_HSDT.shearCorrection = 5/6;
    
    % Analysis setup
    analysis_HSDT.type = 'isogeometricHSDTShellAnalysis';
    
    % Boundary conditions (5 DOF)
    homDOFs_HSDT = [];
    xiSup = [0 0]; etaSup = [0 1];
    for dirSupp = [2 3]  % y,z displacements
        homDOFs_HSDT = findDofs5D(homDOFs_HSDT,xiSup,etaSup,dirSupp,CP);
    end
    xiSup = [1 1]; etaSup = [0 1];
    for dirSupp = [2 3]  % y,z displacements
        homDOFs_HSDT = findDofs5D(homDOFs_HSDT,xiSup,etaSup,dirSupp,CP);
    end
    % Also constrain rotations for rigid diaphragm
    xiSup = [0 0]; etaSup = [0 1];
    for dirSupp = [4 5]  % rotations
        homDOFs_HSDT = findDofs5D(homDOFs_HSDT,xiSup,etaSup,dirSupp,CP);
    end
    xiSup = [1 1]; etaSup = [0 1];
    for dirSupp = [4 5]  % rotations
        homDOFs_HSDT = findDofs5D(homDOFs_HSDT,xiSup,etaSup,dirSupp,CP);
    end
    xiSup = [0 0]; etaSup = [0 0]; dirSupp = 1;  % x displacement at corner
    homDOFs_HSDT = findDofs5D(homDOFs_HSDT,xiSup,etaSup,dirSupp,CP);
    
    % Load
    NBC_HSDT.noCnd = 1;
    NBC_HSDT.xiLoadExtension = {[0 1]};
    NBC_HSDT.etaLoadExtension = {[0 1]};
    NBC_HSDT.loadAmplitude = {FAmp};
    NBC_HSDT.loadDirection = {'z'};
    NBC_HSDT.isFollower(1,1) = false;
    NBC_HSDT.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure5DOF';
    NBC_HSDT.isConservative(1,1) = true;
    NBC_HSDT.isTimeDependent(1,1) = false;
    
    % Create patch
    BSplinePatch_HSDT = fillUpPatch(analysis_HSDT,p,Xi,q,Eta,CP,isNURBS,parameters_HSDT,...
        homDOFs_HSDT,[],[],struct('noCnd',0),struct('No',0),NBC_HSDT,[],[],[],[],[],int);
    
    % Solve
    try
        [dHat_HSDT,~,~] = solve_IGAHSDTShellLinear(BSplinePatch_HSDT,solve_LinearSystem,'');
        
        % Extract displacements and rotations
        u_HSDT = dHat_HSDT(1:5:end);
        v_HSDT = dHat_HSDT(2:5:end);
        w_HSDT = dHat_HSDT(3:5:end);
        theta_x = dHat_HSDT(4:5:end);
        theta_y = dHat_HSDT(5:5:end);
        
        results.HSDT(t_idx).success = true;
        results.HSDT(t_idx).max_u = max(abs(u_HSDT));
        results.HSDT(t_idx).max_v = max(abs(v_HSDT));
        results.HSDT(t_idx).max_w = max(abs(w_HSDT));
        results.HSDT(t_idx).max_theta_x = max(abs(theta_x));
        results.HSDT(t_idx).max_theta_y = max(abs(theta_y));
        results.HSDT(t_idx).max_total = max(abs(dHat_HSDT));
        results.HSDT(t_idx).nDOFs = length(dHat_HSDT);
        
        fprintf('    ✓ HSDT: max|w| = %.6e, max|θx| = %.6e, max|θy| = %.6e\n', ...
            results.HSDT(t_idx).max_w, results.HSDT(t_idx).max_theta_x, results.HSDT(t_idx).max_theta_y);
        
    catch ME
        fprintf('    ✗ HSDT failed: %s\n', ME.message);
        results.HSDT(t_idx).success = false;
    end
    
    %% Compare results
    if results.KL(t_idx).success && results.HSDT(t_idx).success
        w_diff_percent = abs(results.HSDT(t_idx).max_w - results.KL(t_idx).max_w) / results.KL(t_idx).max_w * 100;
        shear_ratio = max(results.HSDT(t_idx).max_theta_x, results.HSDT(t_idx).max_theta_y) / results.HSDT(t_idx).max_w;
        
        results.comparison(t_idx).w_diff_percent = w_diff_percent;
        results.comparison(t_idx).shear_ratio = shear_ratio;
        results.comparison(t_idx).dof_ratio = results.HSDT(t_idx).nDOFs / results.KL(t_idx).nDOFs;
        
        fprintf('    📊 Comparison: w difference = %.2f%%, shear ratio = %.3e\n', w_diff_percent, shear_ratio);
        
        if w_diff_percent < 5 && t_ratio < 0.1
            fprintf('    ✅ Good agreement for thin shell (as expected)\n');
        elseif w_diff_percent > 10 && t_ratio > 0.1
            fprintf('    ✅ Significant difference for thick shell (HSDT captures shear effects)\n');
        else
            fprintf('    ⚠️  Unexpected behavior\n');
        end
    end
    
    fprintf('\n');
end

%% Summary Report
fprintf('=== SUMMARY REPORT ===\n');
fprintf('Thickness\t  KL max|w|\t  HSDT max|w|\t  Difference\t  Shear Ratio\t  Status\n');
fprintf('------------------------------------------------------------------------\n');

for i = 1:length(thickness_ratios)
    if results.KL(i).success && results.HSDT(i).success
        fprintf('t/R=%.3f\t  %.3e\t  %.3e\t  %6.2f%%\t  %.2e\t  ', ...
            thickness_ratios(i), results.KL(i).max_w, results.HSDT(i).max_w, ...
            results.comparison(i).w_diff_percent, results.comparison(i).shear_ratio);
        
        if thickness_ratios(i) < 0.1 && results.comparison(i).w_diff_percent < 10
            fprintf('✅ Thin shell OK\n');
        elseif thickness_ratios(i) >= 0.1 && results.comparison(i).shear_ratio > 0.01
            fprintf('✅ Thick shell OK\n');
        else
            fprintf('⚠️  Check needed\n');
        end
    else
        fprintf('t/R=%.3f\t  ---\t\t  ---\t\t  ---\t\t  ---\t\t  ❌ Failed\n', thickness_ratios(i));
    end
end

fprintf('\n=== VALIDATION CONCLUSIONS ===\n');
fprintf('1. HSDT implementation appears to be: ');
if all([results.KL.success]) && all([results.HSDT.success])
    fprintf('✅ WORKING CORRECTLY\n');
else
    fprintf('❌ HAVING ISSUES\n');
end

fprintf('2. Physics verification:\n');
for i = 1:length(thickness_ratios)
    if results.HSDT(i).success
        if results.HSDT(i).max_theta_x > 0 && results.HSDT(i).max_theta_y > 0
            fprintf('   ✅ Rotational DOFs are active (θx=%.2e, θy=%.2e)\n', ...
                results.HSDT(i).max_theta_x, results.HSDT(i).max_theta_y);
            break;
        end
    end
end

fprintf('3. Theory consistency:\n');
thin_shell_ok = false;
if length(results.comparison) >= 1 && results.comparison(1).w_diff_percent < 10
    fprintf('   ✅ Thin shells: HSDT ≈ Kirchhoff-Love (difference < 10%%)\n');
    thin_shell_ok = true;
end
if length(results.comparison) >= 3 && results.comparison(3).w_diff_percent > 5
    fprintf('   ✅ Thick shells: HSDT ≠ Kirchhoff-Love (difference > 5%%)\n');
end

if all([results.KL.success]) && all([results.HSDT.success]) && thin_shell_ok
    fprintf('\n🎉 HSDT implementation is VALIDATED! 🎉\n');
else
    fprintf('\n⚠️  HSDT implementation needs review ⚠️\n');
end

fprintf('\nTest completed.\n');