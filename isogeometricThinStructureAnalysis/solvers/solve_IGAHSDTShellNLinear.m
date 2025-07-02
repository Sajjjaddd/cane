function [dHat, CPHistory, resHistory, isConverged, BSplinePatch, minElSize, ...
          centerDisplacementHistory, loadHistory] = ...
    solve_IGAHSDTShellNLinear...
    (BSplinePatch, propNLinearAnalysis, solve_LinearSystem, ...
    plot_IGANonlinear, graph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
% 
% Returns the displacement field, convergence history, and other results for 
% nonlinear HSDT shell analysis with 5 DOF per control point using
% Newton-Raphson iteration with load stepping
%
%                Input :
%         BSplinePatch : Structure containing HSDT shell patch information
%  propNLinearAnalysis : Structure containing nonlinear analysis properties:
%                        .method : 'newtonRaphson'
%                        .noLoadSteps : Number of load steps
%                        .eps : Convergence tolerance
%                        .maxIter : Maximum iterations per load step
%                        .conservativeLoad : Load conservativeness flag
%   solve_LinearSystem : Function handle for linear system solver
%    plot_IGANonlinear : Plotting function handle (optional)
%                graph : Graphics properties structure
%               outMsg : Output message flag
%
%               Output :
%                 dHat : Final displacement field [5*noCPs x 1]
%            CPHistory : Control point displacement history
%           resHistory : Residual norm history
%          isConverged : Convergence flag for each load step
%         BSplinePatch : Updated patch structure
%             minElSize : Minimum element size
% centerDisplacementHistory : Center point w-displacement history
%          loadHistory : Load factor history for each converged step
%
% Function layout:
%
% 0. Read input and initialize
%
% 1. Setup nonlinear analysis parameters
%
% 2. Load stepping loop
% ->
%    2i. Initialize load step
%
%   2ii. Newton-Raphson iteration loop
%   ->
%        2ii.1. Compute tangent stiffness and residual
%
%        2ii.2. Apply boundary conditions
%
%        2ii.3. Solve linear system
%
%        2ii.4. Update displacement field
%
%        2ii.5. Check convergence
%
%        2ii.6. Update iteration counter
%   <-
%   2iii. Store load step results
% <-
%
% 3. Finalize and output results
%
%% Function main body

if strcmp(outMsg, 'outputEnabled')
    fprintf('_________________________________________________________\n');
    fprintf('#########################################################\n');
    fprintf('Nonlinear analysis for an isogeometric HSDT shell\n');
    fprintf('with 5 DOF per control point has been initiated\n');
    fprintf('Method: %s\n', propNLinearAnalysis.method);
    fprintf('Load steps: %d\n', propNLinearAnalysis.noLoadSteps);
    fprintf('_________________________________________________________\n\n');
    tic;
end

%% 0. Read input and initialize

% Define analysis type
analysis.type = 'isogeometricHSDTShellAnalysis';

% Get patch properties
CP = BSplinePatch.CP;
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));
numDOFs = 5 * numCPs_xi * numCPs_eta;

% Create DOF numbering (5 DOF per control point)
BSplinePatch.DOFNumbering = zeros(numCPs_xi, numCPs_eta, 5);
k = 1;
for cpj = 1:numCPs_eta
    for cpi = 1:numCPs_xi
        BSplinePatch.DOFNumbering(cpi, cpj, 1) = k;     % u
        BSplinePatch.DOFNumbering(cpi, cpj, 2) = k + 1; % v
        BSplinePatch.DOFNumbering(cpi, cpj, 3) = k + 2; % w
        BSplinePatch.DOFNumbering(cpi, cpj, 4) = k + 3; % θx
        BSplinePatch.DOFNumbering(cpi, cpj, 5) = k + 4; % θy
        k = k + 5;
    end
end

% Boundary conditions
homDOFs = BSplinePatch.homDOFs;
freeDOFs = 1:numDOFs;
freeDOFs(ismember(freeDOFs, homDOFs)) = [];

% Initialize solution vectors
dHat = zeros(numDOFs, 1);
dHatIncrement = zeros(numDOFs, 1);

%% 1. Setup nonlinear analysis parameters

noLoadSteps = propNLinearAnalysis.noLoadSteps;
tolerance = propNLinearAnalysis.eps;
maxIterations = propNLinearAnalysis.maxIter;

% Load stepping
loadFactors = linspace(0, 1, noLoadSteps + 1);
loadFactors = loadFactors(2:end);  % Remove zero load factor

% History storage
CPHistory = zeros(numDOFs, noLoadSteps);
resHistory = zeros(maxIterations, noLoadSteps);
isConverged = false(noLoadSteps, 1);

% Center point tracking for load-displacement curve
centerCP_xi = ceil(numCPs_xi / 2);
centerCP_eta = ceil(numCPs_eta / 2);
centerDOF_w = BSplinePatch.DOFNumbering(centerCP_xi, centerCP_eta, 3); % w-displacement
centerDisplacementHistory = zeros(noLoadSteps, 1);
loadHistory = zeros(noLoadSteps, 1);

% Dummy variables for compatibility
dHatSaved = 'undefined';
dHatDot = 'undefined';
dHatDotSaved = 'undefined';
connections = 'undefined';
propCoupling = 'undefined';
propStrDynamics = 'undefined';
isReferenceUpdated = false;
noWeakDBCCnd = 0;
t = 0;
tab = '  ';

%% 2. Load stepping loop

for iLoadStep = 1:noLoadSteps
    
    currentLoadFactor = loadFactors(iLoadStep);
    
    if strcmp(outMsg, 'outputEnabled')
        fprintf('\n--- Load Step %d/%d (Load Factor: %.4f) ---\n', ...
                iLoadStep, noLoadSteps, currentLoadFactor);
    end
    
    %% 2i. Initialize load step
    
    % Reset iteration counter
    iNLinearIter = 0;
    isLoadStepConverged = false;
    residualNorm = inf;
    
    %% 2ii. Newton-Raphson iteration loop
    
    while ~isLoadStepConverged && iNLinearIter < maxIterations
        
        iNLinearIter = iNLinearIter + 1;
        
        if strcmp(outMsg, 'outputEnabled')
            fprintf('    Iteration %d: ', iNLinearIter);
        end
        
        %% 2ii.1. Compute tangent stiffness and residual
        
        [KTangent, residualVct, ~, ~, minElSize] = ...
            computeTangentStiffMtxResVctIGAHSDTShellNLinear...
            ([], [], dHat, dHatSaved, dHatDot, dHatDotSaved, ...
            BSplinePatch, connections, propCoupling, currentLoadFactor, ...
            1, iLoadStep, iNLinearIter, noWeakDBCCnd, t, propStrDynamics, ...
            isReferenceUpdated, tab, '');
        
        %% 2ii.2. Apply boundary conditions
        
        % Apply homogeneous boundary conditions to system
        KTangentFree = KTangent(freeDOFs, freeDOFs);
        residualFree = residualVct(freeDOFs);
        
        %% 2ii.3. Solve linear system
        
        if rcond(KTangentFree) < 1e-12
            warning('Tangent stiffness matrix is ill-conditioned (rcond = %.2e)', rcond(KTangentFree));
        end
        
        deltaU_free = solve_LinearSystem(KTangentFree, -residualFree);
        
        %% 2ii.4. Update displacement field
        
        % Initialize increment
        dHatIncrement = zeros(numDOFs, 1);
        dHatIncrement(freeDOFs) = deltaU_free;
        
        % Update total displacement
        dHat = dHat + dHatIncrement;
        
        %% 2ii.5. Check convergence
        
        residualNorm = norm(residualFree);
        displacementNorm = norm(deltaU_free);
        
        % Store residual history
        resHistory(iNLinearIter, iLoadStep) = residualNorm;
        
        % Convergence criteria
        if residualNorm < tolerance
            isLoadStepConverged = true;
            if strcmp(outMsg, 'outputEnabled')
                fprintf('CONVERGED (Residual: %.3e)\n', residualNorm);
            end
        else
            if strcmp(outMsg, 'outputEnabled')
                fprintf('Residual: %.3e, Disp increment: %.3e\n', ...
                        residualNorm, displacementNorm);
            end
        end
        
        %% 2ii.6. Update iteration counter and check divergence
        
        if residualNorm > 1e6
            warning('Analysis appears to be diverging (residual = %.2e)', residualNorm);
            break;
        end
        
    end  % End Newton-Raphson iteration loop
    
    %% 2iii. Store load step results
    
    isConverged(iLoadStep) = isLoadStepConverged;
    CPHistory(:, iLoadStep) = dHat;
    
    % Store center displacement and current load for load-displacement curve
    centerDisplacementHistory(iLoadStep) = dHat(centerDOF_w);
    loadHistory(iLoadStep) = currentLoadFactor;
    
    if ~isLoadStepConverged
        warning('Load step %d did not converge after %d iterations', ...
                iLoadStep, maxIterations);
    end
    
    % Optional plotting during iterations
    if ~strcmp(plot_IGANonlinear, 'undefined') && strcmp(outMsg, 'outputEnabled')
        % Call plotting function if provided
        % plot_IGANonlinear(BSplinePatch, dHat, graph, iLoadStep);
    end
    
    % Print load step summary
    if strcmp(outMsg, 'outputEnabled')
        fprintf('    Load Step %d: Load Factor = %.4f, Center w = %.6e\n', ...
                iLoadStep, currentLoadFactor, centerDisplacementHistory(iLoadStep));
    end
    
end  % End load stepping loop

%% 3. Finalize and output results

% Extract final displacement components for analysis
numCPs = numCPs_xi * numCPs_eta;
u_final = dHat(1:5:end);      % u-displacements
v_final = dHat(2:5:end);      % v-displacements  
w_final = dHat(3:5:end);      % w-displacements
theta_x_final = dHat(4:5:end); % θx-rotations
theta_y_final = dHat(5:5:end); % θy-rotations

if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('\n=== Nonlinear HSDT Analysis Summary ===\n');
    fprintf('Total computation time: %.2f seconds\n', computationalTime);
    fprintf('Converged load steps: %d/%d\n', sum(isConverged), noLoadSteps);
    
    fprintf('\nFinal displacement ranges:\n');
    fprintf('u: [%.6e, %.6e]\n', min(u_final), max(u_final));
    fprintf('v: [%.6e, %.6e]\n', min(v_final), max(v_final));
    fprintf('w: [%.6e, %.6e]\n', min(w_final), max(w_final));
    fprintf('θx: [%.6e, %.6e] rad\n', min(theta_x_final), max(theta_x_final));
    fprintf('θy: [%.6e, %.6e] rad\n', min(theta_y_final), max(theta_y_final));
    
    fprintf('\nMaximum values:\n');
    fprintf('|w|_max: %.6e\n', max(abs(w_final)));
    fprintf('|θx|_max: %.6e rad (%.3f deg)\n', max(abs(theta_x_final)), rad2deg(max(abs(theta_x_final))));
    fprintf('|θy|_max: %.6e rad (%.3f deg)\n', max(abs(theta_y_final)), rad2deg(max(abs(theta_y_final))));
    
    fprintf('\n_____________Nonlinear HSDT Analysis Ended______________\n');
    fprintf('#######################################################\n\n');
    
    % Plot load-displacement curve for center point
    if strcmp(outMsg, 'outputEnabled')
        plotLoadDisplacementCurve(centerDisplacementHistory, loadHistory, ...
                                  BSplinePatch, graph);
    end
end

end


function plotLoadDisplacementCurve(centerDisplacementHistory, loadHistory, ...
                                   BSplinePatch, graph)
%% Plot load-displacement curve for center point of shell
%
% This function plots the load factor vs center displacement curve
% showing the nonlinear behavior during load stepping

% Get load amplitude for scaling
if isfield(BSplinePatch, 'FGamma')
    totalLoad = norm(BSplinePatch.FGamma);
else
    totalLoad = 1.0; % Default scaling
end

% Convert to physical units
if isfield(BSplinePatch, 'NBC') && ~isempty(BSplinePatch.NBC.loadAmplitude)
    loadAmplitude = BSplinePatch.NBC.loadAmplitude{1};
else
    loadAmplitude = -90; % Default for Scordelis-Lo roof
end

% Create new figure
figure_handle = figure('Name', 'Load-Displacement Curve (Center Point)', ...
                      'Position', [100, 100, 800, 600]);

% Plot load-displacement curve
actualLoads = loadHistory * loadAmplitude;
displacements_mm = centerDisplacementHistory * 1000; % Convert to mm

plot(abs(displacements_mm), abs(actualLoads), 'o-', ...
     'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'blue', ...
     'Color', 'blue');

% Formatting
grid on;
xlabel('Center Point |w| Displacement [mm]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Applied Load Magnitude [N/m²]', 'FontSize', 12, 'FontWeight', 'bold');
title('Nonlinear Load-Displacement Curve (HSDT Shell Center)', ...
      'FontSize', 14, 'FontWeight', 'bold');

% Add annotations
text(0.7*max(abs(displacements_mm)), 0.3*max(abs(actualLoads)), ...
     sprintf('Final displacement: %.3f mm\nFinal load: %.1f N/m²', ...
             abs(displacements_mm(end)), abs(actualLoads(end))), ...
     'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Enhance plot appearance
set(gca, 'FontSize', 11);
set(gca, 'LineWidth', 1.2);
box on;

% Add gradient to show loading direction
hold on;
for i = 1:length(displacements_mm)-1
    plot(abs(displacements_mm(i:i+1)), abs(actualLoads(i:i+1)), ...
         'Color', [0, 0, 1-0.1*i/length(displacements_mm)], 'LineWidth', 2);
end
hold off;

% Display curve characteristics
fprintf('\n=== Load-Displacement Curve Analysis ===\n');
fprintf('Center point coordinates: CP(%d,%d)\n', ...
        ceil(size(BSplinePatch.CP,1)/2), ceil(size(BSplinePatch.CP,2)/2));
fprintf('Maximum displacement: %.6f mm\n', max(abs(displacements_mm)));
fprintf('Maximum load: %.1f N/m²\n', max(abs(actualLoads)));

% Check for nonlinear behavior
if length(displacements_mm) > 1
    % Linear stiffness from first two points
    linearStiffness = (actualLoads(2) - actualLoads(1)) / ...
                     (displacements_mm(2) - displacements_mm(1));
    
    % Final stiffness from last two points  
    finalStiffness = (actualLoads(end) - actualLoads(end-1)) / ...
                    (displacements_mm(end) - displacements_mm(end-1));
    
    stiffnessRatio = finalStiffness / linearStiffness;
    
    fprintf('Initial stiffness: %.2e N/mm/m²\n', linearStiffness);
    fprintf('Final stiffness: %.2e N/mm/m²\n', finalStiffness);
    fprintf('Stiffness ratio (final/initial): %.3f\n', stiffnessRatio);
    
    if stiffnessRatio < 0.8
        fprintf('*** Significant stiffness softening detected ***\n');
    elseif stiffnessRatio > 1.2
        fprintf('*** Significant stiffness hardening detected ***\n');
    else
        fprintf('Moderate nonlinear stiffness change\n');
    end
end

fprintf('==========================================\n\n');

end