function plotLoadDisplacementHistory(centerDisplacementHistory, loadHistory, ...
                                     isConverged, BSplinePatch, propNLinearAnalysis, ...
                                     compareLinear, linearCenterDisplacement)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Plots comprehensive load-displacement curves for center point of shell
% during nonlinear HSDT analysis, showing convergence behavior and
% comparison with linear analysis if available.
%
%                Input :
% centerDisplacementHistory : Center point w-displacement for each load step
%          loadHistory : Load factor history for each load step  
%          isConverged : Convergence flag for each load step
%         BSplinePatch : Patch structure containing load information
%  propNLinearAnalysis : Nonlinear analysis properties
%        compareLinear : Flag to compare with linear results (optional)
% linearCenterDisplacement : Linear center displacement for comparison (optional)
%
%               Output :
%                     : Creates comprehensive load-displacement plots
%
% Function layout :
%
% 1. Process input data and extract load information
%
% 2. Create main load-displacement curve plot
%
% 3. Add convergence indicators and annotations
%
% 4. Optional: Add linear comparison
%
% 5. Create secondary analysis plots
%
%% Function main body

%% 1. Process input data

% Set default values for optional inputs
if nargin < 6
    compareLinear = false;
end
if nargin < 7
    linearCenterDisplacement = 0;
end

% Extract load amplitude
if isfield(BSplinePatch, 'NBC') && ~isempty(BSplinePatch.NBC.loadAmplitude)
    loadAmplitude = BSplinePatch.NBC.loadAmplitude{1};
else
    loadAmplitude = -90; % Default for Scordelis-Lo roof
end

% Convert to physical units
actualLoads = loadHistory * loadAmplitude;
displacements_mm = centerDisplacementHistory * 1000; % Convert to mm

% Separate converged and non-converged steps
convergedLoads = actualLoads(isConverged);
convergedDisplacements = displacements_mm(isConverged);
nonConvergedLoads = actualLoads(~isConverged);
nonConvergedDisplacements = displacements_mm(~isConverged);

%% 2. Create main load-displacement curve plot

figure('Name', 'Comprehensive Load-Displacement Analysis', ...
       'Position', [50, 50, 1200, 800]);

% Main plot
subplot(2, 2, [1 2]);
hold on;

% Plot converged steps
if ~isempty(convergedLoads)
    plot(abs(convergedDisplacements), abs(convergedLoads), 'o-', ...
         'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'green', ...
         'Color', 'blue', 'DisplayName', 'Converged Steps');
end

% Plot non-converged steps (if any)
if ~isempty(nonConvergedLoads)
    plot(abs(nonConvergedDisplacements), abs(nonConvergedLoads), 'x', ...
         'LineWidth', 2, 'MarkerSize', 10, 'Color', 'red', ...
         'DisplayName', 'Non-Converged Steps');
end

% Add linear comparison if requested
if compareLinear && linearCenterDisplacement ~= 0
    linearDisp_mm = abs(linearCenterDisplacement) * 1000;
    linearLoad = abs(loadAmplitude);
    plot([0, linearDisp_mm], [0, linearLoad], '--', ...
         'LineWidth', 2, 'Color', 'black', 'DisplayName', 'Linear Response');
    
    % Mark linear solution point
    plot(linearDisp_mm, linearLoad, 'd', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'yellow', ...
         'MarkerEdgeColor', 'black', 'DisplayName', 'Linear Solution');
end

% Formatting
grid on;
xlabel('Center Point |w| Displacement [mm]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Applied Load Magnitude [N/m²]', 'FontSize', 12, 'FontWeight', 'bold');
title('Nonlinear Load-Displacement Curve (HSDT Shell Center)', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast');

% Add load step annotations
for i = 1:length(actualLoads)
    if isConverged(i)
        text(abs(displacements_mm(i)), abs(actualLoads(i)), sprintf('%d', i), ...
             'FontSize', 8, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'Color', 'blue');
    end
end

hold off;

%% 3. Convergence analysis subplot

subplot(2, 2, 3);
bar(1:length(isConverged), isConverged, 'FaceColor', [0.2 0.8 0.2]);
xlabel('Load Step');
ylabel('Converged (1=Yes, 0=No)');
title('Convergence History');
grid on;
ylim([-0.1, 1.1]);

% Add convergence percentage
convergenceRate = sum(isConverged) / length(isConverged) * 100;
text(0.7*length(isConverged), 0.8, sprintf('%.1f%% Converged', convergenceRate), ...
     'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'white');

%% 4. Stiffness analysis subplot

subplot(2, 2, 4);
if length(convergedDisplacements) > 1
    % Calculate incremental stiffness
    stiffness = zeros(length(convergedDisplacements)-1, 1);
    for i = 1:length(stiffness)
        deltaLoad = convergedLoads(i+1) - convergedLoads(i);
        deltaDisp = convergedDisplacements(i+1) - convergedDisplacements(i);
        if abs(deltaDisp) > 1e-12
            stiffness(i) = deltaLoad / deltaDisp;
        end
    end
    
    plot(2:length(convergedDisplacements), stiffness, 's-', ...
         'LineWidth', 2, 'MarkerSize', 6, 'Color', 'purple');
    xlabel('Load Step');
    ylabel('Incremental Stiffness [N/mm/m²]');
    title('Stiffness Evolution');
    grid on;
    
    % Add trend line
    if length(stiffness) > 2
        hold on;
        p = polyfit(2:length(convergedDisplacements), stiffness', 1);
        trendLine = polyval(p, 2:length(convergedDisplacements));
        plot(2:length(convergedDisplacements), trendLine, '--', ...
             'LineWidth', 1.5, 'Color', 'red');
        
        if p(1) < 0
            text(0.7*length(convergedDisplacements), 0.7*max(stiffness), ...
                 'Stiffness Softening', 'FontSize', 9, 'Color', 'red', ...
                 'FontWeight', 'bold');
        else
            text(0.7*length(convergedDisplacements), 0.7*max(stiffness), ...
                 'Stiffness Hardening', 'FontSize', 9, 'Color', 'blue', ...
                 'FontWeight', 'bold');
        end
        hold off;
    end
else
    text(0.5, 0.5, 'Insufficient data for stiffness analysis', ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
    xlim([0 1]);
    ylim([0 1]);
end

%% 5. Display comprehensive analysis results

fprintf('\n=== Comprehensive Load-Displacement Analysis ===\n');
fprintf('Total load steps: %d\n', length(isConverged));
fprintf('Converged steps: %d (%.1f%%)\n', sum(isConverged), convergenceRate);
fprintf('Non-converged steps: %d\n', sum(~isConverged));

if ~isempty(convergedDisplacements)
    fprintf('\nFinal converged state:\n');
    fprintf('  Displacement: %.6f mm\n', abs(convergedDisplacements(end)));
    fprintf('  Load: %.1f N/m²\n', abs(convergedLoads(end)));
    
    % Compare with linear if available
    if compareLinear && linearCenterDisplacement ~= 0
        linearDisp_mm = abs(linearCenterDisplacement) * 1000;
        nonlinearityFactor = abs(convergedDisplacements(end)) / linearDisp_mm;
        fprintf('\nLinear vs Nonlinear comparison:\n');
        fprintf('  Linear displacement: %.6f mm\n', linearDisp_mm);
        fprintf('  Nonlinear displacement: %.6f mm\n', abs(convergedDisplacements(end)));
        fprintf('  Nonlinearity factor: %.3f\n', nonlinearityFactor);
        
        if nonlinearityFactor > 1.1
            fprintf('  *** Significant stiffness softening (%.1f%% increase) ***\n', ...
                    (nonlinearityFactor-1)*100);
        elseif nonlinearityFactor < 0.9
            fprintf('  *** Significant stiffness hardening (%.1f%% decrease) ***\n', ...
                    (1-nonlinearityFactor)*100);
        else
            fprintf('  Small nonlinear effects (%.1f%% change)\n', ...
                    abs(nonlinearityFactor-1)*100);
        end
    end
end

fprintf('===============================================\n\n');

% Enhance overall plot appearance
sgtitle(sprintf('Nonlinear HSDT Analysis: %d Load Steps, %.1f%% Convergence', ...
                length(isConverged), convergenceRate), ...
        'FontSize', 16, 'FontWeight', 'bold');

end