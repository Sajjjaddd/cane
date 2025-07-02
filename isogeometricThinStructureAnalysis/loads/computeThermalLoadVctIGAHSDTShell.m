function [FThermal] = computeThermalLoadVctIGAHSDTShell...
    (FThermal, BSplinePatch, xiLoadExtension, etaLoadExtension, ...
    temperatureData, thermalLoadType, isConservative, t, intDomain, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns the thermal load vector for HSDT shell analysis with 5 DOF per
% control point, supporting both constant and variable through-thickness
% temperature distributions.
%
%                Input :
%             FThermal : Existing thermal load vector to add to
%         BSplinePatch : The B-Spline patch structure for HSDT
%     xiLoadExtension : Load extension in xi-direction [xi_start, xi_end]
%    etaLoadExtension : Load extension in eta-direction [eta_start, eta_end]
%      temperatureData : Temperature data structure:
%                        .T0 : Reference temperature [K]
%                        .T : Temperature field function or constant [K]
%                        .alpha : Thermal expansion coefficient [1/K]
%                        .throughThickness : 'constant' or 'variable'
%                        .zCoords : z-coordinates for variable temperature (optional)
%                        .TProfile : Temperature profile function of z (optional)
%     thermalLoadType : Type of thermal load ('uniform', 'gradient', 'custom')
%       isConservative : Flag for conservative thermal loads
%                    t : Time parameter (for time-dependent temperature)
%            intDomain : Integration domain structure
%               outMsg : Output message flag
%
%               Output :
%             FThermal : Updated thermal load vector [5*noCPs x 1]
%
% Function layout :
%
% 1. Read input and initialize
%
% 2. Setup integration parameters
%
% 3. Loop over thermal load domain
% ->
%    3i. Compute temperature field at current point
%    3ii. Calculate thermal strains
%    3iii. Compute thermal force contributions
%    3iv. Assemble to global thermal load vector
% <-
%
% 4. Apply thermal expansion effects
%
%% Function main body

if strcmp(outMsg, 'outputEnabled')
    fprintf('_______________________________________________________\n');
    fprintf('#######################################################\n');
    fprintf('Computation of thermal load vector for HSDT shell\n');
    fprintf('analysis has been initiated\n');
    fprintf('_______________________________________________________\n\n');
end

%% 1. Read input and initialize

% Get patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
parameters = BSplinePatch.parameters;

% Material properties
E = parameters.E;
nue = parameters.nue;
thickness = parameters.t;

% Thermal properties
T0 = temperatureData.T0;          % Reference temperature
alpha = temperatureData.alpha;    % Thermal expansion coefficient

% Determine temperature field type
if isfield(temperatureData, 'throughThickness')
    throughThicknessType = temperatureData.throughThickness;
else
    throughThicknessType = 'constant';
end

% Get current temperature
if isa(temperatureData.T, 'function_handle')
    getCurrentTemp = temperatureData.T;
    isTemperatureFunction = true;
else
    currentTemp = temperatureData.T;
    isTemperatureFunction = false;
end

%% 2. Setup integration parameters

% Integration domain bounds
xiStart = xiLoadExtension(1);
xiEnd = xiLoadExtension(2);
etaStart = etaLoadExtension(1);
etaEnd = etaLoadExtension(2);

% Element-wise integration
xiSpan = findspan(Xi, xiStart, p);
etaSpan = findspan(Eta, etaStart, q);
xiSpanEnd = findspan(Xi, xiEnd, p);
etaSpanEnd = findspan(Eta, etaEnd, q);

% Number of Gauss points
if strcmp(intDomain.type, 'default')
    numGP_xi = p + 1;
    numGP_eta = q + 1;
else
    numGP_xi = intDomain.xiNGP;
    numGP_eta = intDomain.etaNGP;
end

[xiGP, xiGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_xi);
[etaGP, etaGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_eta);

%% 3. Loop over thermal load domain

% Loop over elements in the load domain
for j = etaSpan:etaSpanEnd
    for i = xiSpan:xiSpanEnd
        
        % Check if element exists
        if Xi(i+1) ~= Xi(i) && Eta(j+1) ~= Eta(j)
            
            % Element bounds
            xiLeft = max(Xi(i), xiStart);
            xiRight = min(Xi(i+1), xiEnd);
            etaLeft = max(Eta(j), etaStart);
            etaRight = min(Eta(j+1), etaEnd);
            
            % Skip if no overlap
            if xiLeft >= xiRight || etaLeft >= etaRight
                continue;
            end
            
            % Element freedom table for HSDT (5 DOF per CP)
            EFT = zeros(1, 5*(p+1)*(q+1));
            k = 1;
            for cpj = j-q:j
                for cpi = i-p:i
                    EFT(k) = BSplinePatch.DOFNumbering(cpi, cpj, 1);     % u
                    EFT(k+1) = BSplinePatch.DOFNumbering(cpi, cpj, 2);   % v
                    EFT(k+2) = BSplinePatch.DOFNumbering(cpi, cpj, 3);   % w
                    EFT(k+3) = BSplinePatch.DOFNumbering(cpi, cpj, 4);   % θx
                    EFT(k+4) = BSplinePatch.DOFNumbering(cpi, cpj, 5);   % θy
                    k = k + 5;
                end
            end
            
            % Element thermal load vector
            FThermalEl = zeros(length(EFT), 1);
            
            % Integration domain transformation
            detJxiu = (xiRight - xiLeft) * (etaRight - etaLeft) / 4;
            
            %% 3i-3iv. Gauss point integration
            for iGP_eta = 1:numGP_eta
                for iGP_xi = 1:numGP_xi
                    
                    % Map to physical domain
                    xi = (xiLeft + xiRight + xiGP(iGP_xi) * (xiRight - xiLeft)) / 2;
                    eta = (etaLeft + etaRight + etaGP(iGP_eta) * (etaRight - etaLeft)) / 2;
                    
                    % Check if point is within load domain
                    if xi < xiStart || xi > xiEnd || eta < etaStart || eta > etaEnd
                        continue;
                    end
                    
                    % Compute basis functions and derivatives
                    nDrvBasis = 1;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (i, p, xi, Xi, j, q, eta, Eta, CP, isNURBS, nDrvBasis);
                    
                    % Compute base vectors
                    nDrvBaseVct = 1;
                    [dG1, dG2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (i, p, j, q, CP, nDrvBaseVct, dR);
                    
                    % Surface normal and area element
                    G3Tilde = cross(dG1(:,1), dG2(:,1));
                    dA = norm(G3Tilde);
                    
                    %% 3i. Compute temperature field at current point
                    if isTemperatureFunction
                        % Get physical coordinates
                        x = zeros(3, 1);
                        for iCP = 1:(p+1)*(q+1)
                            x = x + dR(iCP, 1) * CP(i-p+floor((iCP-1)/(p+1)), ...
                                                   j-q+mod(iCP-1, q+1)+1, 1:3)';
                        end
                        T_current = getCurrentTemp(x(1), x(2), x(3), t);
                    else
                        T_current = currentTemp;
                    end
                    
                    %% 3ii. Calculate thermal strains
                    [thermalStrains] = computeThermalStrains...
                        (T_current, T0, alpha, thickness, throughThicknessType, temperatureData);
                    
                    %% 3iii. Compute thermal force contributions
                    [FThermalGP] = computeThermalForceContributions...
                        (thermalStrains, dR, dG1, dG2, parameters, p, q);
                    
                    %% 3iv. Add to element thermal load vector
                    weightGP = dA * detJxiu * xiGW(iGP_xi) * etaGW(iGP_eta);
                    FThermalEl = FThermalEl + FThermalGP * weightGP;
                    
                end
            end
            
            % Assemble element contribution to global thermal load vector
            FThermal(EFT) = FThermal(EFT) + FThermalEl;
            
        end
    end
end

%% 4. Final processing

if strcmp(outMsg, 'outputEnabled')
    thermalLoadMagnitude = norm(FThermal);
    fprintf('Thermal load vector computation completed\n');
    fprintf('Temperature: %.1f K (ΔT = %.1f K)\n', ...
            temperatureData.T, temperatureData.T - T0);
    fprintf('Thermal expansion coefficient: %.2e /K\n', alpha);
    fprintf('Through-thickness type: %s\n', throughThicknessType);
    fprintf('Thermal load magnitude: %.6e\n', thermalLoadMagnitude);
    fprintf('_______________________________________________________\n\n');
end

end


function [thermalStrains] = computeThermalStrains...
    (T_current, T0, alpha, thickness, throughThicknessType, temperatureData)
%% Compute thermal strains for HSDT

% Temperature difference
deltaT = T_current - T0;

switch throughThicknessType
    case 'constant'
        % Constant temperature through thickness
        % Only membrane thermal strains
        epsilon_thermal_membrane = alpha * deltaT * [1; 1; 0]; % [εxx; εyy; γxy]
        epsilon_thermal_bending = [0; 0; 0];  % No bending thermal strains
        
    case 'variable'
        % Variable temperature through thickness
        if isfield(temperatureData, 'TProfile') && isa(temperatureData.TProfile, 'function_handle')
            % Custom temperature profile T(z) where z ∈ [-t/2, t/2]
            TProfile = temperatureData.TProfile;
            
            % Integrate membrane thermal strain (average through thickness)
            zCoords = linspace(-thickness/2, thickness/2, 21);
            T_profile_vals = zeros(size(zCoords));
            for iz = 1:length(zCoords)
                T_profile_vals(iz) = TProfile(zCoords(iz));
            end
            
            % Membrane strain (average)
            T_avg = trapz(zCoords, T_profile_vals) / thickness;
            deltaT_avg = T_avg - T0;
            epsilon_thermal_membrane = alpha * deltaT_avg * [1; 1; 0];
            
            % Bending strain (moment of temperature)
            T_moment = trapz(zCoords, T_profile_vals .* zCoords);
            T_moment_normalized = T_moment / (thickness^3/12);
            epsilon_thermal_bending = alpha * T_moment_normalized * [1; 1; 0];
            
        else
            % Linear temperature variation (simplified)
            % T(z) = T_mid + (T_top - T_bottom) * z / thickness
            if isfield(temperatureData, 'T_top') && isfield(temperatureData, 'T_bottom')
                T_top = temperatureData.T_top;
                T_bottom = temperatureData.T_bottom;
                T_mid = (T_top + T_bottom) / 2;
                
                % Membrane thermal strain
                deltaT_avg = T_mid - T0;
                epsilon_thermal_membrane = alpha * deltaT_avg * [1; 1; 0];
                
                % Bending thermal strain
                deltaT_gradient = (T_top - T_bottom) / thickness;
                epsilon_thermal_bending = alpha * deltaT_gradient * thickness/2 * [1; 1; 0];
            else
                % Fallback to constant
                epsilon_thermal_membrane = alpha * deltaT * [1; 1; 0];
                epsilon_thermal_bending = [0; 0; 0];
            end
        end
        
    otherwise
        error('Unknown through-thickness temperature type: %s', throughThicknessType);
end

% Package thermal strains
thermalStrains.membrane = epsilon_thermal_membrane;
thermalStrains.bending = epsilon_thermal_bending;
thermalStrains.shear = [0; 0]; % No thermal shear strains for HSDT

end


function [FThermalGP] = computeThermalForceContributions...
    (thermalStrains, dR, dG1, dG2, parameters, p, q)
%% Compute thermal force contributions at Gauss point

% Material matrices
E = parameters.E;
nue = parameters.nue;
t = parameters.t;

% Membrane material matrix
Dm = E * t / (1 - nue^2) * ...
    [1    nue  0
     nue  1    0
     0    0    (1-nue)/2];

% Bending material matrix  
Db = E * t^3 / (12 * (1 - nue^2)) * ...
    [1    nue  0
     nue  1    0
     0    0    (1-nue)/2];

% Number of DOFs per element
numDOFsEl = 5 * (p + 1) * (q + 1);
FThermalGP = zeros(numDOFsEl, 1);

% Current base vectors
g1 = dG1(:, 1);
g2 = dG2(:, 1);

% Compute contravariant base vectors
gContra = [g1 g2]' * [g1 g2];
gContraInv = inv(gContra);
g1_contra = gContraInv(1,1) * g1 + gContraInv(1,2) * g2;
g2_contra = gContraInv(2,1) * g1 + gContraInv(2,2) * g2;

% Local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface([g1 g2], [g1_contra g2_contra]);

% Transformation matrix
TContra2LC = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell(eLC, [g1_contra g2_contra]);

% Thermal stresses
thermalStresses_membrane = Dm * thermalStrains.membrane;
thermalStresses_bending = Db * thermalStrains.bending;

% Compute B-operator matrices
[BMembrane, BBending] = computeBOperatorMatricesHSDT(dR, g1, g2, p, q);

% Transform to local Cartesian system
BMembraneLC = TContra2LC * BMembrane;
BBendingLC = TContra2LC * BBending;

% Thermal force contributions (negative of internal forces due to thermal stresses)
FThermalGP = FThermalGP - BMembraneLC' * thermalStresses_membrane;
FThermalGP = FThermalGP - BBendingLC' * thermalStresses_bending;

end


function [BMembrane, BBending] = computeBOperatorMatricesHSDT(dR, g1, g2, p, q)
%% Compute B-operator matrices for membrane and bending

numCPsEl = (p + 1) * (q + 1);
numDOFsEl = 5 * numCPsEl;

BMembrane = zeros(3, numDOFsEl);
BBending = zeros(3, numDOFsEl);

% Loop over all DOFs
for iDOF = 1:numDOFsEl
    % Get control point and DOF type
    k = ceil(iDOF/5);
    dir = iDOF - 5*(k - 1);
    
    switch dir
        case 1 % u-displacement
            BMembrane(1, iDOF) = dR(k, 2);      % ∂u/∂ξ contribution to εξξ
            BMembrane(3, iDOF) = 0.5 * dR(k, 4); % 0.5*∂u/∂η contribution to γξη
            
        case 2 % v-displacement
            BMembrane(2, iDOF) = dR(k, 4);      % ∂v/∂η contribution to εηη
            BMembrane(3, iDOF) = 0.5 * dR(k, 2); % 0.5*∂v/∂ξ contribution to γξη
            
        case 3 % w-displacement
            % No contribution to membrane or bending strains
            
        case 4 % θx-rotation
            BBending(1, iDOF) = -dR(k, 2);      % -∂θx/∂ξ contribution to κξξ
            BBending(3, iDOF) = -0.5 * dR(k, 4); % -0.5*∂θx/∂η contribution to κξη
            
        case 5 % θy-rotation
            BBending(2, iDOF) = -dR(k, 4);      % -∂θy/∂η contribution to κηη
            BBending(3, iDOF) = -0.5 * dR(k, 2); % -0.5*∂θy/∂ξ contribution to κξη
    end
end

end