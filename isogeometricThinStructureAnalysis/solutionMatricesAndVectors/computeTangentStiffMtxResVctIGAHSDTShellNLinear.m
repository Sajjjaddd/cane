function [tangentStiffMtx, residualVct, BSplinePatch, propCoupling, minElArea] = ...
    computeTangentStiffMtxResVctIGAHSDTShellNLinear...
    (KConstant, tanMtxLoad, dHat, dHatSaved, dHatDot, dHatDotSaved, ...
    BSplinePatch, connections, propCoupling, loadFactor, noPatch, ...
    noTimeStep, iNLinearIter, noWeakDBCCnd, t, propStrDynamics, ...
    isReferenceUpdated, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns the tangent stiffness matrix and the residual vector for the 
% nonlinear HSDT shell formulation with 5 DOF per control point and 
% geometric nonlinearity.
%
% HSDT DOF ordering per control point: [u, v, w, θx, θy]
%
%                 Input : 
%             KConstant : Dummy variable for this function
%            tanMtxLoad : Dummy variable for this function
%                  dHat : The displacement field of the current iteration
%             dHatSaved : The displacement field of the previous time step
%               dHatDot : The velocity field (dummy for static analysis)
%          dHatDotSaved : The velocity field of previous time step (dummy)
%          BSplinePatch : The B-Spline patch structure with HSDT information
%           connections : Dummy variable for this function
%          propCoupling : Dummy variable for this function
%            loadFactor : Load factor for continuation methods
%               noPatch : Dummy variable for this function
%            noTimeStep : Current time step number
%          iNLinearIter : Current nonlinear iteration number
%          noWeakDBCCnd : Dummy variable for this function
%                     t : Time instance
%       propStrDynamics : Dummy variable for this function
%    isReferenceUpdated : Flag for reference configuration update
%                   tab : Tabulation for output formatting
%                outMsg : Output message flag
%
%                Output :
%       tangentStiffMtx : The tangent stiffness matrix
%           residualVct : The residual force vector
%          BSplinePatch : Updated patch structure
%          propCoupling : Dummy output
%             minElArea : Minimum element area in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Compute material matrices
%
% 2. Initialize output arrays
%
% 3. Loop over elements
% ->
%    3i. Create Element Freedom Table
%
%   3ii. Initialize element matrices
%
%  3iii. Loop over Gauss Points
%  ->
%        3iii.1. Compute basis functions and derivatives
%
%        3iii.2. Compute current geometry (deformed configuration)
%
%        3iii.3. Compute current displacement gradients and rotations
%
%        3iii.4. Compute strains in current configuration
%
%        3iii.5. Compute stresses from strains
%
%        3iii.6. Compute tangent stiffness matrix components
%
%        3iii.7. Compute residual force components
%
%        3iii.8. Add contributions to element matrices
%  <-
%   3iv. Assemble element contributions to global matrices
% <-
%
% 4. Apply load factor to external forces
%
% 5. Compute residual vector
%
%% Function main body

%% 0. Read input

% Check input
isBSplinePatchCell = false;
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input');
    else
        BSplinePatch = BSplinePatch{1};
    end
    isBSplinePatchCell = true;
end

% Extract patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
parameters = BSplinePatch.parameters;
int = BSplinePatch.int;
DOFNumbering = BSplinePatch.DOFNumbering;
NBC = BSplinePatch.NBC;

% Geometric properties
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Initialize minimum element area
tolerance = 1e-4;
if abs(CP(1, 1, 1) - CP(numCPs_xi, 1, 1)) >= tolerance
    minElArea = abs(CP(1, 1, 1) - CP(numCPs_xi, 1, 1));
else
    minElArea = CP(1, 1, 1) - CP(1, numCPs_eta, 1);
end

% Total number of DOFs (5 per control point)
numDOFs = 5 * numCPs_xi * numCPs_eta;

% Element DOFs
numDOFsEl = 5 * (p + 1) * (q + 1);

%% 1. Compute material matrices

% Membrane material matrix
Dm = parameters.E * parameters.t / (1 - parameters.nue^2) * ...
    [1              parameters.nue 0
     parameters.nue 1              0
     0              0              (1 - parameters.nue)/2];

% Bending material matrix
Db = parameters.E * parameters.t^3 / (12 * (1 - parameters.nue^2)) * ...
    [1              parameters.nue 0
     parameters.nue 1              0
     0              0              (1 - parameters.nue)/2];

% Shear material matrix
G = parameters.E / (2 * (1 + parameters.nue));
if isfield(parameters, 'shearCorrection')
    shearCorrection = parameters.shearCorrection;
else
    shearCorrection = 5/6;
end
Ds = shearCorrection * G * parameters.t * [1 0; 0 1];

%% 2. Initialize output arrays
tangentStiffMtx = zeros(numDOFs, numDOFs);
residualVct = zeros(numDOFs, 1);

% Integration rule
if strcmp(int.type, 'default')
    numGP_xi = p + 1;
    numGP_eta = q + 1;
elseif strcmp(int.type, 'user')
    numGP_xi = int.xiNGP;
    numGP_eta = int.etaNGP;
end

[xiGP, xiGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_xi);
[etaGP, etaGW] = getGaussPointsAndWeightsOverUnitDomain(numGP_eta);

%% 3. Loop over elements
for j = q + 1:numKnots_eta - q - 1
    for i = p + 1:numKnots_xi - p - 1
        if Xi(i + 1) ~= Xi(i) && Eta(j + 1) ~= Eta(j)
            
            %% 3i. Create Element Freedom Table
            EFT = zeros(1, numDOFsEl);
            k = 1;
            for cpj = j - q:j
                for cpi = i - p:i
                    EFT(k) = DOFNumbering(cpi, cpj, 1);     % u
                    EFT(k + 1) = DOFNumbering(cpi, cpj, 2); % v
                    EFT(k + 2) = DOFNumbering(cpi, cpj, 3); % w
                    EFT(k + 3) = DOFNumbering(cpi, cpj, 4); % θx
                    EFT(k + 4) = DOFNumbering(cpi, cpj, 5); % θy
                    k = k + 5;
                end
            end
            
            % Extract element displacement vector
            dHatEl = dHat(EFT);
            
            %% 3ii. Initialize element matrices
            KeTotal = zeros(numDOFsEl, numDOFsEl);
            FeResidual = zeros(numDOFsEl, 1);
            
            % Jacobian for parametric transformation
            detJxiu = (Xi(i + 1) - Xi(i)) * (Eta(j + 1) - Eta(j)) / 4;
            
            elementArea = 0;
            
            %% 3iii. Loop over Gauss Points
            for iGP_eta = 1:numGP_eta
                for iGP_xi = 1:numGP_xi
                    
                    %% 3iii.1. Compute basis functions and derivatives
                    xi = (Xi(i + 1) + Xi(i) + xiGP(iGP_xi) * (Xi(i + 1) - Xi(i))) / 2;
                    eta = (Eta(j + 1) + Eta(j) + etaGP(iGP_eta) * (Eta(j + 1) - Eta(j))) / 2;
                    
                    nDrvBasis = 2;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (i, p, xi, Xi, j, q, eta, Eta, CP, isNURBS, nDrvBasis);
                    
                    %% 3iii.2. Compute current geometry (deformed configuration)
                    nDrvBaseVct = 1;
                    [dG1, dG2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (i, p, j, q, CP, nDrvBaseVct, dR);
                    
                    % Current configuration includes displacements
                    [dg1, dg2] = computeCurrentBaseVectorsHSDT...
                        (i, p, j, q, CP, dHatEl, nDrvBaseVct, dR);
                    
                    % Reference and current surface normals
                    G3Tilde = cross(dG1(:, 1), dG2(:, 1));
                    g3Tilde = cross(dg1(:, 1), dg2(:, 1));
                    
                    dA_ref = norm(G3Tilde);
                    dA_cur = norm(g3Tilde);
                    
                    G3 = G3Tilde / dA_ref;
                    g3 = g3Tilde / dA_cur;
                    
                    %% 3iii.3. Compute current displacement gradients and rotations
                    [strainsCurrent, rotationsCurrent] = computeCurrentStrainsRotationsHSDT...
                        (dR, dg1, dg2, g3, dHatEl, p, q);
                    
                    %% 3iii.4. Compute strains in current configuration
                    % Membrane strains (Green-Lagrange)
                    strainsMembrane = strainsCurrent.membrane;
                    
                    % Bending strains (curvature changes)
                    strainsBending = strainsCurrent.bending;
                    
                    % Shear strains
                    strainsShear = strainsCurrent.shear;
                    
                    %% 3iii.5. Compute stresses from strains
                    stressesMembrane = Dm * strainsMembrane;
                    stressesBending = Db * strainsBending;
                    stressesShear = Ds * strainsShear;
                    
                    %% 3iii.6. Compute tangent stiffness matrix components
                    [KeTangent] = computeElementTangentStiffnessHSDT...
                        (p, q, dR, dg1, dg2, g3, dHatEl, Dm, Db, Ds, ...
                        stressesMembrane, stressesBending, stressesShear);
                    
                    %% 3iii.7. Compute residual force components
                    [FeInternal] = computeElementInternalForcesHSDT...
                        (p, q, dR, dg1, dg2, g3, stressesMembrane, stressesBending, stressesShear);
                    
                    %% 3iii.8. Add contributions to element matrices
                    weightGP = dA_cur * detJxiu * xiGW(iGP_xi) * etaGW(iGP_eta);
                    
                    KeTotal = KeTotal + KeTangent * weightGP;
                    FeResidual = FeResidual + FeInternal * weightGP;
                    
                    elementArea = elementArea + weightGP;
                end
            end
            
            %% 3iv. Assemble element contributions to global matrices
            tangentStiffMtx(EFT, EFT) = tangentStiffMtx(EFT, EFT) + KeTotal;
            residualVct(EFT) = residualVct(EFT) + FeResidual;
            
            % Update minimum element area
            if elementArea < minElArea
                minElArea = elementArea;
            end
        end
    end
end

%% 4. Apply load factor to external forces
if isfield(BSplinePatch, 'FGamma')
    FExternal = loadFactor * BSplinePatch.FGamma;
else
    FExternal = zeros(numDOFs, 1);
end

%% 5. Compute residual vector
residualVct = residualVct - FExternal;

% Output formatting
if strcmp(outMsg, 'outputEnabled')
    fprintf('%s Tangent stiffness matrix and residual vector computation completed\n', tab);
    fprintf('%s Nonlinear iteration: %d, Load factor: %.4f\n', tab, iNLinearIter, loadFactor);
    fprintf('%s Residual norm: %.6e\n', tab, norm(residualVct));
end

% Handle cell output
if isBSplinePatchCell
    BSplinePatch = {BSplinePatch};
end

end