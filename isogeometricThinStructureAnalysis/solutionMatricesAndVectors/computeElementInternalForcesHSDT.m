function [FeInternal] = computeElementInternalForcesHSDT...
    (p, q, dR, dg1, dg2, g3, stressesMembrane, stressesBending, stressesShear)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns the element internal force vector for HSDT nonlinear analysis
% corresponding to the current stress state
%
%                Input :
%                  p,q : Polynomial degrees
%                   dR : Basis functions and derivatives
%            dg1, dg2 : Current covariant base vectors and derivatives
%                  g3 : Current surface normal vector (unit vector)
%      stressesMembrane : Current membrane stresses [3x1] - [Nξξ; Nηη; Nξη]
%      stressesBending : Current bending stresses [3x1] - [Mξξ; Mηη; Mξη]
%       stressesShear : Current shear stresses [2x1] - [Qξ; Qη]
%
%               Output :
%           FeInternal : Element internal force vector [5*(p+1)*(q+1) x 1]
%                       corresponding to current stress state
%
% Function layout :
%
% 1. Initialize element force vector
%
% 2. Compute B-operator matrices in current configuration
%
% 3. Compute internal forces from virtual work principle
%
% 4. Assemble total internal force vector
%
%% Function main body

%% 1. Initialize element force vector

% Number of control points and DOFs in element
numCPsEl = (p + 1) * (q + 1);
numDOFsEl = 5 * numCPsEl;

% Initialize internal force vector
FeInternal = zeros(numDOFsEl, 1);

%% 2. Compute B-operator matrices in current configuration

% Current base vectors
g1 = dg1(:, 1);
g2 = dg2(:, 1);

% Compute contravariant base vectors
gContra = [g1 g2]' * [g1 g2];
gContraInv = inv(gContra);
g1_contra = gContraInv(1,1) * g1 + gContraInv(1,2) * g2;
g2_contra = gContraInv(2,1) * g1 + gContraInv(2,2) * g2;

% Local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface([g1 g2], [g1_contra g2_contra]);

% Transformation matrix
TContra2LC = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell(eLC, [g1_contra g2_contra]);

% Initialize B-operator matrices
BMembrane = zeros(3, numDOFsEl);
BBending = zeros(3, numDOFsEl);
BShear = zeros(2, numDOFsEl);

% Loop over all DOFs to compute B-operator matrices
for iDOF = 1:numDOFsEl
    % Compute local node number and DOF direction
    k = ceil(iDOF/5);
    dir = iDOF - 5*(k - 1);
    
    % Initialize B-operator entries
    BMembrane(:, iDOF) = 0;
    BBending(:, iDOF) = 0; 
    BShear(:, iDOF) = 0;
    
    % DOF contributions based on type
    switch dir
        case 1 % u-displacement
            BMembrane(1, iDOF) = dR(k, 2) * g1(1);  % ∂u/∂ξ contribution to εξξ
            BMembrane(3, iDOF) = 0.5 * dR(k, 4) * g2(1); % 0.5*∂u/∂η contribution to γξη
            
        case 2 % v-displacement
            BMembrane(2, iDOF) = dR(k, 4) * g2(2);  % ∂v/∂η contribution to εηη
            BMembrane(3, iDOF) = 0.5 * dR(k, 2) * g1(2); % 0.5*∂v/∂ξ contribution to γξη
            
        case 3 % w-displacement
            BShear(1, iDOF) = dR(k, 2);  % ∂w/∂ξ contribution to γξz
            BShear(2, iDOF) = dR(k, 4);  % ∂w/∂η contribution to γηz
            
        case 4 % θx-rotation
            BBending(1, iDOF) = -dR(k, 2);  % -∂θx/∂ξ contribution to κξξ
            BBending(3, iDOF) = -0.5 * dR(k, 4); % -0.5*∂θx/∂η contribution to κξη
            BShear(1, iDOF) = dR(k, 1);  % θx contribution to γξz
            
        case 5 % θy-rotation
            BBending(2, iDOF) = -dR(k, 4);  % -∂θy/∂η contribution to κηη  
            BBending(3, iDOF) = -0.5 * dR(k, 2); % -0.5*∂θy/∂ξ contribution to κξη
            BShear(2, iDOF) = dR(k, 1);  % θy contribution to γηz
    end
end

% Transform B-operator matrices to local Cartesian system
BMembraneLC = TContra2LC * BMembrane;
BBendingLC = TContra2LC * BBending;
BShearLC = TContra2LC(1:2, 1:2) * BShear;

%% 3. Compute internal forces from virtual work principle

% Internal forces from virtual work: F_int = B^T * stress
FeMembrane = BMembraneLC' * stressesMembrane;
FeBending = BBendingLC' * stressesBending;
FeShear = BShearLC' * stressesShear;

%% 4. Assemble total internal force vector

FeInternal = FeMembrane + FeBending + FeShear;

end