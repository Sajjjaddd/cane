function [KeTangent] = computeElementTangentStiffnessHSDT...
    (p, q, dR, dg1, dg2, g3, dHatEl, Dm, Db, Ds, ...
    stressesMembrane, stressesBending, stressesShear)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns the element tangent stiffness matrix for HSDT nonlinear analysis
% including material and geometric stiffness contributions
%
%                Input :
%                  p,q : Polynomial degrees
%                   dR : Basis functions and derivatives
%            dg1, dg2 : Current covariant base vectors and derivatives
%                  g3 : Current surface normal vector (unit vector)
%               dHatEl : Element displacement vector [5*(p+1)*(q+1) x 1]
%            Dm,Db,Ds : Material matrices (membrane, bending, shear)
%      stressesMembrane : Current membrane stresses [3x1]
%      stressesBending : Current bending stresses [3x1]
%       stressesShear : Current shear stresses [2x1]
%
%               Output :
%           KeTangent : Element tangent stiffness matrix
%                       [5*(p+1)*(q+1) x 5*(p+1)*(q+1)]
%
% Function layout :
%
% 1. Initialize element matrix
%
% 2. Compute material stiffness matrix (linear part)
%
% 3. Compute geometric stiffness matrix (nonlinear part)
%
% 4. Assemble total tangent stiffness matrix
%
%% Function main body

%% 1. Initialize element matrix

% Number of control points and DOFs in element
numCPsEl = (p + 1) * (q + 1);
numDOFsEl = 5 * numCPsEl;

% Initialize tangent stiffness matrix
KeTangent = zeros(numDOFsEl, numDOFsEl);

%% 2. Compute material stiffness matrix (linear part)

% This is similar to the linear HSDT stiffness but computed in current configuration
KeMaterial = computeMaterialStiffnessHSDT(p, q, dR, dg1, dg2, g3, Dm, Db, Ds);

%% 3. Compute geometric stiffness matrix (nonlinear part)

% Geometric stiffness arises from stress-displacement coupling
KeGeometric = computeGeometricStiffnessHSDT...
    (p, q, dR, dg1, dg2, g3, stressesMembrane, stressesBending, stressesShear);

%% 4. Assemble total tangent stiffness matrix

KeTangent = KeMaterial + KeGeometric;

end


function [KeMaterial] = computeMaterialStiffnessHSDT(p, q, dR, dg1, dg2, g3, Dm, Db, Ds)
%% Sub-function: Material stiffness matrix in current configuration

% Number of DOFs
numCPsEl = (p + 1) * (q + 1);
numDOFsEl = 5 * numCPsEl;

% Initialize B-operator matrices
BMembrane = zeros(3, numDOFsEl);
BBending = zeros(3, numDOFsEl);
BShear = zeros(2, numDOFsEl);

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

% Compute material stiffness contributions
KeMembrane = BMembraneLC' * Dm * BMembraneLC;
KeBending = BBendingLC' * Db * BBendingLC;
KeShear = BShearLC' * Ds * BShearLC;

% Total material stiffness
KeMaterial = KeMembrane + KeBending + KeShear;

end


function [KeGeometric] = computeGeometricStiffnessHSDT...
    (p, q, dR, dg1, dg2, g3, stressesMembrane, stressesBending, stressesShear)
%% Sub-function: Geometric stiffness matrix

% Number of DOFs
numCPsEl = (p + 1) * (q + 1);
numDOFsEl = 5 * numCPsEl;

% Initialize geometric stiffness matrix
KeGeometric = zeros(numDOFsEl, numDOFsEl);

% Current base vectors  
g1 = dg1(:, 1);
g2 = dg2(:, 1);

% For HSDT, geometric stiffness mainly comes from membrane stress effects
% This is a simplified implementation - full geometric stiffness requires
% careful treatment of stress-displacement coupling

% Extract stress components
Nxi = stressesMembrane(1);    % Membrane stress in ξ-direction
Neta = stressesMembrane(2);   % Membrane stress in η-direction
Nxieta = stressesMembrane(3); % Membrane shear stress

% Loop over DOF pairs for geometric stiffness
for iDOF = 1:numDOFsEl
    for jDOF = 1:numDOFsEl
        
        % Get local node numbers and directions
        i_node = ceil(iDOF/5);
        i_dir = iDOF - 5*(i_node - 1);
        j_node = ceil(jDOF/5);
        j_dir = jDOF - 5*(j_node - 1);
        
        % Geometric stiffness mainly affects displacement DOFs (1,2,3)
        if i_dir <= 3 && j_dir <= 3
            
            % Simplified geometric stiffness contribution
            % Full implementation would require detailed stress-gradient coupling
            Kg_contrib = 0;
            
            % Contributions from membrane stresses to geometric stiffness
            if i_dir == j_dir
                Kg_contrib = Kg_contrib + Nxi * dR(i_node, 2) * dR(j_node, 2); % ξ-direction
                Kg_contrib = Kg_contrib + Neta * dR(i_node, 4) * dR(j_node, 4); % η-direction
                Kg_contrib = Kg_contrib + Nxieta * (dR(i_node, 2) * dR(j_node, 4) + ...
                                                    dR(i_node, 4) * dR(j_node, 2)); % mixed
            end
            
            KeGeometric(iDOF, jDOF) = Kg_contrib;
        end
    end
end

end