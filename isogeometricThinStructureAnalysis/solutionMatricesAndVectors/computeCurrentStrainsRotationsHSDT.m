function [strainsCurrent, rotationsCurrent] = computeCurrentStrainsRotationsHSDT...
    (dR, dg1, dg2, g3, dHatEl, p, q)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns the current strains and rotations for HSDT nonlinear analysis
% with 5 DOF per control point in the deformed configuration
%
%                Input :
%                   dR : Basis functions and derivatives
%                 dg1 : Current covariant base vector g1 and derivatives
%                 dg2 : Current covariant base vector g2 and derivatives  
%                  g3 : Current surface normal vector (unit vector)
%               dHatEl : Element displacement vector [5*(p+1)*(q+1) x 1]
%                       DOF ordering: [u, v, w, θx, θy] per control point
%                  p,q : Polynomial degrees
%
%               Output :
%       strainsCurrent : Structure containing current strains
%                        .membrane [3x1] - membrane strains [εξξ; εηη; γξη]
%                        .bending [3x1] - bending strains [κξξ; κηη; κξη] 
%                        .shear [2x1] - shear strains [γξz; γηz]
%      rotationsCurrent : Structure containing current rotations
%                         .theta_x - θx rotation values at Gauss point
%                         .theta_y - θy rotation values at Gauss point
%
% Function layout :
%
% 1. Extract rotation components from element vector
%
% 2. Compute current rotations at Gauss point
%
% 3. Compute membrane strains (Green-Lagrange)
%
% 4. Compute bending strains (curvature changes)
%
% 5. Compute shear strains (transverse shear)
%
%% Function main body

%% 1. Extract rotation components from element vector

% Number of control points in element
numCPsEl = (p + 1) * (q + 1);

% Initialize rotation arrays
theta_x_el = zeros(numCPsEl, 1);  % θx-rotations
theta_y_el = zeros(numCPsEl, 1);  % θy-rotations

% Extract rotations from element vector (5 DOF per control point)
for iCP = 1:numCPsEl
    theta_x_el(iCP) = dHatEl(5*(iCP-1) + 4);  % θx-rotation
    theta_y_el(iCP) = dHatEl(5*(iCP-1) + 5);  % θy-rotation
end

%% 2. Compute current rotations at Gauss point

% Interpolate rotations at Gauss point
theta_x_GP = 0;
theta_y_GP = 0;

% Derivatives of rotations at Gauss point
dtheta_x_dxi = 0;  dtheta_x_deta = 0;
dtheta_y_dxi = 0;  dtheta_y_deta = 0;

for iCP = 1:numCPsEl
    % Rotation values
    theta_x_GP = theta_x_GP + dR(iCP, 1) * theta_x_el(iCP);
    theta_y_GP = theta_y_GP + dR(iCP, 1) * theta_y_el(iCP);
    
    % Rotation derivatives
    dtheta_x_dxi = dtheta_x_dxi + dR(iCP, 2) * theta_x_el(iCP);   % ∂θx/∂ξ
    dtheta_x_deta = dtheta_x_deta + dR(iCP, 4) * theta_x_el(iCP); % ∂θx/∂η
    
    dtheta_y_dxi = dtheta_y_dxi + dR(iCP, 2) * theta_y_el(iCP);   % ∂θy/∂ξ  
    dtheta_y_deta = dtheta_y_deta + dR(iCP, 4) * theta_y_el(iCP); % ∂θy/∂η
end

%% 3. Compute membrane strains (Green-Lagrange)

% Current covariant base vectors
g1 = dg1(:, 1);
g2 = dg2(:, 1);

% Current covariant metric tensor components
g11 = g1' * g1;
g22 = g2' * g2;
g12 = g1' * g2;

% Reference metric tensor (unit metric in parametric space)
% For NURBS surfaces, this needs to be computed from reference base vectors
% For simplicity, assuming orthogonal parametrization: G11=G22=1, G12=0
G11 = 1.0;  % Could be computed as G1'*G1 if reference base vectors available
G22 = 1.0;  % Could be computed as G2'*G2 if reference base vectors available  
G12 = 0.0;  % Could be computed as G1'*G2 if reference base vectors available

% Green-Lagrange strains in parametric coordinates
% E_ξξ = (g11 - G11)/2, E_ηη = (g22 - G22)/2, E_ξη = (g12 - G12)/2
strain_mem_xi_xi = 0.5 * (g11 - G11);
strain_mem_eta_eta = 0.5 * (g22 - G22);
strain_mem_xi_eta = 0.5 * (g12 - G12);

% Engineering shear strain (factor of 2)
gamma_mem_xi_eta = 2.0 * strain_mem_xi_eta;

strainsMembrane = [strain_mem_xi_xi; strain_mem_eta_eta; gamma_mem_xi_eta];

%% 4. Compute bending strains (curvature changes)

% HSDT bending strains are directly related to rotation gradients
% κξξ = -∂θx/∂ξ, κηη = -∂θy/∂η, κξη = -(∂θx/∂η + ∂θy/∂ξ)
strain_bend_xi_xi = -dtheta_x_dxi;
strain_bend_eta_eta = -dtheta_y_deta;
strain_bend_xi_eta = -(dtheta_x_deta + dtheta_y_dxi);

strainsBending = [strain_bend_xi_xi; strain_bend_eta_eta; strain_bend_xi_eta];

%% 5. Compute shear strains (transverse shear)

% Compute displacement gradients for w (out-of-plane displacement)
dw_dxi = 0;  dw_deta = 0;

for iCP = 1:numCPsEl
    w_el = dHatEl(5*(iCP-1) + 3);  % w-displacement
    dw_dxi = dw_dxi + dR(iCP, 2) * w_el;   % ∂w/∂ξ
    dw_deta = dw_deta + dR(iCP, 4) * w_el; % ∂w/∂η
end

% HSDT shear strains: γξz = ∂w/∂ξ + θx, γηz = ∂w/∂η + θy
strain_shear_xi_z = dw_dxi + theta_x_GP;
strain_shear_eta_z = dw_deta + theta_y_GP;

strainsShear = [strain_shear_xi_z; strain_shear_eta_z];

%% Assemble output structures

strainsCurrent.membrane = strainsMembrane;
strainsCurrent.bending = strainsBending;
strainsCurrent.shear = strainsShear;

rotationsCurrent.theta_x = theta_x_GP;
rotationsCurrent.theta_y = theta_y_GP;
rotationsCurrent.dtheta_x_dxi = dtheta_x_dxi;
rotationsCurrent.dtheta_x_deta = dtheta_x_deta;
rotationsCurrent.dtheta_y_dxi = dtheta_y_dxi;
rotationsCurrent.dtheta_y_deta = dtheta_y_deta;

end