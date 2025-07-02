function [dg1, dg2] = computeCurrentBaseVectorsHSDT...
    (i, p, j, q, CP, dHatEl, nDrvBaseVct, dR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns the covariant base vectors of the current (deformed) configuration
% for HSDT shell analysis with 5 DOF per control point
%
%                Input :
%                  i,p : Knot span and polynomial degree in xi-direction
%                  j,q : Knot span and polynomial degree in eta-direction
%                   CP : Control Point coordinates and weights
%               dHatEl : Element displacement vector [5*(p+1)*(q+1) x 1]
%                       DOF ordering: [u, v, w, θx, θy] per control point
%           nDrvBaseVct : Number of derivatives to compute for base vectors
%                   dR : Basis functions and derivatives
%
%               Output :
%                 dg1 : Covariant base vector g1 and its derivatives in 
%                       current configuration [3 x (nDrvBaseVct+1)]
%                 dg2 : Covariant base vector g2 and its derivatives in
%                       current configuration [3 x (nDrvBaseVct+1)]
%
% Function layout :
%
% 1. Compute reference base vectors
%
% 2. Extract displacement components from element vector
%
% 3. Compute current base vectors = reference + displacement gradients
%
%% Function main body

%% 1. Compute reference base vectors

% Compute reference base vectors using existing function
[dG1, dG2] = computeBaseVectorsAndDerivativesForBSplineSurface...
    (i, p, j, q, CP, nDrvBaseVct, dR);

%% 2. Extract displacement components from element vector

% Number of control points in element
numCPsEl = (p + 1) * (q + 1);

% Initialize displacement arrays
u_el = zeros(numCPsEl, 1);  % u-displacements
v_el = zeros(numCPsEl, 1);  % v-displacements  
w_el = zeros(numCPsEl, 1);  % w-displacements

% Extract displacements from element vector (5 DOF per control point)
for iCP = 1:numCPsEl
    u_el(iCP) = dHatEl(5*(iCP-1) + 1);  % u-displacement
    v_el(iCP) = dHatEl(5*(iCP-1) + 2);  % v-displacement
    w_el(iCP) = dHatEl(5*(iCP-1) + 3);  % w-displacement
    % dHatEl(5*(iCP-1) + 4) = θx-rotation (not used for base vectors)
    % dHatEl(5*(iCP-1) + 5) = θy-rotation (not used for base vectors)
end

%% 3. Compute current base vectors = reference + displacement gradients

% Initialize current base vectors
dg1 = zeros(3, nDrvBaseVct + 1);
dg2 = zeros(3, nDrvBaseVct + 1);

% Loop over derivative orders
for iDrv = 1:nDrvBaseVct + 1
    
    % Initialize displacement gradients
    du_dxi = 0; du_deta = 0;
    dv_dxi = 0; dv_deta = 0;
    dw_dxi = 0; dw_deta = 0;
    
    % Determine derivative indices for basis functions
    if iDrv == 1
        % 0th derivative (base vectors)
        dR_xi_idx = 2;   % ∂R/∂ξ
        dR_eta_idx = 4;  % ∂R/∂η
    elseif iDrv == 2 && nDrvBaseVct >= 1
        % 1st derivative (base vector derivatives)
        dR_xi_idx = 3;   % ∂²R/∂ξ²
        dR_eta_idx = 5;  % ∂²R/∂ξ∂η
    else
        continue;  % Higher derivatives not implemented
    end
    
    % Compute displacement gradients
    for iCP = 1:numCPsEl
        du_dxi = du_dxi + dR(iCP, dR_xi_idx) * u_el(iCP);
        du_deta = du_deta + dR(iCP, dR_eta_idx) * u_el(iCP);
        
        dv_dxi = dv_dxi + dR(iCP, dR_xi_idx) * v_el(iCP);
        dv_deta = dv_deta + dR(iCP, dR_eta_idx) * v_el(iCP);
        
        dw_dxi = dw_dxi + dR(iCP, dR_xi_idx) * w_el(iCP);
        dw_deta = dw_deta + dR(iCP, dR_eta_idx) * w_el(iCP);
    end
    
    % Current base vectors = Reference base vectors + displacement gradients
    % g1 = G1 + [∂u/∂ξ; ∂v/∂ξ; ∂w/∂ξ]
    % g2 = G2 + [∂u/∂η; ∂v/∂η; ∂w/∂η]
    dg1(:, iDrv) = dG1(:, iDrv) + [du_dxi; dv_dxi; dw_dxi];
    dg2(:, iDrv) = dG2(:, iDrv) + [du_deta; dv_deta; dw_deta];
end

end