function [G3, dA, GabCovariant, GabContravariant, GContravariant, eLC, BV] = ...
    computeG3AndCotraBaseVectorsCovConMetricsAndCurvatureTensor ...
    (GCovariant, dGCovariant)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the covariant, contravariant and local Cartesian base vectors,
% the second derivatives of the position vector i.e. the first derivatives
% of the covariant base vectors as well as the curvature coefficients for a
% 3D NURBS surface.
%
%             Input :
%        gCovariant : The covariant base vectors
%       dGCovariant : The derivatives of the covariant base vectorsghts
%
%            Output :
%                G3 : The normal to the NURBS surface base vector
%                dA : The area spanned by the covariant base vectors
%      GabCovariant : The covariant metric coefficients
%  GabContravariant : The contravariant metric coefficients
%    GContravariant : The contravariant base vectors
%               eLC : The local Cartesian basis
%                BV : The curvature coefficients in a Voigt notation
%
% Functiona layout :
%
% 1. Compute the normal vector to the surface vector G3Tilde which is not normalized
%
% 2. Compute the area pf the element spanned by g1 and g2
%
% 3. Compute the normalized orthogonal to the surface base vector
%
% 4. Compute the covariant metric coefficient tensor
%
% 5. Compute the contravariant metric coefficient tensor
%
% 6. Compute the contravariant base vectors
%
% 7. Compute the local cartesian base vectors
%
% 8. Compute the curvature vector
%
%% Function main body

%% 1. Compute the normal vector to the surface vector G3Tilde which is not normalized
G3Tilde = cross(GCovariant(:,1),GCovariant(:,2));

%% 2. Compute the area pf the element spanned by g1 and g2
dA = norm(G3Tilde);

%% 3. Compute the normalized orthogonal to the surface base vector
G3 = G3Tilde/dA;

%% 4. Compute the covariant metric coefficient tensor
GabCovariant = GCovariant'*GCovariant; 

%% 5. Compute the contravariant metric coefficient tensor
GabContravariant = inv(GabCovariant);

%% 6. Compute the contravariant base vectors
GContravariant = GCovariant*GabContravariant';

%% 7. Compute the local cartesian base vectors
eLC(:,1) = GCovariant(:,1)/sqrt(GCovariant(:,1)'*GCovariant(:,1));
eLC(:,2) = GContravariant(:,2)/sqrt(GContravariant(:,2)'*GContravariant(:,2));
eLC(:,3) = G3;

%% 8. Compute the curvature vector
BV = dGCovariant'*G3;

% curvature tensor
% bab = [bv(1) bv(3); bv(3) bv(2)];

end
