function [homDOFs, supportInfo] = create3DSupports_HSDT...
    (homDOFs, supportConfigs, CP, DOFNumbering, supportType, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Creates comprehensive 3D supports for HSDT shell analysis with full 
% 5 DOF support including rotations. Supports various support types
% including fixed, pinned, roller, and custom constraint combinations.
%
%                Input :
%               homDOFs : Existing homogeneous DOF array to append to
%        supportConfigs : Array of support configuration structures:
%                         .location : Support location specification
%                           .type : 'point', 'edge', 'corner', 'region'
%                           .xi : xi coordinates [xi_start, xi_end] or xi_point
%                           .eta : eta coordinates [eta_start, eta_end] or eta_point
%                         .constraints : DOF constraint specification
%                           .u : true/false for u-displacement constraint
%                           .v : true/false for v-displacement constraint  
%                           .w : true/false for w-displacement constraint
%                           .thetax : true/false for θx-rotation constraint
%                           .thetay : true/false for θy-rotation constraint
%                         .values : Constraint values (for inhomogeneous BC)
%                         .stiffness : Support stiffness (for elastic supports)
%                    CP : Control point coordinates and weights
%          DOFNumbering : DOF numbering array [noCPs_xi x noCPs_eta x 5]
%           supportType : Global support type override ('fixed', 'pinned', etc.)
%                outMsg : Output message flag
%
%               Output :
%               homDOFs : Updated homogeneous DOF array
%           supportInfo : Detailed support information structure
%
% Function layout :
%
% 1. Input validation and initialization
%
% 2. Process each support configuration
% ->
%    2i. Determine support location and affected control points
%    2ii. Apply constraint specifications for all 5 DOF
%    2iii. Handle special support types (fixed, pinned, etc.)
%    2iv. Store support information
% <-
%
% 3. Generate comprehensive support summary
%
%% Function main body

if strcmp(outMsg, 'outputEnabled')
    fprintf('_________________________________________________________\n');
    fprintf('#########################################################\n');
    fprintf('Creating 3D supports for HSDT shell analysis\n');
    fprintf('Full 5 DOF support: [u, v, w, θx, θy]\n');
    fprintf('_________________________________________________________\n\n');
end

%% 1. Input validation and initialization

% Initialize support information
supportInfo = struct();
supportInfo.totalSupports = 0;
supportInfo.constrainedDOFs = [0, 0, 0, 0, 0]; % Count for each DOF type
supportInfo.supports = [];

% Get mesh dimensions
numCPs_xi = size(CP, 1);
numCPs_eta = size(CP, 2);
totalCPs = numCPs_xi * numCPs_eta;

% Validate DOF numbering
if size(DOFNumbering, 3) ~= 5
    error('DOF numbering must have 5 components for HSDT analysis');
end

% Initialize constraints tracking
constraintMatrix = false(totalCPs, 5); % Track which DOFs are constrained

%% 2. Process each support configuration

for iSupport = 1:length(supportConfigs)
    config = supportConfigs(iSupport);
    
    if strcmp(outMsg, 'outputEnabled')
        fprintf('Processing support %d/%d...\n', iSupport, length(supportConfigs));
    end
    
    %% 2i. Determine support location and affected control points
    
    [affectedCPs] = findAffectedControlPoints(config.location, numCPs_xi, numCPs_eta);
    
    %% 2ii. Apply constraint specifications for all 5 DOF
    
    % Determine constraints based on configuration or support type
    if isfield(config, 'constraints')
        constraints = config.constraints;
    else
        constraints = getDefaultConstraints(supportType);
    end
    
    % Override with global support type if specified
    if nargin >= 5 && ~isempty(supportType)
        constraints = getDefaultConstraints(supportType);
    end
    
    %% 2iii. Apply constraints to affected control points
    
    newDOFs = [];
    constraintCount = [0, 0, 0, 0, 0];
    
    for iCP = 1:length(affectedCPs)
        cpIndex = affectedCPs(iCP);
        [cpi, cpj] = ind2sub([numCPs_xi, numCPs_eta], cpIndex);
        
        % Apply u-displacement constraint
        if isfield(constraints, 'u') && constraints.u
            newDOF = DOFNumbering(cpi, cpj, 1);
            if ~ismember(newDOF, homDOFs)
                newDOFs = [newDOFs, newDOF];
                constraintMatrix(cpIndex, 1) = true;
                constraintCount(1) = constraintCount(1) + 1;
            end
        end
        
        % Apply v-displacement constraint
        if isfield(constraints, 'v') && constraints.v
            newDOF = DOFNumbering(cpi, cpj, 2);
            if ~ismember(newDOF, homDOFs)
                newDOFs = [newDOFs, newDOF];
                constraintMatrix(cpIndex, 2) = true;
                constraintCount(2) = constraintCount(2) + 1;
            end
        end
        
        % Apply w-displacement constraint
        if isfield(constraints, 'w') && constraints.w
            newDOF = DOFNumbering(cpi, cpj, 3);
            if ~ismember(newDOF, homDOFs)
                newDOFs = [newDOFs, newDOF];
                constraintMatrix(cpIndex, 3) = true;
                constraintCount(3) = constraintCount(3) + 1;
            end
        end
        
        % Apply θx-rotation constraint
        if isfield(constraints, 'thetax') && constraints.thetax
            newDOF = DOFNumbering(cpi, cpj, 4);
            if ~ismember(newDOF, homDOFs)
                newDOFs = [newDOFs, newDOF];
                constraintMatrix(cpIndex, 4) = true;
                constraintCount(4) = constraintCount(4) + 1;
            end
        end
        
        % Apply θy-rotation constraint
        if isfield(constraints, 'thetay') && constraints.thetay
            newDOF = DOFNumbering(cpi, cpj, 5);
            if ~ismember(newDOF, homDOFs)
                newDOFs = [newDOFs, newDOF];
                constraintMatrix(cpIndex, 5) = true;
                constraintCount(5) = constraintCount(5) + 1;
            end
        end
    end
    
    % Add new DOFs to homogeneous constraint list
    homDOFs = [homDOFs, newDOFs];
    
    %% 2iv. Store support information
    
    support = struct();
    support.id = iSupport;
    support.location = config.location;
    support.constraints = constraints;
    support.affectedCPs = affectedCPs;
    support.constrainedDOFs = newDOFs;
    support.constraintCount = constraintCount;
    
    if isfield(config, 'description')
        support.description = config.description;
    else
        support.description = sprintf('Support %d', iSupport);
    end
    
    supportInfo.supports = [supportInfo.supports; support];
    supportInfo.totalSupports = supportInfo.totalSupports + 1;
    supportInfo.constrainedDOFs = supportInfo.constrainedDOFs + constraintCount;
    
    if strcmp(outMsg, 'outputEnabled')
        fprintf('  %s: %d CPs, [%d %d %d %d %d] DOFs constrained\n', ...
                support.description, length(affectedCPs), constraintCount);
    end
end

%% 3. Generate comprehensive support summary

% Sort homogeneous DOFs
homDOFs = sort(unique(homDOFs));

% Calculate support statistics
supportInfo.totalConstrainedDOFs = length(homDOFs);
supportInfo.constraintMatrix = constraintMatrix;
supportInfo.totalFreeDOFs = 5 * totalCPs - supportInfo.totalConstrainedDOFs;

if strcmp(outMsg, 'outputEnabled')
    fprintf('\n=== 3D Support Summary ===\n');
    fprintf('Total supports created: %d\n', supportInfo.totalSupports);
    fprintf('Total control points: %d\n', totalCPs);
    fprintf('Total DOFs: %d (5 per control point)\n', 5 * totalCPs);
    fprintf('Constrained DOFs: %d\n', supportInfo.totalConstrainedDOFs);
    fprintf('Free DOFs: %d\n', supportInfo.totalFreeDOFs);
    fprintf('\nDOF constraint breakdown:\n');
    fprintf('  u-displacements: %d\n', supportInfo.constrainedDOFs(1));
    fprintf('  v-displacements: %d\n', supportInfo.constrainedDOFs(2));
    fprintf('  w-displacements: %d\n', supportInfo.constrainedDOFs(3));
    fprintf('  θx-rotations: %d\n', supportInfo.constrainedDOFs(4));
    fprintf('  θy-rotations: %d\n', supportInfo.constrainedDOFs(5));
    fprintf('_________________________________________________________\n\n');
end

end


function [affectedCPs] = findAffectedControlPoints(location, numCPs_xi, numCPs_eta)
%% Find control points affected by support

switch location.type
    case 'point'
        % Single control point
        xi_idx = round(location.xi * (numCPs_xi - 1)) + 1;
        eta_idx = round(location.eta * (numCPs_eta - 1)) + 1;
        xi_idx = max(1, min(numCPs_xi, xi_idx));
        eta_idx = max(1, min(numCPs_eta, eta_idx));
        affectedCPs = sub2ind([numCPs_xi, numCPs_eta], xi_idx, eta_idx);
        
    case 'corner'
        % Corner point (specific corner)
        switch lower(location.corner)
            case 'bottom-left'
                affectedCPs = sub2ind([numCPs_xi, numCPs_eta], 1, 1);
            case 'bottom-right'
                affectedCPs = sub2ind([numCPs_xi, numCPs_eta], numCPs_xi, 1);
            case 'top-left'
                affectedCPs = sub2ind([numCPs_xi, numCPs_eta], 1, numCPs_eta);
            case 'top-right'
                affectedCPs = sub2ind([numCPs_xi, numCPs_eta], numCPs_xi, numCPs_eta);
            otherwise
                error('Unknown corner specification: %s', location.corner);
        end
        
    case 'edge'
        % Edge support
        affectedCPs = [];
        
        if isfield(location, 'xi') && length(location.xi) == 2
            % xi-direction edge
            xi_start = round(location.xi(1) * (numCPs_xi - 1)) + 1;
            xi_end = round(location.xi(2) * (numCPs_xi - 1)) + 1;
            xi_indices = xi_start:xi_end;
            
            if isfield(location, 'eta_fixed')
                eta_idx = round(location.eta_fixed * (numCPs_eta - 1)) + 1;
                for xi_idx = xi_indices
                    affectedCPs = [affectedCPs, sub2ind([numCPs_xi, numCPs_eta], xi_idx, eta_idx)];
                end
            end
        end
        
        if isfield(location, 'eta') && length(location.eta) == 2
            % eta-direction edge
            eta_start = round(location.eta(1) * (numCPs_eta - 1)) + 1;
            eta_end = round(location.eta(2) * (numCPs_eta - 1)) + 1;
            eta_indices = eta_start:eta_end;
            
            if isfield(location, 'xi_fixed')
                xi_idx = round(location.xi_fixed * (numCPs_xi - 1)) + 1;
                for eta_idx = eta_indices
                    affectedCPs = [affectedCPs, sub2ind([numCPs_xi, numCPs_eta], xi_idx, eta_idx)];
                end
            end
        end
        
    case 'region'
        % Rectangular region
        xi_start = round(location.xi(1) * (numCPs_xi - 1)) + 1;
        xi_end = round(location.xi(2) * (numCPs_xi - 1)) + 1;
        eta_start = round(location.eta(1) * (numCPs_eta - 1)) + 1;
        eta_end = round(location.eta(2) * (numCPs_eta - 1)) + 1;
        
        affectedCPs = [];
        for xi_idx = xi_start:xi_end
            for eta_idx = eta_start:eta_end
                affectedCPs = [affectedCPs, sub2ind([numCPs_xi, numCPs_eta], xi_idx, eta_idx)];
            end
        end
        
    otherwise
        error('Unknown support location type: %s', location.type);
end

% Ensure indices are within bounds
affectedCPs = unique(affectedCPs);
affectedCPs = affectedCPs(affectedCPs <= numCPs_xi * numCPs_eta);

end


function [constraints] = getDefaultConstraints(supportType)
%% Get default constraint configuration for support type

switch lower(supportType)
    case 'fixed'
        % Fixed support: all DOFs constrained
        constraints.u = true;
        constraints.v = true;
        constraints.w = true;
        constraints.thetax = true;
        constraints.thetay = true;
        
    case 'pinned'
        % Pinned support: displacements fixed, rotations free
        constraints.u = true;
        constraints.v = true;
        constraints.w = true;
        constraints.thetax = false;
        constraints.thetay = false;
        
    case 'roller_x'
        % Roller in x-direction: u free, v,w fixed, rotations free
        constraints.u = false;
        constraints.v = true;
        constraints.w = true;
        constraints.thetax = false;
        constraints.thetay = false;
        
    case 'roller_y'
        % Roller in y-direction: v free, u,w fixed, rotations free
        constraints.u = true;
        constraints.v = false;
        constraints.w = true;
        constraints.thetax = false;
        constraints.thetay = false;
        
    case 'roller_z'
        % Roller in z-direction: w free, u,v fixed, rotations free
        constraints.u = true;
        constraints.v = true;
        constraints.w = false;
        constraints.thetax = false;
        constraints.thetay = false;
        
    case 'simply_supported'
        % Simply supported: w fixed, u,v,rotations free
        constraints.u = false;
        constraints.v = false;
        constraints.w = true;
        constraints.thetax = false;
        constraints.thetay = false;
        
    case 'clamped_edge'
        % Clamped edge: w and one rotation fixed
        constraints.u = false;
        constraints.v = false;
        constraints.w = true;
        constraints.thetax = true;
        constraints.thetay = false;
        
    case 'symmetry_x'
        % Symmetry about x-axis: v and θy constrained
        constraints.u = false;
        constraints.v = true;
        constraints.w = false;
        constraints.thetax = false;
        constraints.thetay = true;
        
    case 'symmetry_y'
        % Symmetry about y-axis: u and θx constrained
        constraints.u = true;
        constraints.v = false;
        constraints.w = false;
        constraints.thetax = true;
        constraints.thetay = false;
        
    case 'custom'
        % Custom constraints (to be specified elsewhere)
        constraints.u = false;
        constraints.v = false;
        constraints.w = false;
        constraints.thetax = false;
        constraints.thetay = false;
        
    otherwise
        error('Unknown support type: %s', supportType);
end

end