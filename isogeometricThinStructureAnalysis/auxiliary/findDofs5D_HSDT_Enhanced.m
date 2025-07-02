function [homDOFs, inhomDOFs, valuesInhomDOFs, supportInfo] = findDofs5D_HSDT_Enhanced...
    (BSplinePatch, boundaryConditions, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Enhanced DOF finding function for HSDT shell analysis with comprehensive
% boundary condition support for all 5 DOF including detailed support
% configurations and automatic constraint validation.
%
%                Input :
%         BSplinePatch : B-Spline patch structure containing control points
%                        and DOF numbering for 5 DOF per control point
%   boundaryConditions : Array of boundary condition structures:
%                        .type : Boundary condition type
%                          'dirichlet' - Prescribed displacements/rotations
%                          'neumann' - Prescribed forces/moments (for future)
%                          'support' - Support with constraint specifications
%                        .location : Location specification (see create3DSupports_HSDT)
%                        .constraints : Constraint specification for supports
%                          .u, .v, .w, .thetax, .thetay : true/false
%                        .values : Prescribed values for Dirichlet BC
%                          .u, .v, .w, .thetax, .thetay : scalar values
%                        .supportType : Predefined support type ('fixed', 'pinned', etc.)
%               outMsg : Output message flag
%
%               Output :
%              homDOFs : Array of homogeneous (zero) constraint DOF numbers
%            inhomDOFs : Array of inhomogeneous (non-zero) constraint DOF numbers
%       valuesInhomDOFs : Values for inhomogeneous constraints
%          supportInfo : Comprehensive support information structure
%
% Function layout :
%
% 1. Input validation and initialization
%
% 2. Process boundary conditions
% ->
%    2i. Homogeneous constraints (supports with zero displacement/rotation)
%    2ii. Inhomogeneous constraints (prescribed non-zero values)
%    2iii. Validate constraint compatibility
% <-
%
% 3. Generate detailed constraint summary
%
%% Function main body

if strcmp(outMsg, 'outputEnabled')
    fprintf('_________________________________________________________\n');
    fprintf('#########################################################\n');
    fprintf('Enhanced 5-DOF Boundary Condition Processing for HSDT\n');
    fprintf('DOF ordering: [u, v, w, θx, θy] per control point\n');
    fprintf('_________________________________________________________\n\n');
end

%% 1. Input validation and initialization

% Initialize output arrays
homDOFs = [];
inhomDOFs = [];
valuesInhomDOFs = [];

% Initialize support information
supportInfo = struct();
supportInfo.totalBCs = 0;
supportInfo.homogeneousConstraints = 0;
supportInfo.inhomogeneousConstraints = 0;
supportInfo.constrainedDOFs = [0, 0, 0, 0, 0]; % Count for each DOF type
supportInfo.boundaryConditions = [];

% Get mesh and DOF information
CP = BSplinePatch.CP;
DOFNumbering = BSplinePatch.DOFNumbering;
numCPs_xi = size(CP, 1);
numCPs_eta = size(CP, 2);
totalCPs = numCPs_xi * numCPs_eta;
totalDOFs = 5 * totalCPs;

% Validate DOF numbering
if size(DOFNumbering, 3) ~= 5
    error('DOF numbering must have 5 components for HSDT analysis');
end

% Initialize constraint tracking
constraintMatrix = false(totalCPs, 5); % Track which DOFs are constrained
constraintValues = zeros(totalCPs, 5);  % Track constraint values

if strcmp(outMsg, 'outputEnabled')
    fprintf('Mesh information:\n');
    fprintf('  Control points: %d x %d = %d\n', numCPs_xi, numCPs_eta, totalCPs);
    fprintf('  Total DOFs: %d (5 per control point)\n', totalDOFs);
    fprintf('  Processing %d boundary conditions...\n\n', length(boundaryConditions));
end

%% 2. Process boundary conditions

for iBC = 1:length(boundaryConditions)
    bc = boundaryConditions(iBC);
    
    if strcmp(outMsg, 'outputEnabled')
        fprintf('Processing BC %d/%d: Type = %s\n', iBC, length(boundaryConditions), bc.type);
    end
    
    switch lower(bc.type)
        case 'support'
            %% 2i. Process support boundary conditions
            
            % Use create3DSupports_HSDT for support processing
            if isfield(bc, 'supportType') && ~isempty(bc.supportType)
                supportType = bc.supportType;
            else
                supportType = 'custom';
            end
            
            % Create temporary support config
            tempConfig.location = bc.location;
            if isfield(bc, 'constraints')
                tempConfig.constraints = bc.constraints;
            end
            if isfield(bc, 'description')
                tempConfig.description = bc.description;
            else
                tempConfig.description = sprintf('Support BC %d', iBC);
            end
            
            % Process support
            [newHomDOFs, tempSupportInfo] = create3DSupports_HSDT...
                ([], tempConfig, CP, DOFNumbering, supportType, 'outputDisabled');
            
            % Add to homogeneous constraints
            homDOFs = [homDOFs, newHomDOFs];
            
            % Update constraint tracking
            for iSupport = 1:length(tempSupportInfo.supports)
                support = tempSupportInfo.supports(iSupport);
                for iCP = 1:length(support.affectedCPs)
                    cpIndex = support.affectedCPs(iCP);
                    if support.constraints.u
                        constraintMatrix(cpIndex, 1) = true;
                    end
                    if support.constraints.v
                        constraintMatrix(cpIndex, 2) = true;
                    end
                    if support.constraints.w
                        constraintMatrix(cpIndex, 3) = true;
                    end
                    if support.constraints.thetax
                        constraintMatrix(cpIndex, 4) = true;
                    end
                    if support.constraints.thetay
                        constraintMatrix(cpIndex, 5) = true;
                    end
                end
            end
            
            supportInfo.constrainedDOFs = supportInfo.constrainedDOFs + tempSupportInfo.constrainedDOFs;
            supportInfo.homogeneousConstraints = supportInfo.homogeneousConstraints + length(newHomDOFs);
            
        case 'dirichlet'
            %% 2ii. Process Dirichlet boundary conditions (prescribed values)
            
            % Find affected control points
            [affectedCPs] = findAffectedControlPoints(bc.location, numCPs_xi, numCPs_eta);
            
            newInhomDOFs = [];
            newInhomValues = [];
            constraintCount = [0, 0, 0, 0, 0];
            
            for iCP = 1:length(affectedCPs)
                cpIndex = affectedCPs(iCP);
                [cpi, cpj] = ind2sub([numCPs_xi, numCPs_eta], cpIndex);
                
                % Process u-displacement constraint
                if isfield(bc.values, 'u') && ~isnan(bc.values.u)
                    dofNum = DOFNumbering(cpi, cpj, 1);
                    if ~ismember(dofNum, [homDOFs, inhomDOFs])
                        newInhomDOFs = [newInhomDOFs, dofNum];
                        newInhomValues = [newInhomValues, bc.values.u];
                        constraintMatrix(cpIndex, 1) = true;
                        constraintValues(cpIndex, 1) = bc.values.u;
                        constraintCount(1) = constraintCount(1) + 1;
                    end
                end
                
                % Process v-displacement constraint
                if isfield(bc.values, 'v') && ~isnan(bc.values.v)
                    dofNum = DOFNumbering(cpi, cpj, 2);
                    if ~ismember(dofNum, [homDOFs, inhomDOFs])
                        newInhomDOFs = [newInhomDOFs, dofNum];
                        newInhomValues = [newInhomValues, bc.values.v];
                        constraintMatrix(cpIndex, 2) = true;
                        constraintValues(cpIndex, 2) = bc.values.v;
                        constraintCount(2) = constraintCount(2) + 1;
                    end
                end
                
                % Process w-displacement constraint
                if isfield(bc.values, 'w') && ~isnan(bc.values.w)
                    dofNum = DOFNumbering(cpi, cpj, 3);
                    if ~ismember(dofNum, [homDOFs, inhomDOFs])
                        newInhomDOFs = [newInhomDOFs, dofNum];
                        newInhomValues = [newInhomValues, bc.values.w];
                        constraintMatrix(cpIndex, 3) = true;
                        constraintValues(cpIndex, 3) = bc.values.w;
                        constraintCount(3) = constraintCount(3) + 1;
                    end
                end
                
                % Process θx-rotation constraint
                if isfield(bc.values, 'thetax') && ~isnan(bc.values.thetax)
                    dofNum = DOFNumbering(cpi, cpj, 4);
                    if ~ismember(dofNum, [homDOFs, inhomDOFs])
                        newInhomDOFs = [newInhomDOFs, dofNum];
                        newInhomValues = [newInhomValues, bc.values.thetax];
                        constraintMatrix(cpIndex, 4) = true;
                        constraintValues(cpIndex, 4) = bc.values.thetax;
                        constraintCount(4) = constraintCount(4) + 1;
                    end
                end
                
                % Process θy-rotation constraint
                if isfield(bc.values, 'thetay') && ~isnan(bc.values.thetay)
                    dofNum = DOFNumbering(cpi, cpj, 5);
                    if ~ismember(dofNum, [homDOFs, inhomDOFs])
                        newInhomDOFs = [newInhomDOFs, dofNum];
                        newInhomValues = [newInhomValues, bc.values.thetay];
                        constraintMatrix(cpIndex, 5) = true;
                        constraintValues(cpIndex, 5) = bc.values.thetay;
                        constraintCount(5) = constraintCount(5) + 1;
                    end
                end
            end
            
            % Add to inhomogeneous constraints
            inhomDOFs = [inhomDOFs, newInhomDOFs];
            valuesInhomDOFs = [valuesInhomDOFs, newInhomValues];
            
            supportInfo.constrainedDOFs = supportInfo.constrainedDOFs + constraintCount;
            supportInfo.inhomogeneousConstraints = supportInfo.inhomogeneousConstraints + length(newInhomDOFs);
            
        case 'neumann'
            %% 2iii. Process Neumann boundary conditions (for future implementation)
            
            if strcmp(outMsg, 'outputEnabled')
                fprintf('  Warning: Neumann BC processing not yet implemented\n');
            end
            
        otherwise
            error('Unknown boundary condition type: %s', bc.type);
    end
    
    %% Store boundary condition information
    
    bcInfo = struct();
    bcInfo.id = iBC;
    bcInfo.type = bc.type;
    bcInfo.location = bc.location;
    if isfield(bc, 'description')
        bcInfo.description = bc.description;
    else
        bcInfo.description = sprintf('%s BC %d', bc.type, iBC);
    end
    
    supportInfo.boundaryConditions = [supportInfo.boundaryConditions; bcInfo];
    supportInfo.totalBCs = supportInfo.totalBCs + 1;
    
    if strcmp(outMsg, 'outputEnabled')
        fprintf('  %s processed successfully\n', bcInfo.description);
    end
end

%% 3. Generate detailed constraint summary

% Remove duplicates and sort
homDOFs = sort(unique(homDOFs));
[inhomDOFs, uniqueIdx] = unique(inhomDOFs);
valuesInhomDOFs = valuesInhomDOFs(uniqueIdx);

% Calculate final statistics
supportInfo.constraintMatrix = constraintMatrix;
supportInfo.constraintValues = constraintValues;
supportInfo.totalConstrainedDOFs = length(homDOFs) + length(inhomDOFs);
supportInfo.totalFreeDOFs = totalDOFs - supportInfo.totalConstrainedDOFs;

% Validate constraints
[isValid, validationMsg] = validateConstraints(constraintMatrix, numCPs_xi, numCPs_eta);
supportInfo.isValid = isValid;
supportInfo.validationMessage = validationMsg;

if strcmp(outMsg, 'outputEnabled')
    fprintf('\n=== Enhanced 5-DOF Boundary Condition Summary ===\n');
    fprintf('Total boundary conditions: %d\n', supportInfo.totalBCs);
    fprintf('Total DOFs: %d\n', totalDOFs);
    fprintf('Constrained DOFs: %d (Homogeneous: %d, Inhomogeneous: %d)\n', ...
            supportInfo.totalConstrainedDOFs, length(homDOFs), length(inhomDOFs));
    fprintf('Free DOFs: %d\n', supportInfo.totalFreeDOFs);
    fprintf('\nConstraint breakdown by DOF type:\n');
    fprintf('  u-displacements: %d\n', supportInfo.constrainedDOFs(1));
    fprintf('  v-displacements: %d\n', supportInfo.constrainedDOFs(2));
    fprintf('  w-displacements: %d\n', supportInfo.constrainedDOFs(3));
    fprintf('  θx-rotations: %d\n', supportInfo.constrainedDOFs(4));
    fprintf('  θy-rotations: %d\n', supportInfo.constrainedDOFs(5));
    
    fprintf('\nConstraint validation: ');
    if isValid
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
        fprintf('  %s\n', validationMsg);
    end
    fprintf('_________________________________________________________\n\n');
end

end


function [affectedCPs] = findAffectedControlPoints(location, numCPs_xi, numCPs_eta)
%% Find control points affected by boundary condition
% (Reuse from create3DSupports_HSDT)

switch location.type
    case 'point'
        % Single control point
        xi_idx = round(location.xi * (numCPs_xi - 1)) + 1;
        eta_idx = round(location.eta * (numCPs_eta - 1)) + 1;
        xi_idx = max(1, min(numCPs_xi, xi_idx));
        eta_idx = max(1, min(numCPs_eta, eta_idx));
        affectedCPs = sub2ind([numCPs_xi, numCPs_eta], xi_idx, eta_idx);
        
    case 'corner'
        % Corner point
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
        % Edge boundary condition
        affectedCPs = [];
        
        if isfield(location, 'side')
            % Named edge specification
            switch lower(location.side)
                case 'bottom'
                    for xi_idx = 1:numCPs_xi
                        affectedCPs = [affectedCPs, sub2ind([numCPs_xi, numCPs_eta], xi_idx, 1)];
                    end
                case 'top'
                    for xi_idx = 1:numCPs_xi
                        affectedCPs = [affectedCPs, sub2ind([numCPs_xi, numCPs_eta], xi_idx, numCPs_eta)];
                    end
                case 'left'
                    for eta_idx = 1:numCPs_eta
                        affectedCPs = [affectedCPs, sub2ind([numCPs_xi, numCPs_eta], 1, eta_idx)];
                    end
                case 'right'
                    for eta_idx = 1:numCPs_eta
                        affectedCPs = [affectedCPs, sub2ind([numCPs_xi, numCPs_eta], numCPs_xi, eta_idx)];
                    end
                otherwise
                    error('Unknown edge specification: %s', location.side);
            end
        else
            % Parametric edge specification (reuse from create3DSupports_HSDT)
            if isfield(location, 'xi') && length(location.xi) == 2
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
        error('Unknown location type: %s', location.type);
end

% Ensure indices are within bounds
affectedCPs = unique(affectedCPs);
affectedCPs = affectedCPs(affectedCPs <= numCPs_xi * numCPs_eta);

end


function [isValid, validationMsg] = validateConstraints(constraintMatrix, numCPs_xi, numCPs_eta)
%% Validate constraint configuration

isValid = true;
validationMsg = '';

% Check for rigid body motion prevention
totalCPs = numCPs_xi * numCPs_eta;
constrainedCount = sum(constraintMatrix, 1);

% HSDT shell requires sufficient constraints to prevent rigid body motion
% Minimum requirements:
% - At least 3 points with w constrained (prevent translation in z)
% - At least 2 points with u constrained (prevent translation in x)  
% - At least 2 points with v constrained (prevent translation in y)
% - Some rotational constraints to prevent rigid body rotations

if constrainedCount(3) < 3  % w constraints
    isValid = false;
    validationMsg = [validationMsg, 'Insufficient w constraints (need >=3). '];
end

if constrainedCount(1) < 2  % u constraints
    isValid = false;
    validationMsg = [validationMsg, 'Insufficient u constraints (need >=2). '];
end

if constrainedCount(2) < 2  % v constraints
    isValid = false;
    validationMsg = [validationMsg, 'Insufficient v constraints (need >=2). '];
end

% Check for over-constraints in small structures
if totalCPs <= 4
    totalConstraints = sum(constraintMatrix(:));
    if totalConstraints > 0.8 * 5 * totalCPs
        validationMsg = [validationMsg, 'Possible over-constraint in small structure. '];
    end
end

if isValid
    validationMsg = 'All constraint validation checks passed';
end

end