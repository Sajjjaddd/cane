function [boundaryConditions] = createCommonSupports_HSDT(supportType, varargin)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Convenience function for creating common support configurations for HSDT 
% shell analysis. Provides pre-defined support patterns for typical structural
% problems including Scordelis-Lo roof, simply supported plates, etc.
%
%                Input :
%          supportType : String specifying the common support configuration:
%                        'scordelis_lo' - Standard Scordelis-Lo roof supports
%                        'simply_supported_plate' - Simply supported rectangular plate
%                        'clamped_plate' - Fully clamped rectangular plate
%                        'cantilever_beam' - Cantilever beam support
%                        'three_point_support' - Three point support (minimum)
%                        'four_point_support' - Four corner point supports
%                        'fixed_edge' - One edge completely fixed
%                        'pinned_edge' - One edge pinned
%                        'symmetry_quarter' - Quarter model with symmetry BCs
%                        'custom' - Custom support specification
%             varargin : Additional parameters depending on support type
%
%               Output :
%   boundaryConditions : Array of boundary condition structures ready for
%                        use with findDofs5D_HSDT_Enhanced
%
% Function layout :
%
% 1. Input validation
%
% 2. Create support configurations based on type
%
% 3. Return boundary condition structures
%
%% Function main body

% Parse input arguments
p = inputParser;
addRequired(p, 'supportType', @ischar);
addParameter(p, 'description', '', @ischar);
addParameter(p, 'customConstraints', [], @isstruct);
addParameter(p, 'locations', [], @iscell);
addParameter(p, 'values', [], @isstruct);
parse(p, supportType, varargin{:});

% Initialize boundary conditions array
boundaryConditions = [];

%% Create support configurations based on type

switch lower(supportType)
    
    case 'scordelis_lo'
        %% Scordelis-Lo Roof Standard Supports
        % Two curved edges are free, two straight edges have supports
        
        % Support 1: Left edge (xi=0) - diaphragm support
        bc1.type = 'support';
        bc1.location.type = 'edge';
        bc1.location.side = 'left';
        bc1.supportType = 'custom';
        bc1.constraints.u = false;  % Free longitudinal movement
        bc1.constraints.v = true;   % Constrain transverse movement
        bc1.constraints.w = true;   % Constrain vertical movement
        bc1.constraints.thetax = false;  % Free rotation about x
        bc1.constraints.thetay = false; % Free rotation about y
        bc1.description = 'Scordelis-Lo left edge diaphragm';
        
        % Support 2: Right edge (xi=1) - diaphragm support
        bc2.type = 'support';
        bc2.location.type = 'edge';
        bc2.location.side = 'right';
        bc2.supportType = 'custom';
        bc2.constraints.u = false;  % Free longitudinal movement
        bc2.constraints.v = true;   % Constrain transverse movement
        bc2.constraints.w = true;   % Constrain vertical movement
        bc2.constraints.thetax = false;  % Free rotation about x
        bc2.constraints.thetay = false; % Free rotation about y
        bc2.description = 'Scordelis-Lo right edge diaphragm';
        
        boundaryConditions = [bc1, bc2];
        
    case 'simply_supported_plate'
        %% Simply Supported Rectangular Plate
        % All edges have w constrained, all other DOFs free
        
        % Bottom edge
        bc1.type = 'support';
        bc1.location.type = 'edge';
        bc1.location.side = 'bottom';
        bc1.supportType = 'simply_supported';
        bc1.description = 'Simply supported bottom edge';
        
        % Top edge
        bc2.type = 'support';
        bc2.location.type = 'edge';
        bc2.location.side = 'top';
        bc2.supportType = 'simply_supported';
        bc2.description = 'Simply supported top edge';
        
        % Left edge
        bc3.type = 'support';
        bc3.location.type = 'edge';
        bc3.location.side = 'left';
        bc3.supportType = 'simply_supported';
        bc3.description = 'Simply supported left edge';
        
        % Right edge
        bc4.type = 'support';
        bc4.location.type = 'edge';
        bc4.location.side = 'right';
        bc4.supportType = 'simply_supported';
        bc4.description = 'Simply supported right edge';
        
        boundaryConditions = [bc1, bc2, bc3, bc4];
        
    case 'clamped_plate'
        %% Fully Clamped Rectangular Plate
        % All edges have all DOFs constrained
        
        % Bottom edge
        bc1.type = 'support';
        bc1.location.type = 'edge';
        bc1.location.side = 'bottom';
        bc1.supportType = 'fixed';
        bc1.description = 'Clamped bottom edge';
        
        % Top edge
        bc2.type = 'support';
        bc2.location.type = 'edge';
        bc2.location.side = 'top';
        bc2.supportType = 'fixed';
        bc2.description = 'Clamped top edge';
        
        % Left edge
        bc3.type = 'support';
        bc3.location.type = 'edge';
        bc3.location.side = 'left';
        bc3.supportType = 'fixed';
        bc3.description = 'Clamped left edge';
        
        % Right edge
        bc4.type = 'support';
        bc4.location.type = 'edge';
        bc4.location.side = 'right';
        bc4.supportType = 'fixed';
        bc4.description = 'Clamped right edge';
        
        boundaryConditions = [bc1, bc2, bc3, bc4];
        
    case 'cantilever_beam'
        %% Cantilever Beam
        % One edge fully fixed, all other edges free
        
        bc1.type = 'support';
        bc1.location.type = 'edge';
        bc1.location.side = 'left';
        bc1.supportType = 'fixed';
        bc1.description = 'Cantilever fixed end';
        
        boundaryConditions = bc1;
        
    case 'three_point_support'
        %% Minimum Three Point Support
        % Three corner points to prevent rigid body motion
        
        % Bottom-left corner - fixed
        bc1.type = 'support';
        bc1.location.type = 'corner';
        bc1.location.corner = 'bottom-left';
        bc1.supportType = 'fixed';
        bc1.description = 'Three-point support corner 1';
        
        % Bottom-right corner - pinned
        bc2.type = 'support';
        bc2.location.type = 'corner';
        bc2.location.corner = 'bottom-right';
        bc2.supportType = 'pinned';
        bc2.description = 'Three-point support corner 2';
        
        % Top-left corner - w only
        bc3.type = 'support';
        bc3.location.type = 'corner';
        bc3.location.corner = 'top-left';
        bc3.supportType = 'simply_supported';
        bc3.description = 'Three-point support corner 3';
        
        boundaryConditions = [bc1, bc2, bc3];
        
    case 'four_point_support'
        %% Four Corner Point Supports
        % All corners pinned
        
        % Bottom-left corner
        bc1.type = 'support';
        bc1.location.type = 'corner';
        bc1.location.corner = 'bottom-left';
        bc1.supportType = 'pinned';
        bc1.description = 'Four-point support corner 1';
        
        % Bottom-right corner
        bc2.type = 'support';
        bc2.location.type = 'corner';
        bc2.location.corner = 'bottom-right';
        bc2.supportType = 'pinned';
        bc2.description = 'Four-point support corner 2';
        
        % Top-left corner
        bc3.type = 'support';
        bc3.location.type = 'corner';
        bc3.location.corner = 'top-left';
        bc3.supportType = 'pinned';
        bc3.description = 'Four-point support corner 3';
        
        % Top-right corner
        bc4.type = 'support';
        bc4.location.type = 'corner';
        bc4.location.corner = 'top-right';
        bc4.supportType = 'pinned';
        bc4.description = 'Four-point support corner 4';
        
        boundaryConditions = [bc1, bc2, bc3, bc4];
        
    case 'fixed_edge'
        %% One Edge Completely Fixed
        
        bc1.type = 'support';
        bc1.location.type = 'edge';
        bc1.location.side = 'bottom';
        bc1.supportType = 'fixed';
        bc1.description = 'Fixed edge support';
        
        boundaryConditions = bc1;
        
    case 'pinned_edge'
        %% One Edge Pinned
        
        bc1.type = 'support';
        bc1.location.type = 'edge';
        bc1.location.side = 'bottom';
        bc1.supportType = 'pinned';
        bc1.description = 'Pinned edge support';
        
        boundaryConditions = bc1;
        
    case 'symmetry_quarter'
        %% Quarter Model with Symmetry Boundary Conditions
        % Two edges have symmetry conditions
        
        % Symmetry about x-axis (bottom edge)
        bc1.type = 'support';
        bc1.location.type = 'edge';
        bc1.location.side = 'bottom';
        bc1.supportType = 'symmetry_x';
        bc1.description = 'Symmetry about x-axis';
        
        % Symmetry about y-axis (left edge)
        bc2.type = 'support';
        bc2.location.type = 'edge';
        bc2.location.side = 'left';
        bc2.supportType = 'symmetry_y';
        bc2.description = 'Symmetry about y-axis';
        
        boundaryConditions = [bc1, bc2];
        
    case 'custom'
        %% Custom Support Configuration
        % Use provided constraints and locations
        
        if isempty(p.Results.locations)
            error('Custom support type requires locations parameter');
        end
        
        locations = p.Results.locations;
        customConstraints = p.Results.customConstraints;
        
        for i = 1:length(locations)
            bc.type = 'support';
            bc.location = locations{i};
            
            if ~isempty(customConstraints)
                bc.constraints = customConstraints;
                bc.supportType = 'custom';
            else
                bc.supportType = 'fixed';  % Default to fixed
            end
            
            if ~isempty(p.Results.description)
                bc.description = sprintf('%s %d', p.Results.description, i);
            else
                bc.description = sprintf('Custom support %d', i);
            end
            
            if i == 1
                boundaryConditions = bc;
            else
                boundaryConditions = [boundaryConditions, bc];
            end
        end
        
    otherwise
        error('Unknown support type: %s', supportType);
end

% Add global description if provided
if ~isempty(p.Results.description)
    for i = 1:length(boundaryConditions)
        if ~isfield(boundaryConditions(i), 'description') || isempty(boundaryConditions(i).description)
            boundaryConditions(i).description = sprintf('%s %d', p.Results.description, i);
        end
    end
end

end