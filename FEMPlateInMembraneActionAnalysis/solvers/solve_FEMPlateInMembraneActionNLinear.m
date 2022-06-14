function [dHat,FComplete,minElSize] = ...
    solve_FEMPlateInMembraneActionNLinear ...
    (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    propNBC, computeBodyForces, parameters, solve_LinearSystem, ...
    propNLinearAnalysis, propGaussInt, propOutput, caseName, pathToOutput, ...
    outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displacement field, the complete force vector and the minimum
% element area size for a plate in membrane action problem solved with the
% classical Finite Element Method for the nonlinear case.
%
%              Input :
%       propAnalysis : Structure containing information about the analysis
%                           .type : The analysis type
%             strMsh : Nodes and elements in the mesh
%               dHat : Initial conditions
%            homDOFs : The global numbering of the nodes where homogeneous
%                      Dirichlet boundary conditions are applied 
%          inhomDOFs : The global numbering of the nodes where 
%                      inhomogeneous Dirichlet boundary conditions are 
%                      applied
%    valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                      Dirichlet boundary conditions are applied
%            propNBC : Structure containing information on the Neumann
%                      boundary conditions
%                         .nodes : The nodes where Neumann boundary 
%                                  conditions are applied
%                      .loadType : The type of the load for each Neumann 
%                                  node
%                     .fctHandle : The function handle for each Neumann 
%                                  node for the computation of the load 
%                                  vector (these functions are unde the 
%                                  folder load)
%  computeBodyForces : Function handle to body force vector computation
%         parameters : Problem specific technical parameters
% solve_LinearSystem : Function handle to the solution of the 
%    nLinearAnalysis :     .scheme : The employed nonlinear scheme
%                       .tolerance : The residual tolerance
%                         .maxIter : The maximum number of the nonlinear 
%                                    iterations
%       propGaussInt : On the numerical integration (quadrature)
%                       .type : 'default', 'manual'
%                       .noGP : Number of Gauss Points
%         propOutput : Structure containing information on writting the
%                      results for postprocessing,
%                                .isOutput : Flag on whether the results 
%                                            to be written out
%                       .writeOutputToFile : Function handle to the
%                                            writting out of the results
%                           .VTKResultFile : Specifies the name of the
%                                            VTK result file from which
%                                            the simulation to be restarted
%                                            If it is specified as 
%                                            'undefined' the simulation 
%                                            starts from time TStart
%           caseName : The name of the case in the inputGiD case folder
%       pathToOutput : Define the path to where to write out the results
%             outMsg : On outputting information
%
%             Output :
%               dHat : The nodal displacement field
%          FComplete : The complete force vector
%          minElSize : The minimum element area size in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Find the prescribed and the free DOFs of the system
%
% 2. Solve the nonlinear equation system
%
% 3. Postprocessing
%
% 4. Write out the results into a file
%
% 5. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________________\n');
    fprintf('###################################################################\n');
    fprintf('Computation of the displacement field for a geometrically nonlinear\n');
    fprintf('plate in membrane action problem has been initiated\n');
    fprintf('___________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of nodes in the mesh
numNodes = length(strMsh.nodes(:,1));

% Number of DOFs in the mesh
numDOFs = 2*numNodes;

% GLobal DOF numbering
DOFNumbering = 1:numDOFs;

% Assign dummy variables
uSaved = 'undefined';
uDot = 'undefined';
uDDot = 'undefined';
uDotSaved = 'undefined';
uDDotSaved = 'undefined';
% propStrDynamics.timeDependence = 'STEADY_STATE';
propStrDynamics = 'undefined';
uMeshALE = 'undefined';
dDot = 'undefined';
dDDot = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
precompStiffMtx = 'undefined';
precomResVct = 'undefined';
numTimeStep = 0;

% Title for the output file
title = 'geometrically linear steady-state plane stress analysis';

% Initialize time
t = 0;

% Define tabulation for the output in the command window
tab = '\t';

%% 1. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs,inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = DOFNumbering;
freeDOFs(ismember(freeDOFs,prescribedDoFs)) = [];

%% 2. Compute load vector corresponding to a conservative load
F = computeLoadVctFEMPlateInMembraneAction...
    (strMsh,propAnalysis,propNBC,t,propGaussInt,'');

%% 3. Solve the linear equation system
[dHat,~,FComplete,minElSize] = solve_FEMNLinearSystem ...
    (propAnalysis, uSaved, uDotSaved, uDDotSaved, strMsh, F, ...
    computeBodyForces, parameters, dHat, uDot, uDDot, massMtx, dampMtx, ...
    precompStiffMtx, precomResVct, ...
    @computeTangentStiffMtxResVctFEMPlateInMembraneActionCST, ...
    DOFNumbering, freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    uMeshALE, solve_LinearSystem, propStrDynamics, t, ...
    propNLinearAnalysis, propGaussInt, tab, outMsg);

%% 4. Write out the results into a file
if isfield(propOutput, 'isOutput')
    if isa(propOutput.isOutput, 'logical')
        if propOutput.isOutput
            if isfield(propOutput, 'writeOutputToFile')
                if isa(propOutput.writeOutputToFile, 'function_handle')
                    fprintf('>> Writting out the results to "%s"\n',strcat(pathToOutput, caseName,'/'));
                    DOF4Output = [1:2:numDOFs - 1
                                  2:2:numDOFs];
                    writeOutputFEMPlateInMembraneActionToVTK(propAnalysis, ...
                        propNLinearAnalysis, propStrDynamics, strMsh, ...
                        parameters, dHat, dDot, dDDot, DOF4Output, ...
                        caseName, pathToOutput, title, numTimeStep);
                    else
                    error('Variable propVTK.writeOutputToFile should define a function handle');
                end
            else
                error('Structure propVTK should define variable writeOutputToFile');
            end
        end
    else
        error('Structure propVTK should define boolean isOutput');
    end
end

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nNonlinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('_______________________Nonlinear Analysis Ended____________________\n');
    fprintf('####################################################################\n\n\n');
end

end
