function [displacement,lagrange] = solveSignoriniLagrange_2...
    (mesh,homDBC,contactNodes,F, segments,materialProperties,analysis,maxIteration,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%                  -------------------
%                  Fabien Pean
%                  Andreas Hauso
%                  Georgios Koroniotis
%
%% Function documentation
%
% Returns the displacement field and the Lagrange multipliers corresponding 
% to a plain stress/strain analysis for the given mesh and geometry 
% together with its Dirichlet and Neumann boundary conditions and the 
% contact constraints for multiple rigids walls by applying the Lagrange
% multiplier method.
% 
%              Input :
%               mesh : Elements and nodes of the mesh
%             homDBC : Vector of the prescribed DoFs (by global numbering)
%       contactNodes : structure containing the global numbering of the
%                      canditate contact nodes
%                  F : Global load vector
%           segments : Matrix with the coordinates of two wall determining
%                      points, for every segment j=1..n
% materialProperties : The material properties of the structure
%           analysis : Structure about the analysis type
%       maxIteration : Maximum number of iterations
%             outMsg : write 'outputEnabled' to show output during runtime
%
%             Output :
%       displacement : The resulting displacement field
%           lagrange : .multipliers  : values of the Lagrange multipliers
%                    : .active_nodes : node numbers of the active nodes
%
% Function layout :
%
% 0. Remove fully constrained nodes
%
% 1. Compute the gap function
%
% 2. Compute the master stiffness matrix of the structure
%
% 3. Reduce the system according to the given constraints
%|->
% 4. Solve the system according to Lagrange multipliers
%   4.1 Assemble to the complete displacement vector
%   4.2 Detect active nodes
%   4.3 Rebuild system if new active nodes found
%   4.4 Relax the system until ONLY valid Lagrange multipliers are computed
%   |->
%   4.4.1 Compute the displacement and Lagrange multipliers
%   4.4.2 Detect and delete non-valid Lagrange multipliers and related rows
%   <-|
%<-|
% 5. compute the complete load vector and verify the results
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    if strcmp(analysis.type,'planeStress')
        fprintf('Plane stress analysis has been initiated \n');
    elseif strcmp(analysis.type,'planeStrain')
        fprintf('Plane strain analysis has been initiated \n');
    end
    fprintf('\n');
end

%% 0. Remove fully constrained nodes

% Remove fully constrained nodes from the tests
% from number of indices until 1
for i=size(contactNodes.indices,1):-1:1
    % Determine how many Dirichlet conditions correspond to the node:  
    nodeHasDirichlet=ismember(floor((homDBC+1)/2),contactNodes.indices(i));
    numberOfDirichlet=length(nodeHasDirichlet(nodeHasDirichlet==1));
    % If the 2D node has at least two Dirichlet conditions exclude it from the contact canditates :  
    if (numberOfDirichlet>=2)
       contactNodes.indices(i)=[];
    end
end
contactNodes.positions=mesh.nodes(contactNodes.indices,:);

%% 1. Compute the gap function

% Compute normal, parallel and position vectors of the segments
segments = buildSegmentsData(segments);

% Compute for all nodes the specific gap and save it in the field .gap
contactNodes = computeGapFunction(contactNodes,segments);

%% 2. Compute the master stiffness matrix of the structure
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Computing master stiffness matrix... \n');
end

% Get number of DOFs in original system
nDOFs = length(F);
% Master global stiffness matrix
K = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,materialProperties,analysis);

%% 3. Reduce the system according to the given constraints
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Reducing the system according to the constraints... \n');
end

K_red = K;
F_red = F;

% Remove constrained DOFs (homDBCs)
K_red(:,homDBC) = [];
K_red(homDBC,:) = [];
F_red(homDBC) = [];

%% 4. Solve the system according to Lagrange multipliers

% Initialize the displacement vector
displacement_red = zeros(length(F_red),1);

% Initial values for the iteration
isCndMain = true;
activeNodes = [];
nActiveNodes = 0;
it = 0;

% Counts the number of equations which are solved during iteration   
equations_counter = 0;
    
% Iterate until no more invalid Lagrange multipliers AND no new active nodes
% are added in the pool AND max number of iterations in not reached
while(isCndMain && it<maxIteration)    

    %% 4.1 Assemble to the complete displacement vector
    displacement = buildFullDisplacement(nDOFs,homDBC,displacement_red);

    %% 4.2 Detect active nodes
    activeNodes_tmp = detectActiveDOFs(contactNodes,displacement,segments);

    if(isequaln(activeNodes_tmp,activeNodes) && it~=0)
        isCndMain = false;
    end

    %% 4.3 Rebuild system if new active nodes found
    if(~isempty(activeNodes_tmp) && isCndMain)

        % Update number of active nodes
        activeNodes = activeNodes_tmp;
        nActiveNodes = length(activeNodes);

        % Build constraint matrix C and rhs F
        C  = buildConstraintMatrix(nDOFs,contactNodes,activeNodes,segments);
        F_red = buildRHS(F,contactNodes,activeNodes,segments);

        % Build master system matrix
        K_red = [K,C;C',zeros(size(C,2))];

        % Reduce the system according to the BCs
        K_red(:,homDBC) = [];
        K_red(homDBC,:) = [];
        F_red(homDBC) = [];
    end
    
    %% 4.4 Relax the system until ONLY valid Lagrange multipliers are computed
    isCndLagrange = true;

    while(isCndLagrange && it<maxIteration)

        it = it + 1;
        isCndLagrange = false;

        %% 4.4.1 Compute the displacement and Lagrange multipliers
        
        % count the number of total solved equations
        equations_counter = equations_counter + length(K_red);
        
        if strcmp(outMsg,'outputEnabled')
            fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(K_red),cond(K_red));
        end
        
        % Solve using the backslash operator
        displacement_red = K_red\F_red;
        
        %% 4.4.2 Detect and delete non-valid Lagrange multipliers and related rows
        
        % Lagrange Multipliers indices
        lmDOFsIndices = length(displacement_red)-nActiveNodes+1:length(displacement_red);
        
        if max(displacement_red(lmDOFsIndices)) > 0
            isCndLagrange = true;
            isCndMain = true;
        end
        
        % Find the indices of only non-compressive Lagrange Multipliers
        lmDOFsIndices = lmDOFsIndices(displacement_red(lmDOFsIndices)>=0);
       
        % Delete non-valid Lagrange multipliers and related rows
        K_red(:,lmDOFsIndices) = [];
        K_red(lmDOFsIndices,:) = [];
        F_red(lmDOFsIndices) = [];
        
        % Delete active nodes related to non-valid Lagrange multipliers
        activeNodes(lmDOFsIndices-length(displacement_red)+nActiveNodes) = [];
        
        % Update the number of active nodes
        nActiveNodes=length(activeNodes);

    end % end of inner while loop

end % end of main while loop

%% 5. Compute the complete load vector and verify the results

% Select and save node numbers of active nodes
allContactNodes=[];
for j=1:segments.number
    allContactNodes=[allContactNodes;contactNodes.indices];
end  

% Build full displacement vector
displacement = buildFullDisplacement(nDOFs,homDBC,displacement_red);

% Keep only lagrange multipliers of the active nodes
lagrangeIndices = length(displacement_red)-nActiveNodes+1:length(displacement_red);
lagrange.multipliers = displacement_red(lagrangeIndices);
lagrange.active_nodes = allContactNodes(activeNodes);

%% 6. Print info
if strcmp(outMsg,'outputEnabled')
    % energy of the structure
    energy = displacement'*K*displacement;
    % output
    fprintf('\n');
    fprintf('Output informations...\n');
    fprintf('\t Constraints solved in %d iterations. A total of %d equations were solved. \n',it,equations_counter);
    fprintf('\t %d active nodes found.\n',length(lagrange.active_nodes));
    fprintf('\t #DOF: %d \n\t Energy norm of the structure: %4.2f\n',nDOFs,energy);
    if it >= maxIteration
        fprintf('\t Max number of iterations of has been reached !! Not Converged !!\n');
    end
end

end