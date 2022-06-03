%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Plane stress analysis for a rectangular plate subject to uniform
%        pressure on its top edge
%
% Date : 19.02.2014
%
%% Preamble

% Clear memory
clear;

% Clear the command window
clc;

% Close all generated windows
close all;
 
%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the low order basis functions
addpath('../../basisFunctions/');

% Add the equation system solvers
addpath('../../equationSystemSolvers/');

% Add all the efficient computation functions
addpath('../../efficientComputation/');

% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMPlateInMembraneActionAnalysis/';
caseName = 'cantileverBeamPlaneStress';
% caseName = 'PlateWithAHolePlaneStress';
% caseName = 'PlateWithMultipleHolesPlaneStress';
% caseName = 'InfinitePlateWithAHolePlaneStress';
% caseName = 'unitTest_curvedPlateTipShearPlaneStress';
% caseName = 'NACA2412_AoA5_CSD';
% caseName = 'turek_csd';

% Parse the data from the GiD input file
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, propStrDynamics, propGaussInt] = ...
    parse_StructuralModelFromGid(pathToCase, caseName, 'outputEnabled');

%% GUI

% On the body forces
computeBodyForces = @computeConstantVerticalStructureBodyForceVct;

% Choose solver for the linear equation system
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Output properties
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMPlateInMembraneActionToVTK;
propVTK.VTKResultFile = 'undefined';

% Define traction vector
propNBC.tractionLoadVct = [0; -1e1; 0];

% Initialize graphics index
propGraph.index = 1;

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Visualization of the configuration
F = computeLoadVctFEMPlateInMembraneAction...
    (strMsh, propAnalysis, propNBC, 0, propGaussInt,'');
segmentsContact = [];
propGraph.index = plot_referenceConfigurationFEMPlateInMembraneAction ...
    (strMsh, propAnalysis, F, homDOFs, segmentsContact, propGraph, 'outputEnabled');

%% Initialize solution
numNodes = length(strMsh.nodes(:,1));
numDOFs = 2*numNodes;
dHat = zeros(numDOFs,1);

%% Solve the plate in membrane action problem
[dHat, FComplete, minElSize] = ...
    solve_FEMPlateInMembraneActionNLinear...
    (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, ...
    computeBodyForces, parameters, solve_LinearSystem, propNLinearAnalysis, ...
    propGaussInt, propVTK, caseName, pathToOutput, 'outputEnabled');

%% Postprocessing
% graph.visualization.geometry = 'reference_and_current';
% resultant = 'stress';
% component = 'y';
% nodeIDs_active = [];
% contactSegments = [];
% graph.index = plot_currentConfigurationAndResultants ...
%     (propAnalysis, strMsh, homDOFs, dHat, nodeIDs_active, ...
%     contactSegments, parameters, resultant, component, graph);

%% END OF THE SCRIPT
