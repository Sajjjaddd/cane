%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Run unitTest cases for the following type of analysis,
%        - Isogeometric beam analysis
%        - Isogeometric membrane analysis
%        - Isogeometric Kichhoff-Love shell analysis
%        - Finite element formulation for plate in membrane action analysis
%        - Stabilized finite element analysis for the Navier-Stokes problem
%
% Date : 12.02.2015
%
%% Clear memnory and command window
clc;
clear;
 
%% Includes

% Add transient analysis functions
addpath('../../transientAnalysis/');

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Add all fast computation functions
addpath('../../efficientComputation/');

% Add the equation system solvers
addpath('../../equationSystemSolvers/');

% Add the classical finite element basis functions
addpath('../../basisFunctions/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');

% Add all functions related to the isogeometric beam analysis
addpath('../../isogeometricBeamAnalysis/stiffnessMatrices/',...
        '../../isogeometricBeamAnalysis/graphics/',...
        '../../isogeometricBeamAnalysis/loads/',...
        '../../isogeometricBeamAnalysis/postprocessing/',...
        '../../isogeometricBeamAnalysis/solvers/',...
        '../../isogeometricBeamAnalysis/math/',...
        '../../isogeometricBeamAnalysis/auxiliary/',...
        '../../isogeometricBeamAnalysis/errorComputation/');

% Add all functions related to the isogeometric Kirchhoff-Love shell formulation
addpath('../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../isogeometricThinStructureAnalysis/graphicsMultipatches/',...
        '../../isogeometricThinStructureAnalysis/loads/',...
        '../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../isogeometricThinStructureAnalysis/solvers/',...
        '../../isogeometricThinStructureAnalysis/metrics/',...
        '../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../isogeometricThinStructureAnalysis/penaltyDecompositionKLShell/',...
        '../../isogeometricThinStructureAnalysis/penaltyDecompositionMembrane/',...
        '../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionKLShell/',...
        '../../isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/',...
        '../../isogeometricThinStructureAnalysis/errorComputation/',...
        '../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../isogeometricThinStructureAnalysis/formFindingAnalysis/');
    
% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/',...
        '../../FEMPlateInMembraneActionAnalysis/initialConditions/');
    
% Add all functions related to the Finite Element Methods for Computational
% Fluid Dynamics problems
addpath('../../FEMComputationalFluidDynamicsAnalysis/solutionMatricesAndVectors/',...
        '../../FEMComputationalFluidDynamicsAnalysis/initialConditions',...
        '../../FEMComputationalFluidDynamicsAnalysis/solvers/',...
        '../../FEMComputationalFluidDynamicsAnalysis/loads/',...
        '../../FEMComputationalFluidDynamicsAnalysis/output/',...
        '../../FEMComputationalFluidDynamicsAnalysis/ALEMotion/',...
        '../../FEMComputationalFluidDynamicsAnalysis/transientAnalysis/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all unit test functions and classes
addpath('../../unitTest/');

% Add all functions for Monte Carlo-related functionality
addpath('../../unitTestSupportFunctions/');
    
%% Global variables for the unit tests
isLight = true;

%% Import the matlab test suite
import matlab.unittest.TestSuite;
import matlab.unittest.selectors.HasName;
import matlab.unittest.constraints.EndsWithSubstring;

%% Run the unit test cases for the quadratures
suiteQuadrature = TestSuite.fromClass(?testQuadrature);
resultQuadrature = run(suiteQuadrature);

%% Run the unit test cases for the utilities functions
suiteClassUtilityFunctions = TestSuite.fromClass(?testUtilityFunctions);
resultIGAUtilityFunctions = run(suiteClassUtilityFunctions);

%% Run the unit test cases for the isogeometric beam analysis
suiteClassIGABeam = TestSuite.fromClass(?testIGABeamAnalysis);
resultIGABeam = run(suiteClassIGABeam);

%% Run the unit test cases for the isogeometric membrane analysis
suiteClassIGAMembrane = TestSuite.fromClass(?testIGAMembraneAnalysis);
if isLight
    suiteClassIGAMembrane = suiteClassIGAMembrane.selectIf...
        (HasName('testIGAMembraneAnalysis/testPagewiseComputationTangStiffMtxIGAMembrane') | ...
        HasName('testIGAMembraneAnalysis/testNitscheWeakDBCPlateMembraneAnalysis') | ...
        HasName('testIGAMembraneAnalysis/testTransientNonlinearMembraneAnalysis') | ...
        HasName('testIGAMembraneAnalysis/testTransientNitscheWeakDBCPlateMembraneAnalysis') | ...
        HasName('testIGAMembraneAnalysis/testDDMNitschePlateMembraneAnalysis') | ...
        HasName('testIGAMembraneAnalysis/testFoFiMembraneAnalysis') | ...
        HasName('testIGAMembraneAnalysis/testDDMFourPointSailAnalysis'));
end
resultIGAMembrane = run(suiteClassIGAMembrane);

%% Run the unit test cases for the isogeometric Kichhoff-Love shell analysis
warning('Test testLinearKirchoffLoveShellMultipatchAnalysis needs to be fixed');
suiteClassIGAKLShell = TestSuite.fromClass(?testIGAKirchhoffLoveShellAnalysis);
if isLight
    suiteClassIGAKLShell = suiteClassIGAKLShell.selectIf...
        (HasName('testIGAKirchhoffLoveShellAnalysis/testLinearKirchoffLoveShellAnalysis') | ...
        HasName('testIGAKirchhoffLoveShellAnalysis/testNonlinearKirchoffLoveShellAnalysis') | ...
        HasName('testIGAKirchhoffLoveShellAnalysis/testLinearKirchoffLoveShellMultipatchAnalysis'));
end
resultIGAKLShell = run(suiteClassIGAKLShell);

%% Run the unit test cases for the finite element formulation for plate in membrane action analysis
suiteClassFEMPlateInMembraneAction = TestSuite.fromClass(?testFEMPlateInMembraneActionAnalysis);
if isLight
    suiteClassFEMPlateInMembraneAction = suiteClassFEMPlateInMembraneAction.selectIf...
        (HasName('testFEMPlateInMembraneActionAnalysis/testCurvedPlateInMambraneActionSteadyStateLinear') | ...
        HasName('testFEMPlateInMembraneActionAnalysis/testCantileverBeamTransient'));
end
resultFEMPlateInMembraneAction = run(suiteClassFEMPlateInMembraneAction);

%% Run the unit test cases for stabilized finite element formulation for the Navier-Stokes problem
suiteClassFEM4CFD = TestSuite.fromClass(?testFEMComputationalFluidDynamicsAnalysis);
if isLight
    suiteClassFEM4CFD = suiteClassFEM4CFD.selectIf...
        (HasName('testFEMComputationalFluidDynamicsAnalysis/testFEM4NavierStokesSteadyState2D') | ...
        HasName('testFEMComputationalFluidDynamicsAnalysis/testFEM4NavierStokesTransientBossak3D') | ...
        HasName('testFEMComputationalFluidDynamicsAnalysis/testFEM4NavierStokesSteadyState2DReferencePaper'));
end
resultFEM4CFD = run(suiteClassFEM4CFD);

%% END OF SCRIPT
