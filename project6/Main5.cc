#include "mpi.h"
#include "Definitions.h"
#include "MeshUtilities.h"
#include "Quadrature.h"
#include "PostProcessorVtk.h"
#include "ElementTests.h"

#include "MaterialModelBar1D.h"

#include "TwoNodeBar.h" 
#include "Assembler.h"

const unsigned int         SpatialDimension = 3;
const unsigned int         DegreesOfFreedom = 3;
const unsigned int numberOfQuadraturePoints = 1;

typedef MaterialModels::MaterialModel1DBar                          MaterialModel;
typedef Elements::FiniteBar<MaterialModel,SpatialDimension>         Element;
typedef Elements::Properties                                        ElementProperties;
typedef Element::Node                                               Node;
typedef Element::Vector                                             Vector;
typedef Element::Point                                              Point;
typedef Element::Stress                                             Stress;
typedef Element::Strain                                             Strain;

typedef SingleElementMesh<Element>                                  Mesh;

typedef Element                                                     PhysicalElement;
typedef Assembler<PhysicalElement>                                  PhysicalAssembler;


const unsigned int NumberOfNodesPerElement = Element::NumberOfNodes;


int main(int arc, char *argv[]) {
  
  ignoreUnusedVariables(arc,argv);
  
  // The following lines simply create an output directory
  char sprintfBuffer[500];
  sprintf(sprintfBuffer,"Output_Main5");
  const string outputPath = string(sprintfBuffer);
  printf("Writing files to %s\n", outputPath.c_str());
  const bool createNewDirectories = true;
  Utilities::directoryCreator(outputPath, createNewDirectories, Quiet);
  
  // %%%%%%%%%%%%%%%%%%                                           %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%% Problem 4) (i) Creation of material model %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%                                           %%%%%%%%%%%%%%%%%%
  
  // TODO: Create your materialModel
  
  // ...
  
  // %%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%% Problem 4) (ii) Creation of mesh           %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%%%%%%%%
  
  Mesh mesh;

  const string meshFileName = "crossUnitCube.dat";
  MeshUtilities::readMeshFromFile<Element>(meshFileName,&mesh);
    
  array<size_t, SpatialDimension> numberOfCubesPerSide = {{6-1,6-1,6-1}};
  const double preperiodicSpatialTolerance = 1e-4;
  const double sideLength = 1.0;
  
  Vector unitXVector = Vector::Zero(); unitXVector(0) = sideLength;
  Vector unitYVector = Vector::Zero(); unitYVector(1) = sideLength;
  Vector unitZVector = Vector::Zero(); unitZVector(2) = sideLength;
  
  const array<Vector,SpatialDimension> patternVectors
    = {{unitXVector,unitYVector,unitZVector}};

  MeshUtilities::buildPeriodicMeshFromUnitCell(patternVectors             ,
                                               numberOfCubesPerSide       ,
                                               &mesh                      ,
                                               preperiodicSpatialTolerance);

  size_t numberOfNodes    = mesh._nodes.size();
  size_t numberOfElements = mesh._connectivity.size();

  cout << "Number of nodes in the mesh: "               << numberOfNodes    << endl;
  cout << "Number of beam elements in the mesh: "       << numberOfElements << endl;

  
  //%%%%%%%%%%%%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (iii) Preliminary stuff for elements %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%%%%
  
  // Create the element type and the element properties
  const double barDensity = 1522;
  const double barArea    = 1.0 ;
  ElementProperties elementProperties(barArea,barDensity);
  
  //%%%%%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (iv) Creation of the elements %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%%%%
  
  // TODO: Collect all physical elements
  vector<PhysicalElement> physicalElements; physicalElements.clear();
  for (unsigned int indexElement = 0; indexElement < (mesh._connectivity).size(); indexElement++){

    // ...
  
	}
  
  
  
  //%%%%%%%%%%%%%%%%%%%%                                         %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (v) Creation of an assembler %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                         %%%%%%%%%%%%%%%%%%%%
 
  // TODO: Create an assembler corresponding to your physical elements  
  
  // ...
 
  //%%%%%%%%%%%%%%%%%%%%                                         %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (vi) Eigenmodes and -freqs   %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                         %%%%%%%%%%%%%%%%%%%%
  vector<Vector> currentNodalDisplacement(numberOfNodes,Vector::Zero());
  
  // TODO: Evaluate stiffnessMatrix, consistentMassMatrix
  MatrixXd stiffnessMatrix;       // ...
  MatrixXd consistentMassMatrix;  // ...
  
  // TODO: Properly define the matrix of which we seek the eigenvalues or eigenvectors.
  //       If you prefer your own way of finding evals and evecs - go for it! :)
  MatrixXd matrixOfWhichWeSeekEigenvalues
    = MatrixXd::Identity(numberOfNodes*DegreesOfFreedom,numberOfNodes*DegreesOfFreedom); // ...
 
  Eigen::EigenSolver<MatrixXd> eigenEvaluator;
  eigenEvaluator.compute(matrixOfWhichWeSeekEigenvalues,true);
  
  Matrix<std::complex<double>,-1, 1> eigenValues  = eigenEvaluator.eigenvalues();
  Matrix<std::complex<double>,-1,-1> eigenVectors = eigenEvaluator.eigenvectors();
  
  
  //%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (vii) Output %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%
  
  // Empty essential BCs
  vector<EssentialBoundaryCondition> emptyEssentialBoundaryConditions;
  emptyEssentialBoundaryConditions.clear();
  
  // TODO: Show that we in fact recover six zero-energy modes and plot the first three eigenmodes
  // NOTE: If you don't like the for-loop we gave you, feel free to change it to whichever form
  //       you prefer. It's simply a guide and by no means >the< only solution.
  // HINT: this is just to show how you can grab evecs
  Matrix<std::complex<double>,-1,1> someEigenMode = eigenVectors.col(0); 
  // HINT: this is simply to illustrate how you can access their real part
  double someEntryOfTheEigenMode = someEigenMode(0,0).real(); 
  ignoreUnusedVariables(someEntryOfTheEigenMode);
  
  for (size_t indexMode = 0; indexMode < 3; indexMode++){
    
    // ...
    
    // Paraview Output
    
    // TODO: Define elementStresses
    const vector<Stress> elementStresses(numberOfElements,Stress::Zero()); // ...
    
    PostProcessors::Vtk::NamedArray<double> vtkStresses;
    vtkStresses._title="Bar Stresses";
    vtkStresses._elementWiseOrNodeWise = PostProcessors::Vtk::ElementWise;

    for (unsigned int elementIndex =0 ; elementIndex <numberOfElements; elementIndex++ )
    { 
      vtkStresses._array.push_back(elementStresses[elementIndex](0));
    }

    PostProcessors::Vtk::NamedArrays<int,double> vtkNamedArrays;
    vtkNamedArrays.addArray(vtkStresses);

    
    char outputFileDesignation[500];                                                     
    sprintf(outputFileDesignation,"%s/Eigenmode%03zu",outputPath.c_str(),indexMode);
    
    PostProcessors::Vtk::makeDeformedMeshFile<Element>( mesh                            ,
                                                        currentNodalDisplacement   ,
                                                        emptyEssentialBoundaryConditions,
                                                        string(outputFileDesignation)   ,
                                                        vtkNamedArrays                  );

  }
  
  return 0;
}

