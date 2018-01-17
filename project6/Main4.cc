#include "mpi.h"
#include "Definitions.h"
#include "MeshUtilities.h"
#include "Quadrature.h"
#include "PostProcessorVtk.h"  
#include "ElementTests.h"

#include "MaterialModelNeoHookean.h"
#include "ElementTypes.h"
#include "Wall.h"

#include "IsoparametricElement.h" 
#include "Assembler.h"
#include "SolverImplicit.h"

const unsigned int         SpatialDimension = 3;
const unsigned int         DegreesOfFreedom = 3;
const unsigned int numberOfQuadraturePoints = 1;

typedef MaterialModels::NeoHookean<SpatialDimension>                MaterialModel;

typedef ElementTypes::Simplex<SpatialDimension>                     ElementType;
typedef Elements::IsoparametricElement<MaterialModel, 
                                       numberOfQuadraturePoints,  
                                       ElementType, 
                                       DegreesOfFreedom>            Element;
typedef Element::Properties                                         ElementProperties;
typedef Element::Node                                               Node;
typedef Element::Vector                                             Vector;
typedef Element::Point                                              Point;
typedef Element::Stress                                             Stress;
typedef Element::Strain                                             Strain;

typedef SingleElementMesh<Element>                                  Mesh;

typedef Element                                                     PhysicalElement;
typedef Elements::ExternalForce::Wall<SpatialDimension,DegreesOfFreedom>
                                                                    ExternalElement;
typedef Assembler<PhysicalElement>                                  PhysicalAssembler;
typedef Assembler<ExternalElement>                                  ExternalAssembler;
typedef SolverImplicitDynamics<PhysicalAssembler,ExternalAssembler
                                                ,PhysicalAssembler> Solver;


const unsigned int NumberOfNodesPerElement = Element::NumberOfNodes;


int main(int arc, char *argv[]) {
  
  ignoreUnusedVariables(arc,argv);
  
    // The following lines simply create an output directory
  char sprintfBuffer[500];
  sprintf(sprintfBuffer,"Output_Main4");
  const string outputPath = string(sprintfBuffer);
  
  printf("Writing files to %s\n", outputPath.c_str());
  const bool createNewDirectories = true;
  Utilities::directoryCreator(outputPath, createNewDirectories, Quiet);
  
  // %%%%%%%%%%%%%%%%%%                                           %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%% Problem 4) (i) Creation of material model %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%                                           %%%%%%%%%%%%%%%%%%
  
  // TODO: Create your materialModel
  
  // ...
  const double materialModelParameter0 = 5.0 * power(10,5);
  const double materialModelParameter1 = 5.0 * power(10,7);
  
  MaterialModel materialModel(materialModelParameter0,materialModelParameter1);
  
  // %%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%% Problem 4) (ii) Creation of mesh           %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%%%%%%%%
  
  Mesh mesh;
  MeshUtilities::readMeshFromFile("SphereTetMesh_002000.dat",&mesh);
  
  const size_t numberOfNodes    = mesh._nodes.size();
  const size_t numberOfElements = mesh._connectivity.size();
  
  
  //%%%%%%%%%%%%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (iii) Preliminary stuff for elements %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%%%%
  
  // Create the element type and the element properties
  const double elementMultiplier = 1522;
  ElementType elementType;
  ElementProperties elementProperties(elementMultiplier);
  
  // Choose a quadrature rule
  const QuadratureRule<Element::SpatialDimension, numberOfQuadraturePoints>
  quadratureRule = Quadrature::buildSimplicialQuadrature<Element::SpatialDimension, numberOfQuadraturePoints>();
  
  ignoreUnusedVariables(elementType);
  ignoreUnusedVariables(quadratureRule);
  
  //%%%%%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (iv) Creation of the elements %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%%%%
  
  // Wall parameter
  const double wallStrength = 1.0*1.0e8;
  Vector wallOriginPosition = Vector::Zero(); wallOriginPosition (2)= +0.0;
  Vector wallNormalDirection= Vector::Zero(); wallNormalDirection(2)= -1.0;
  
  ignoreUnusedVariables(wallStrength);
  
  // TODO: Collect all external elements
  vector<ExternalElement> externalElements; externalElements.clear();
  for (unsigned int indexNode = 0; indexNode < numberOfElements /* TODO: set */; indexNode++){
    // ...
	 // Read out the nodes corresponding no the indexElement'th element
		array<Node,SpatialDimension+1> nodesSimplex;
		for (unsigned int indexNode = 0; indexNode < SpatialDimension+1; indexNode++)
    {
			nodesSimplex[indexNode] = mesh._nodes[mesh._connectivity[indexElement][indexNode]];
		}
    Element simplexElement(nodesSimplex,
                           elementProperties,
                           elementType,
                         & quadratureRule,
                         & materialModel);
    // REMINDER: You can push new elements into a vector via the .push_back option
    externalElements.push_back(simplexElement);
  }
  
  
  // TODO:Collect all physical elements
  vector<PhysicalElement> physicalElements; physicalElements.clear();
  for (unsigned int indexElement = 0; indexElement < numberOfElements /* TODO: set */; indexElement++){
      // ...
	  	 // Read out the nodes corresponding no the indexElement'th element
		array<Node,SpatialDimension+1> nodesSimplex;
		for (unsigned int indexNode = 0; indexNode < SpatialDimension+1; indexNode++)
    {
			nodesSimplex[indexNode] = mesh._nodes[mesh._connectivity[indexElement][indexNode]];
		}
    Element simplexElement(nodesSimplex,
                           elementProperties,
                           elementType,
                         & quadratureRule,
                         & materialModel);
    // REMINDER: You can push new elements into a vector via the .push_back option
    physicalElements.push_back(simplexElement);
	  
	}
  
  
  //%%%%%%%%%%%%%%%%%%%%                                         %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (v) Creation of an assembler %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                         %%%%%%%%%%%%%%%%%%%%
 
  // TODO: Create assemblers corresponding to your physical and external elements
 
  // ...
    PhysicalAssembler physicalAssembler(physicalElements, numberOfNodes);
	ExternalAssembler externalAssembler(externalElements, numberOfNodes);
  
  //%%%%%%%%%%%%%%%%%%%%                                    %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (vi) Creation of solver %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                    %%%%%%%%%%%%%%%%%%%%
 
  // TODO: Create an object of your SolverImplicitDynamics class
 
  // ...
  
  Solver solver(physicalAssembler,externalAssembler,physicalAssembler);
    
  const unsigned int  maxIterations     = 1000 ;
  const double        tolerance         = 1e-4  ;
  
  
 
  //%%%%%%%%%%%%%%%%%%%%                                                        %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (vii) Initialisation %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                        %%%%%%%%%%%%%%%%%%%%
  
  // TODO: Initiate all states that you need for your solver
  vector<Vector> currentNodalDisplacement(numberOfNodes,Vector::Zero());
  // ...
  vector<Vector> currentNodalAcceleration(numberOfNodes,Vector::Zero());
  vector<Vector> currentNodalVelocity(numberOfNodes,Vector::Vector::Zero());
  // TODO: Impose the initial velocity
  // ...
  for (unsigned int nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++) {
      for (unsigned int dofIndex = 0; dofIndex < DegreesOfFreedom; dofIndex++) {
        if (dofIndex == 0){
			currentNodalVelocity[nodeIndex](dofIndex) = 0;
		}
		else if (dofIndex == 1){
			currentNodalVelocity[nodeIndex](dofIndex) = 0;
		}
		else if (dofIndex == 2){
			currentNodalVelocity[nodeIndex](dofIndex) = -10;
		}
      }
    }
  
  // Empty boundary conditions
  vector<EssentialBoundaryCondition> emptyEssentialBoundaryConditions;
  emptyEssentialBoundaryConditions.clear();
  
  //%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 4) (viii) Run %%%%%%%%%%%%%%%%%%%% 
  //%%%%%%%%%%%%%%%%%%%%                       %%%%%%%%%%%%%%%%%%%%
  
  // TODO: Chose the number of loadsteps
  const unsigned int numberOfLoadsteps = 100; // ...
  
  for (unsigned int loadstepIndex = 0; loadstepIndex < numberOfLoadsteps; loadstepIndex++){
    
    if (loadstepIndex % unsigned(1) == 0) {
      printf("\ntimestep %6u (%%%5.1f) at %s\n",
               loadstepIndex,
               100. * loadstepIndex / float(numberOfLoadsteps),
               Utilities::getLocalTimeString().c_str());
    }
    
    // TODO: Call solver
    // ...
    
    vector<Vector> displacements
       = solver.computeNewmarkUpdate(essentialBCs, currentNodalDisplacement,currentNodalVelocity,currentNodalAcceleration, maxIterations, tolerance, true);

    if (!(loadstepIndex%1)){
      
      printf("Giving output at loadstep (%d/%d).\n",loadstepIndex,numberOfLoadsteps);
      
      // TODO: define the following four vectors
      const vector<Stress> nodeStresses    (numberOfNodes   ,Stress::Zero()); // ...
      const vector<Stress> elementStresses (numberOfElements,Stress::Zero()); // ...
      const vector<Strain> nodeStrains     (numberOfNodes   ,Strain::Zero()); // ...
      const vector<Strain> elementStrains  (numberOfElements,Stress::Zero()); // ...
        nodeStresses = assembler.computeElementStresses(displacements);
        elementStresses = assembler.computeNodalStresses (displacements) ;
        nodeStrains = assembler.computeElementStrains(displacements);
        elementStrains= assembler.computeNodalStrains(displacements);
      

      // !!! NOTHING TO BE DONE FROM HERE ONWARDS !!!                                                   
      sprintf(sprintfBuffer,"%s/WallImpact_%03u",outputPath.c_str(),loadstepIndex);
      
      PostProcessors::Vtk::NamedArrays<int,double> vtkNamedArrays;
      PostProcessors::Vtk::makeDeformedMeshFile<Element>( mesh,
                                                          currentNodalDisplacement,
                                                          nodeStresses,
                                                          elementStresses,
                                                          MaterialModel::getStressComponentNames(),
                                                          nodeStrains,
                                                          elementStrains,
                                                          MaterialModel::getStrainComponentNames() ,
                                                          emptyEssentialBoundaryConditions,
                                                          string(sprintfBuffer),
                                                          vtkNamedArrays);
    }
    
  }
  
  return 0;
}
