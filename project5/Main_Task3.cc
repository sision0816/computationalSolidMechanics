// -*- C++ -*-
#include "MeshUtilities.h"
#include "Quadrature.h"
#include "PostProcessorVtk.h"
#include "ElementTests.h"
#include "Assembler.h"
#include "Solver.h"
#include "MaterialModelLinearElastic.h"
#include "MaterialModelNeoHookean.h"
#include "ElementTypes.h"
#include "IsoparametricElement.h" 
#include "mpi.h"
#include "Definitions.h"

const unsigned int         SpatialDimension = 2;
const unsigned int         DegreesOfFreedom = 2;
const unsigned int numberOfQuadraturePoints = 1;

//typedef MaterialModels::MaterialModelLinearElastic<SpatialDimension>  MaterialModel;
typedef MaterialModels::NeoHookean<SpatialDimension>                  MaterialModel;
typedef ElementTypes::Simplex<SpatialDimension>                       ElementType;
typedef Elements::IsoparametricElement<MaterialModel, 
                                       numberOfQuadraturePoints,  
                                       ElementType, 
                                       DegreesOfFreedom>              Element;
typedef Element::Properties                                           ElementProperties;
typedef Element::Node                                                 Node;
typedef Element::Vector                                               Vector;
typedef Element::Point                                                Point;
typedef Element::Stress                                               Stress;
typedef Element::Strain                                               Strain;

typedef SingleElementMesh<Element>                                    Mesh;

const unsigned int                        NumberOfNodesPerElement =   Element::NumberOfNodes;


int main(int arc, char *argv[]) {
  
  ignoreUnusedVariables(arc,argv);

  
  
  
  // %%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%% Problem 3) (ii) Creation of material model %%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%                                            %%%%%%%%%%%%%%%%%%
  
  // ToDo:
  // Implement a linear-elastic material model and the neo-hookian material model the same 
  // way as in project 2. The material constants are given in the project exercise.

  const double materialModelParameter0 = 1.0;
  const double materialModelParameter1 = 4.0;
  
  MaterialModel materialModel(materialModelParameter0,materialModelParameter1);
  

  
  
  // %%%%%%%%%%%%%%%%%%                                       %%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%% Problem 3) (i) Initialization of mesh %%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%                                       %%%%%%%%%%%%%%%%%%%%

  // ToDo: 
  // Create a finite element mesh of the format used in HW1. Don't worry, you don't
  // have to type out everything again. We provide you with a utility function. You
  // simply need to set the right lengths, etc.
  // Check the exercise sheet to find information about the size of the mesh. 
 
  Mesh mesh;

  array<double,2> sideLengths     ;
  sideLengths[0]=1.0; sideLengths[1]=0.4;
  
  array<size_t,2> numberOfElementsPerSide;
  numberOfElementsPerSide[0] = 20; numberOfElementsPerSide[1] = 8;
  
  MeshUtilities::buildRectangularTriangleMesh(sideLengths             ,
                                              numberOfElementsPerSide ,
                                              &mesh                   );
  
  const size_t numberOfNodes    = mesh._nodes.size();
  const size_t numberOfElements = mesh._connectivity.size();
  std::cout << "Number of Nodes in the Mesh: "    << numberOfNodes    << endl;
  std::cout << "Number of Elements in the Mesh: " << numberOfElements << endl;
  
  // Find maximum and minimum position of the mesh. This will be helpful for the 
  // essential boundary conditions later. (Of course, you don't need to use it,
  // if you have an alternative way. It may just be useful, that's all we're saying :))

  Point meshMinPosition;
  Point meshMaxPosition;
  MeshUtilities::findBoundingBoxOfGeneralMesh<Mesh,Element::SpatialDimension>
    (mesh            ,
     &meshMinPosition,
     &meshMaxPosition);
  

  
  
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 3) (ii) Preliminary stuff for elements %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
  
  // Create the element type and the element properties
  const double density = 1.;
  ElementType elementType;
  ElementProperties elementProperties(density);
  
  // Choose a quadrature rule
  const QuadratureRule<Element::SpatialDimension, numberOfQuadraturePoints>
  quadratureRule = Quadrature::buildSimplicialQuadrature<Element::SpatialDimension, numberOfQuadraturePoints>();
  
  
  
  
  //%%%%%%%%%%%%%%%%%%%%                                           %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 3) (iii) Creation of the elements %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                           %%%%%%%%%%%%%%%%%%%%
  
  // Create element
  vector<Element> elements;
  
  for (unsigned int indexElement = 0; indexElement < (mesh._connectivity).size(); indexElement++)
  {
    
    // Read out the nodes corresponding no the indexElement'th element
		array<Node,SpatialDimension+1> nodesSimplex;
		for (unsigned int indexNode = 0; indexNode < SpatialDimension+1; indexNode++)
    {
			nodesSimplex[indexNode] = mesh._nodes[mesh._connectivity[indexElement][indexNode]];
		}
		
    // TODO: Create a simplex element
    // REMINDER: a simplex element wants 1) an array of nodes, 2) element properties,
    //           3) element type, 4) address of quadrature rule, 5) address of
    //           material model for its compiler
    Element simplexElement(nodesSimplex,
                           elementProperties,
                           elementType,
                         & quadratureRule,
                         & materialModel);
    // ...
	
      
		// TODO: Push the newly created Element into elements.
    // REMINDER: You can push new elements into a vector via the .push_back option
    elements.push_back(simplexElement);
    // ...
    
	}

 
 
 
  //%%%%%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 3) (iv) Creation of an assembler %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                          %%%%%%%%%%%%%%%%%%%%
 
  // Create an assembler, which interfaces with all elements to assemble the global energy, nodal forces
  // and stffness matrix by requesting the individual element quantities from the respective elements and
  // assembling them into large (global) vectors and matrices;
  Assembler<Element> assembler(elements, numberOfNodes);
 
  // Delete this later:
  ignoreUnusedVariables(quadratureRule);
  ignoreUnusedVariables(elementType);
 
 
 

  //%%%%%%%%%%%%%%%%%%%%                                                        %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%% Problem 3) (v) Creation of essentialBoundaryConditions %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%                                                        %%%%%%%%%%%%%%%%%%%%


  const unsigned int  numberOfLoadsteps = 10.0; // TODO: Choose as many as you need for convergence
  const double        maxTension        = 0.2; // TODO: Choose as instructed on the project assignment
  vector<Vector>      initialGuess(numberOfNodes, Vector::Zero());
  
  double appliedDeformation             = 0.;
  const double increment = maxTension/(double)numberOfLoadsteps;
  

  for (unsigned int loadstepIndex = 0; loadstepIndex < numberOfLoadsteps; loadstepIndex++) {
  
    appliedDeformation += increment;
    cout << "starting load step " << loadstepIndex << endl;
  
    vector<EssentialBoundaryCondition> essentialBCs;
  
    // ToDo:
    // 1) Loop through all node indices.
    // 2) Find all nodes that are located on the left and fix them and move all 
    //    nodes located on the right side of the cantilever and move them by an incremental displacement. 
    // 3) Use "meshMinPosition" and "meshMaxPosition" to find the respective nodes. 
    // 4) The following example shows the essential boundary condition implementation 
    //    for node 5, with a deformation of 3.0 in its x direction:
    //    essentialBCs.push_back(EssentialBoundaryCondition(5, 0, 3.0));
    // 5) Push back the essential boundary conditions into the vector (using the functionality .push_back)

    // ...
    for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++)
    {
      Point nodePosition = mesh._nodes[nodeIndex]._position;

      if (nodePosition(0) == meshMinPosition(0))
      {
        essentialBCs.push_back(EssentialBoundaryCondition(nodeIndex, 0, 0.0));
        essentialBCs.push_back(EssentialBoundaryCondition(nodeIndex, 1, 0.0));
      }
      if (nodePosition(0) == meshMaxPosition(0))
      {
        essentialBCs.push_back(EssentialBoundaryCondition(nodeIndex, 0, appliedDeformation*sideLengths[0]));
        essentialBCs.push_back(EssentialBoundaryCondition(nodeIndex, 1, 0.0));
      }

    }
    
    // Delete this later:
    ignoreUnusedVariables(loadstepIndex);
  
  
  
  
    //%%%%%%%%%%%%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%% Problem 3) (vi) Creation of an iterative solver %%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%                                                 %%%%%%%%%%%%%%%%%%%%

  
    // ToDo: Use the Newton-Raphson solver that uses the assembler and the
    // essential BCs to find the sought displacement at all nodes).

    Solver<Assembler<Element>,Assembler<Element>> solver(assembler);
    
    const unsigned int  maxIterations     = 10000 ;
    const double        tolerance         = 1e-8  ;
    
    vector<Vector> displacements
      = solver.computeNewtonRaphsonSolution(essentialBCs, initialGuess, maxIterations, tolerance, true);

    // ToDo:
    // Update your initial guess for the next loadstepIndex
    initialGuess = displacements;
    // ...

    
    
    
    //%%%%%%%%%%%%%%%%%%%%                                                  %%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%% Problem 3) (vii) Create Paraview file for output %%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%                                                  %%%%%%%%%%%%%%%%%%%%

    PostProcessors::Vtk::NamedArrays<int,double> vtkNamedArrays;    

    // This function will make a number of files called IsoparametricElementsMesh_xyz.vtu that 
    // can be visualized in paraview
    
    // ToDo:
    // Access the assembler to compute the element/node stresses and element/node strains which can then be 
    // used as an additional output in our vtk file:
    // NOTE: The following functions exist in assembler:
    //       computeElementStresses, computeNodalStresses, computeElementStrains, computeNodalStrains
    const vector<Stress> elementStresses = assembler.computeElementStresses(displacements); // ...
    const vector<Stress> nodeStresses = assembler.computeNodalStresses (displacements)       ; // ...
    const vector<Strain> elementStrains = assembler.computeElementStrains(displacements) ; // ...
    const vector<Strain> nodeStrains = assembler.computeNodalStrains(displacements)       ; // ...

    // !!! NOTHING TO BE DONE FROM HERE ONWARDS !!!
    
    char sprintfBuffer[500];                                                     
    sprintf(sprintfBuffer,"IsoparametricElementsMesh_%03u", loadstepIndex);
    
    PostProcessors::Vtk::makeDeformedMeshFile<Element>( mesh,
                                                        displacements,
                                                        nodeStresses,
                                                        elementStresses,
                                                        MaterialModel::getStressComponentNames(),
                                                        nodeStrains,
                                                        elementStrains,
                                                        MaterialModel::getStrainComponentNames() ,
                                                        essentialBCs,
                                                        string(sprintfBuffer),
                                                        vtkNamedArrays);
  }
  
  return 0;
}
