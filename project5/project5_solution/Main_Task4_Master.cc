#include "Definitions.h"
#include "MeshUtilities.h"
#include "PostProcessorVtk.h"
#include "ElementTests.h"

#include "MaterialModelBar1D.h"
#include "TwoNodeBar.h"
#include "ExternalForces.h"
#include "Assembler.h"
#include "Solver.h"

#define SPATIAL_DIMENSION 2
const unsigned int SpatialDimension = SPATIAL_DIMENSION;
const unsigned int DegreesOfFreedom = SpatialDimension;

typedef MaterialModels::MaterialModel1DBar                  MaterialModel;
typedef Elements::FiniteBar<MaterialModel,SpatialDimension> Element;
typedef Elements::Properties                                ElementProperties;
typedef Element::Point                                      Point;
typedef Element::Vector                                     Vector;
typedef Element::Node                                       Node;
typedef Element::Stress                                     Stress;
typedef Element::Strain                                     Strain;

typedef Element                                             PhysicalElement;
typedef Elements::ExternalForce::ConstantBodyForce<Element> ExternalElement;

typedef SingleElementMesh<Element>                          Mesh;

int main(int arc, char *argv[]) {
  
  ignoreUnusedVariables(arc,argv);
  
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%     1 : PRELIMINARY     %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
 
 
  // Material & Geometric parameters
  const double youngModulus               = 0.1;     // E = 0.1   GPa 
  const double beamArea                   = 0.01;    // A = 0.01  m^2
  const double beamDensity                = 1.0e-4;  // rho = 1.0e5 kgm^-3
  const double gravitationalAcceleration  = 9.81;    // g = 9.81 ms^-2
  
  // Delete this later:
  ignoreUnusedVariables(youngModulus);
  
  Vector bodyForceVector  = Vector::Zero(); 
  bodyForceVector(1)      = -gravitationalAcceleration;

  // Initialization of elemenent properties
  ElementProperties elementProperties(beamArea,beamDensity);
  
  // Initialize the material model (We use the linear elastic one)
  MaterialModel materialModel(youngModulus);

  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%     2 : MESH            %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  Mesh mesh; // This is no longer a finite element mesh but a construction of several beam elements

  const string meshFileName = "meshForProblem4.dat";
  MeshUtilities::readMeshFromFile<Element>(meshFileName,&mesh); // !!! EXCHANGE WITH PHYISCAL ELEMENT FOR CONSISTENCY
    
  size_t numberOfNodes    = mesh._nodes.size();
  size_t numberOfElements = mesh._connectivity.size();
     
  Point meshMinPosition;
  Point meshMaxPosition;
  MeshUtilities::findBoundingBoxOfGeneralMesh<Mesh,SpatialDimension>(mesh            ,
                                                                     &meshMinPosition,
                                                                     &meshMaxPosition);

  cout << "Number of nodes in the mesh: "         << numberOfNodes    << endl;
  cout << "Number of beam elements in the mesh: " << numberOfElements << endl;
  cout << "Min position: " << meshMinPosition.transpose() << endl;
  cout << "Max position: " << meshMaxPosition.transpose() << endl;
  
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%     3 : ASSEMBLER       %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  vector<PhysicalElement> physicalElementVector;
  vector<ExternalElement> externalElementVector;
   
  // The assembler needs to know the physical element as well as the external element
  // 1) Loop through all beam-elements of the mesh
  // 2) Find the nodes of each beam-element (Reminder: mesh has a private member called
  //    _connectivity ...)
  // 3) Create the PhysicalElement and push it back into the physicalElementVector
  //    - if you're unsure about how to create it ... check out the constructor in
  //      TwoNodeBar.h
  // 4) Create the externalElement and push it back into the externalElementVector
  //    - if you're unsure about how to create it ... check out the constructor in
  //      ExternalForces.h

  for (size_t elementIndex = 0; elementIndex < numberOfElements; ++elementIndex) // (Step 1 - given)
  {
    const size_t node1 = mesh._connectivity[elementIndex][0];
    const size_t node2 = mesh._connectivity[elementIndex][1];

    const array<Node, PhysicalElement::NumberOfNodes> elementNodes = {{mesh._nodes[node1], mesh._nodes[node2]}};

    physicalElementVector.push_back(PhysicalElement(  elementNodes     ,
                                                      elementProperties,
                                                      &materialModel   ));
           
    externalElementVector.push_back(ExternalElement(physicalElementVector[elementIndex],
                                                    bodyForceVector                    ));
  }

  Assembler<PhysicalElement> physicalAssembler(physicalElementVector,numberOfNodes);
  Assembler<ExternalElement> externalAssembler(externalElementVector,numberOfNodes);
  
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%     4 : ESS. BCs        %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%

  
  // Check through every node -> fix bottom row
  const double nodePositionSearchTolerance = 1e-6;
  vector<EssentialBoundaryCondition> essentialBCs;
  
  // 1) Loop through all nodes
  // 2) Find all nodes that are located on the bottom and constrain them to zero displacement
  // 3) Use "meshMinPosition" to find the respective nodes. 
  // 4) The following example shows the essential boundary condition implementation  
  //    for node 5, with a deformation of 3.0 in its x direction:
  //    essentialBCs.push_back(EssentialBoundaryCondition(5, 0, 3.0));
  // 5) Push back the essential boundary conditions into the vector essentialBCs
  
  for (unsigned int nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    
    // Find position of current node
    const Point nodePosition = mesh._nodes[nodeIndex]._position;
  
    // Check if within "bottom-row"
    if (fabs(nodePosition(1)-meshMinPosition(1))<nodePositionSearchTolerance){
      // If so, constrain all 
      for (unsigned int indexDOF=0; indexDOF< DegreesOfFreedom; indexDOF++){
        essentialBCs.push_back(EssentialBoundaryCondition(nodeIndex, indexDOF, 0.));
      }
    }
    
  }
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%     5 : SOLVER          %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  Solver<Assembler<PhysicalElement>,Assembler<ExternalElement>> solver(physicalAssembler,externalAssembler);
    
  const unsigned int  maxIterations     = 10000 ;
  const double        tolerance         = 1e-8  ;
  
  vector<Vector> initialGuess(numberOfNodes,Vector::Zero());
  
  vector<Vector> displacements
    = solver.computeNewtonRaphsonSolution(essentialBCs, initialGuess, maxIterations, tolerance, true);
    
    
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%     6 : PARAVIEW        %%%%%%%%%%%%%%%%%%%%%%%%%
  // %%%%%%%%%%%%%%%%%%%%%%%%%%                         %%%%%%%%%%%%%%%%%%%%%%%%%

  // Access the assembler to compute the element stresses which can then be 
  // used as an additional output in our vtk file:
  const vector<Stress> elementStresses =
		  physicalAssembler.computeElementStresses(displacements);

  PostProcessors::Vtk::NamedArray<double> vtkStresses;
  vtkStresses._title="Bar Stresses";
  vtkStresses._elementWiseOrNodeWise = PostProcessors::Vtk::ElementWise;

  for (unsigned int elementIndex =0 ; elementIndex <numberOfElements; elementIndex++ )
  { 
	  vtkStresses._array.push_back(elementStresses[elementIndex](0));
  }

  PostProcessors::Vtk::NamedArrays<int,double> vtkNamedArrays;
  vtkNamedArrays.addArray(vtkStresses);

  // This function will make a file called TwoNodeBarMesh.vtu that can be
  // displayed in paraview
  PostProcessors::Vtk::makeDeformedMeshFile<Element>( mesh                    ,
                                                      displacements           ,
                                                      essentialBCs            ,
                                                      string("TwoNodeBarMesh"),
                                                      vtkNamedArrays          );
                                                      
  printf("Generating Output Paraview File \n");
  return 0;
}
