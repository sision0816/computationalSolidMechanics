#include "Definitions.h"
#include "PostProcessorVtk.h"


int main() {

  // Geometrical properties
  // TODO: Set positive integer for numberOfElementsX and numberOfElementsY
  size_t numberOfElementsX = 8;
  size_t numberOfElementsY = 16;

  cout << "Let's build a mesh with numberOfElementsX=" << numberOfElementsX << " and numberOfElementsY=" << numberOfElementsY << endl;
  cout << "This should give a grant total of " << numberOfElementsX*numberOfElementsY*2 << " elements" << endl;

  // TODO: Set positive side lengths
  double sideLengthX = 0.4;
  double sideLengthY = 1.2;

  cout << "The side lengths are: sideLengthX=" << sideLengthX << ", sideLengthY=" << sideLengthY << endl;

  // Initialization of mesh
  Mesh rectangularTriangleMesh;

  // Incremental side length of one element
  // TODO: based on numberOfElementsX,numberOfElementsY
  //       and sideLengthX,sideLengthY, define dx,dy
  double dx = sideLengthX/double(numberOfElementsX);
  double dy = sideLengthY/double(numberOfElementsY);

  cout << "The incremental side lengths are: dx=" << dx << ", dy=" << dy << endl;

  // TODO: Add all nodes into the public member _nodes of rectangularTriangleMesh.
  //       You may want to use a for-loop or even a nested-for-loop to achieve this.
  size_t    currentNodeId       = 0;
  Vector2d  currentNodeLocation = Vector2d::Zero();

  for (size_t j=0; j<=numberOfElementsY; j++) {
    for (size_t i=0; i<=numberOfElementsX; i++) {

      currentNodeLocation(0) = i * dx;
      currentNodeLocation(1) = j * dy;

      rectangularTriangleMesh._nodes.push_back(Node(currentNodeId, currentNodeLocation));
      ++currentNodeId;
    }
  }

 
  // TODO: The goal of the next section is to create a triangular element as shown in the
  //       assignment sheet and to add it onto the _connectivity public member variable of
  //       rectangularTriangleMesh.
  //       In order to build the entire mesh, this then has to be repeated for all elements, the total  
  //       number of which is determined via numberOfElementsX and numberOfElementsY. You may want to use a 
  //       nested for-loop
  for (size_t j=0; j<numberOfElementsY; j++) {
    for (size_t i=0; i<numberOfElementsX; i++) {

      // Every triangle has 3 nodes
      array<int, 3> connection;

      // NodeID of bottom left corner
      size_t bottomLeftCorner = j*(numberOfElementsX+1)+i;

      // Triangle with 90 DEG corner at bottom right
      connection[0] = bottomLeftCorner;
      connection[1] = bottomLeftCorner+1;
      connection[2] = bottomLeftCorner+(numberOfElementsX+1)+1;
      rectangularTriangleMesh._connectivity.push_back(connection);

      // Triangle with 90 DEG corner at top left
      connection[0] = bottomLeftCorner+(numberOfElementsX+1)+1;
      connection[1] = bottomLeftCorner+(numberOfElementsX+1);
      connection[2] = bottomLeftCorner;
      rectangularTriangleMesh._connectivity.push_back(connection);
    }
  }
  
  // Make VTK file
  PostProcessors::Vtk::makeUndeformedMeshFile(rectangularTriangleMesh,string("Mesh"));

  cout << "Succesfully ran main and created .vtu file" << endl;

  return 0;
}



