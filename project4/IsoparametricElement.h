#ifndef SIMPLICIAL_ELEMENT
#define SIMPLICIAL_ELEMENT

#include "/src/Definitions.h"
#include "/src/Utilities.h"
#include "/src/FiniteDeformationUtilities.h"

namespace Elements {
  
  
template <class MaterialModelType, unsigned int QP, class ElementType, unsigned int DOFs,
          template<class> class NodeType = NodeWithId> // keep ElementType and QP general
class IsoparametricElement {

public:

  struct Properties {
    
    Properties(const double & elementMultiplier = 1.) :
      _elementMultiplier(elementMultiplier) {
    }
    
    const double _elementMultiplier;
  };


  static const unsigned int                    QuadPoints = QP;
  static const unsigned int                 NumberOfNodes = ElementType::NumberOfNodes;
  static const unsigned int              SpatialDimension = ElementType::SpatialDimension;
  static const VTKCellType                    VtkCellType = ElementType::VtkCellType;
  static const unsigned int      NumberOfNodesPerBoundary = ElementType::NumberOfNodesPerBoundary;
  static const unsigned int            NumberOfBoundaries = ElementType::NumberOfBoundaries;
  static const unsigned int              DegreesOfFreedom = DOFs;
  
  typedef Eigen::Matrix<double, DegreesOfFreedom, 1>        Vector;
  typedef typename ElementType::Point                       Point;
  typedef NodeType<Point>                                   Node;
  typedef array<Vector, NumberOfNodes>                      NodalDisplacements;
  typedef array<Vector, NumberOfNodes>                      Forces;
  typedef Matrix<double,
                 NumberOfNodes*DegreesOfFreedom,
                 NumberOfNodes*DegreesOfFreedom>            StiffnessMatrix;

  typedef MaterialModelType                                 MaterialModel;
  typedef typename MaterialModel::DisplacementGradient      DisplacementGradient;
  typedef typename MaterialModel::Strain                    Strain;
  typedef typename MaterialModel::Stress                    Stress;
  typedef typename MaterialModel::TangentMatrix             TangentMatrix;
  typedef Matrix<double,DegreesOfFreedom,SpatialDimension>  StressTensor;
  typedef SingleElementBoundary<Node, NumberOfNodesPerBoundary>    SingleBoundary;
  typedef AllElementBoundaries<SingleBoundary, NumberOfBoundaries> AllBoundaries;


  //TODO: create the constructor
  IsoparametricElement(const array<Node, NumberOfNodes> &           nodes         ,
                       const Properties &                           properties    ,
                       const ElementType &                          elementType   ,
                       const QuadratureRule<SpatialDimension, QP> * quadratureRule,
                       const MaterialModel *                        materialModel ) :
      _properties     (properties)    ,
      _quadratureRule (quadratureRule), 
      _materialModel  (materialModel) {

    ignoreUnusedVariables(nodes);
    ignoreUnusedVariables(elementType);
    // cout<<"the first quadrature point in the reference configuration"<< (*_quadraqureRule)._point[0]
    // TODO: We want to copy the information of nodes in terms of ID's and positions
    //       into _nodeIds and _nodePositions, respectively
    // NOTE: nodes contains the global id's of all NumberOfNodes nodes of THIS element
    // NOTE: i.e. NumberOfNodes is not the mesh's total number of nodes, but the
    //       number of nodes in the element
    for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; ++nodeIndex) {
      _nodeIds[nodeIndex]       = 0;
      _nodeIds[nodeIndex]       = nodes[nodeIndex]._id;
      _nodePositions[nodeIndex] = Point::Zero();
      _nodePositions[nodeIndex] = nodes[nodeIndex]._position;
    }
    
    
    
    //TODO: The aim is to take the shape functions and shapefunction-derivatives as
    //      provided by elementType and to appropriately store them inside
    //      _shapeFunctions and _shapeFunctionDerivatives. Using the jacobian and
    //      the quadrature weights, we are further going to set _weightedJacobian
    double volume = 0.;
    ignoreUnusedVariables(volume);
    
    for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
      
      // TODO: Read in the qpth quadrature point from _quadratureRule
      // NOTE: _quadratureRule is an address to a QuadratureRule, and the class
      //       QuadratureRule has the public member _points which is an array that
      //       stores all QPs' position
      //cout << "HINT: This is the position of the first quadraturePoint: " << endl;
      //cout << _quadratureRule->_points[0] << endl;
      
      const Point quadPoint;
      array <Point, QP> quadPoints;
      quadPoints[qpIndex] = _quadratureRule->_points[qpIndex];
      
      
      // Initialization of jacobian, incl. initial zero'ing
      Matrix<double, SpatialDimension, SpatialDimension> jacobian
        = Matrix<double, SpatialDimension, SpatialDimension>::Zero();
      
      // TODO: Read out the shapeFunctionDerivatives provided by elementType for this
      //       particular quadrature point.
      // NOTE: The reference element is defined in elementType, an object passed to_char_type
      //       this constructor
      array<Point, NumberOfNodes> shapeFunctionDerivatives;
      shapeFunctionDerivatives = elementType.computeShapeFunctionDerivatives(quadPoints[qpIndex]);
      
      // TODO: In a (nested) for-loop, evaluate the entries of the jacobian as dervived in
      //       lecture
      for (unsigned int nodeId = 0; nodeId < NumberOfNodes; ++nodeId) {
        for (unsigned int i = 0; i < SpatialDimension; ++i) {
          for (unsigned int j = 0; j < SpatialDimension; ++j) {
            jacobian(i, j) += _nodePositions[nodeId](j) * shapeFunctionDerivatives[nodeId](i);
          }
        }
      }
      
      // TODO: Define the qpIndex'th entry of _shapeFunctions and _shapeFunctionDerivatives
      // HINT: As you can see from the equations given in the assignment sheet,
      //       there's not really much to do for _shapeFunctions, while there is
      //       a tad bit more to do for _shapeFunctionDerivatives
      // HINT: _shapeFunctions and _shapeFunctionDerivatives are arrays containing arrays
      // HINT: someMatrix*someArrayofMatrices is an allowed operation thanks to the
      //       the std- and Eigen-lib.

      _shapeFunctions[qpIndex] = elementType.computeShapeFunctions(quadPoints[qpIndex]);

      for (unsigned int nodeId = 0; nodeId < NumberOfNodes; ++nodeId)
      {
        _shapeFunctionDerivatives[qpIndex][nodeId] = jacobian.inverse() * shapeFunctionDerivatives[nodeId];
      }
      
      // TODO: Set _weightedJacobian[qpIndex]. This should basically equate to J * w_k on
      //       the HW sheet - w_k are the quadrature weights, which for example for the
      //       constant element (i.e. one quadrature point) would just simply be '1'
      //cout << "HINT: This is the weight associated to the first quadrature point " << endl;
      //cout << _quadratureRule->_weights[0] << endl;
      
      _weightedJacobian[qpIndex] = jacobian.determinant() * _quadratureRule->_weights[qpIndex] ; // ...
      
      // TODO: (Accumulatively) update volume, describing the total "weighted Jacobian"
      volume += _weightedJacobian[qpIndex]; // ...
      
    }
    
    // TODO: Finally, based on the volume we evaluated in the above for-loop, we now
    //       evaluate _nodalWeights
    for (unsigned int nodeIndex = 0; nodeIndex < NumberOfNodes; nodeIndex++) {
      _nodalWeights[nodeIndex] = 1.11; // ...
    }
  }
  
  //TODO: compute displacement gradient u_{i,j} in Voigt notation at all quadrature points
  array<DisplacementGradient, QP>
  computeDispGradsAtGaussPoints(const NodalDisplacements & displacements) const {
    
    // TODO: Initialization of an array displacementGradients of displament
    //       gradients at all quadrature points
    array<DisplacementGradient, QP> displacementGradients;

      // TODO: Based on the discrete set of displacements at the nodes and using
    //       the shape function derivatives, evaluate the displacement gradients
    //       at the Gauss-points
    // NOTE: None of the entries of displacementGradients has been zeroed yet...
	for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
      for (unsigned int i = 0; i < DegreesOfFreedom; i++) {
        for (unsigned int j = 0; j < SpatialDimension; j++) {
            displacementGradients[qpIndex](i*SpatialDimension+j) = 0.0;        
        }
      }
    }
    for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
      for (unsigned int i = 0; i < DegreesOfFreedom; i++) {
        for (unsigned int j = 0; j < SpatialDimension; j++) {
          for (unsigned int nodeId = 0; nodeId < NumberOfNodes; nodeId++) {
            displacementGradients[qpIndex](i*SpatialDimension+j) += displacements[nodeId](i) * _shapeFunctionDerivatives[qpIndex][nodeId](j); // ... (a lot to be changed here!)
          }
        }
      }
    }
    
    // Return
    return displacementGradients;
  }
  
  //TODO: compute the potential energy of the element by summation over all quadrature points
  double
  computeEnergy(const NodalDisplacements &  displacements ) const {
                  
    ignoreUnusedVariables(time);
    ignoreUnusedVariables(displacements);
    
    // TODO: Using the computeDispGradsAtGaussPoints function you just defined above,
    //       obtain the displacement gradients at all quadrature points    
    const array<DisplacementGradient, QP> displacementGradients = computeDispGradsAtGaussPoints(displacements); // ...
    
    // TODO: Sweeping through all quadrature points, evaluate the total energy
    double energy = 0.;
    for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
      energy += _weightedJacobian[qpIndex] * _materialModel->computeEnergy(displacementGradients[qpIndex]); // ...
    }
    
    
    // Return - Don't forget to multiply by _properties._elementMultiplier
    return energy * _properties._elementMultiplier;
  }
  
  //TODO: compute the nodal forces of the element by summation over all quadrature points
  Forces
  computeForces(const NodalDisplacements &  displacements ) const {
    
    ignoreUnusedVariables(time);
    ignoreUnusedVariables(displacements);
    const array<DisplacementGradient, QP> displacementGradients = computeDispGradsAtGaussPoints(displacements);     
    // Initialization of forcesAtNodes
    Forces forcesAtNodes;
	for (size_t nodeId = 0; nodeId < NumberOfNodes; nodeId++){
		forcesAtNodes[nodeId] = Vector::Zero();
	}
    
    // TODO: Based on the displacementGradients at the Q.P., evaluate the forces at all nodes
    // NOTE: None of the entries of forcesAtNodes has been zeroed yet...
    for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
      Stress  stress = _materialModel->computeStress(displacementGradients[qpIndex]);
	  Matrix<double,SpatialDimension,SpatialDimension> matrixStress
      = Utilities::convertTensorFromVoigtToStandard<SpatialDimension>(stress);
      for (size_t nodeId = 0; nodeId < NumberOfNodes; nodeId++){
         
         forcesAtNodes[nodeId] +=  matrixStress * _weightedJacobian[qpIndex] *  _shapeFunctionDerivatives[qpIndex][nodeId];// ...
      }
    }
    
    // Return - don't forget to multiply by _properties._elementMultiplier
    return forcesAtNodes * _properties._elementMultiplier;
  }
  
  //TODO: compute the stiffness matrix of the element by summation over all quadrature points
  StiffnessMatrix
  computeStiffnessMatrix(const NodalDisplacements & displacements ) const {
    
    ignoreUnusedVariables(time);
    ignoreUnusedVariables(displacements);
    
    // Initialize zero-stiffness matrix
    Matrix<double, NumberOfNodes*DegreesOfFreedom, NumberOfNodes*DegreesOfFreedom> stiffness
      = Matrix<double, NumberOfNodes*DegreesOfFreedom, NumberOfNodes*DegreesOfFreedom>::Zero();

    // TODO: Using the computeDispGradsAtGaussPoints function you just defined above,
    //       obtain the displacement gradients at all quadrature points
    const array<DisplacementGradient, QP> displacementGradients = computeDispGradsAtGaussPoints(displacements);
    
    // TODO: Based on the displacementGradients at the Q.P., evaluate the several entries
    //       of the stiffness matrix by sweeping through all Q.P.s and evaluating the local
    //       incremental stiffness matrix
    // HINT: The function convertStandardIndicesToVoigtIndex in the Utilities namespace might
    //       proof useful here (you can find it defined in FiniteDeformationUtilities.h) in
    //       the /src/ directory

    for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
	  TangentMatrix tangentMatrix = _materialModel->computeTangentMatrix(displacementGradients[qpIndex]);
      array<array<array<array<double,SpatialDimension>,SpatialDimension>,SpatialDimension>,SpatialDimension> tangentMatrixAsArray = Utilities::convertFourthOrderTensorFromVoigtToStandard<SpatialDimension>(tangentMatrix);
      for (unsigned int nodeId_a = 0; nodeId_a < NumberOfNodes; nodeId_a++){
        for (unsigned int nodeId_b = 0; nodeId_b < NumberOfNodes; nodeId_b++){
          for (unsigned int i = 0; i < DegreesOfFreedom; i++){
            for (unsigned int n = 0; n < SpatialDimension; n++) {
              for (unsigned int j = 0; j < DegreesOfFreedom; j++){
                for (unsigned int l = 0; l < DegreesOfFreedom; l++){
                  stiffness(nodeId_a * DegreesOfFreedom + i, nodeId_b * SpatialDimension + n) += _weightedJacobian[qpIndex] * tangentMatrixAsArray[i][j][n][l] * _shapeFunctionDerivatives[qpIndex][nodeId_a][j] * _shapeFunctionDerivatives[qpIndex][nodeId_b][l];
                }
              }
            }
          }
        }
      }
    }
    //need to be check
    // Return - don't forget to multiply by _properties._elementMultiplier
    return stiffness * _properties._elementMultiplier;
  }
  
  array<size_t, NumberOfNodes>
  getNodeIds() const {
    return _nodeIds;
  }
  
  array<Vector, NumberOfNodes>
  getNodePositions() const {
    return _nodePositions;
  }

  array<double, QP>
  computeWeightsAtGaussPoints() const {
    return _weightedJacobian;
  }

  array<double, NumberOfNodes>
  computeNodalWeights() const {
    return _nodalWeights;
  }

  //TODO: compute stresses all quadrature points
  array<Stress, QP>
  computeStressesAtGaussPoints(const NodalDisplacements & displacements ,
                               const double               time          ) const {
    
    // TODO: Using the computeDispGradsAtGaussPoints function you just defined above,
    //       obtain the displacement gradients at all quadrature points
    const array<DisplacementGradient, QP> displacementGradients = computeDispGradsAtGaussPoints(displacements);
    
    // TODO: Sweep through all quadrature points, evaluate the stress tensor and save
    //       it in stresses
    array<Stress, QP> stresses;
    for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
      stresses[qpIndex] = Stress::Zero();
      stresses[qpIndex] = _materialModel->computeStress(displacementGradients);// need to be check
    }
    
    // Return
    return stresses;
  }
  
  //TODO: compute strains (i.e., linearized strain tensors OR deformation gradients) all quadrature points
  array<Strain, QP>
  computeStrainsAtGaussPoints(const NodalDisplacements & displacements) const {
    
    // TODO: Using the computeDispGradsAtGaussPoints function you just defined above,
    //       obtain the displacement gradients at all quadrature points
    const array<DisplacementGradient, QP> displacementGradients = computeDispGradsAtGaussPoints(displacements);
    
    // TODO: Sweep through all quadrature points, evaluate the strain tensor and save
    //       it in strains
    array<Stress, QP> strains;
    for (unsigned int qpIndex = 0; qpIndex < QP; ++qpIndex) {
      strains[qpIndex] = Strain::Zero();
      Matrix<double,SpatialDimension,SpatialDimension> standardDisplacementGradient
              = Utilities::convertTensorFromVoigtToStandard<SpatialDimension>(displacementGradients[qpIndex]);
      Matrix<double,SpatialDimension,SpatialDimension> standardStrain   = standardDisplacementGradient + Matrix<double,SpatialDimension,SpatialDimension>::Identity();
      strains[qpIndex] = Utilities::convertTensorFromStandardToVoigt<SpatialDimension>(standardStrain);
    }
    
    // Return
    return strains;
  }
  
  
  private:

  const Properties                              _properties;
  array<size_t, NumberOfNodes>                  _nodeIds;
  array<Point, NumberOfNodes>                   _nodePositions;
  const QuadratureRule<SpatialDimension, QP> *  _quadratureRule; //the weights and points, * means address of quadratureRule
  const MaterialModel *                         _materialModel;
  array<double, QP>                             _weightedJacobian;
  array<array<double, NumberOfNodes>, QP>       _shapeFunctions;
  array<array<Point, NumberOfNodes>, QP>        _shapeFunctionDerivatives;
  array<double, NumberOfNodes>                  _nodalWeights;// can ignore this
  AllBoundaries                                 _allBoundaries;
};

}

#endif //SIMPLICIAL_ELEMENT
