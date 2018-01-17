#ifndef ELEMENT_TWO_NODE_BAR
#define ELEMENT_TWO_NODE_BAR

#include "Definitions.h"
#include "Utilities.h"

namespace Elements {

class Properties {
public:
  double _area, _density;
  Properties(const double area, const double density) :
    _area(area),
    _density(density){
  }
};

template<class MaterialModel, unsigned int Dimension>
class FiniteBar {

public:

  static const unsigned int  QuadPoints       = 1;
  static const unsigned int  NumberOfNodes    = 2;
  static const unsigned int  SpatialDimension = Dimension;
  static const unsigned int  DegreesOfFreedom = Dimension;
  static const VTKCellType   VtkCellType = VTK_LINE;
 
  // Typedef's using std and Eigen Classes
  typedef Matrix<double, SpatialDimension, 1>        Point;
  typedef Matrix<double, SpatialDimension, 1>        Vector;
  typedef Matrix<double,
                 DegreesOfFreedom,
                 DegreesOfFreedom>           				 NodeStiffnessMatrix;
  typedef Matrix<double, 
                 NumberOfNodes*DegreesOfFreedom,
                 NumberOfNodes*DegreesOfFreedom>     StiffnessMatrix;
  typedef Matrix<double,
                 NumberOfNodes*DegreesOfFreedom,
                 NumberOfNodes*DegreesOfFreedom>     MassMatrix;
  typedef array<Vector, NumberOfNodes>               NodalDisplacements;
  typedef array<Vector, NumberOfNodes>               Forces;
  
  // Typedef's derived from the MaterialModel
  typedef typename MaterialModel::Strain             Strain;
  typedef typename MaterialModel::Stress             Stress;
  typedef typename MaterialModel::TangentMatrix      TangentMatrix;  
  
  // Typedef based on the standard Node "NodeWithID" defined in Definitions.h
  typedef NodeWithId<Point>                          Node;

  // Public Members
  array<size_t, NumberOfNodes>  _nodeIds;
  Properties                    _properties;
  Point                         _X0,_X1;
  double                        _volume;
  double                        _lengthUndeformed;
  
  // TODO: in case you might want to store more public members - go for it!
  
  FiniteBar(  const array<Node, NumberOfNodes> &  nodes         ,            
              const Properties &                  properties    ,
              const MaterialModel *               materialModel ) :
      _properties   (properties)   ,
      _materialModel(materialModel){

   
    _volume = (nodes[0]._position-nodes[1]._position).norm()
                *properties._area;
      
    _nodeIds[0] = nodes[0]._id;
    _nodeIds[1] = nodes[1]._id;
  
    _X0     = nodes[0]._position;
    _X1     = nodes[1]._position;
  
    _lengthUndeformed = (_X1-_X0).norm();
  }
  
  MassMatrix
  computeConsistentMassMatrix() const {
      
      MassMatrix consistentMassMatrix;
      consistentMassMatrix.fill(0);
      
      // TODO: Fill in the consistent mass matrix      
      for (unsigned int i = 0; i < NumberOfNodes*DegreesOfFreedom; i++) {
        for (unsigned int j = 0; j < NumberOfNodes*DegreesOfFreedom; j++) {
          if (i==j){
            consistentMassMatrix(i,j) = 2.0;
          } else if ((i%DegreesOfFreedom)==(j%DegreesOfFreedom)) {
            consistentMassMatrix(i,j) = 1.0;
          }
        }
      }
      
      return consistentMassMatrix * _properties._density * _properties._area* _lengthUndeformed / 12.0;
  }
  
  double
  computeEnergy(const NodalDisplacements & displacements) const {
    
    Strain strain         = computeBarStrain(displacements);
    double energyDensity  = _materialModel->computeEnergy(strain);
    
    double energy = energyDensity * _properties._area * _lengthUndeformed;
    
    return energy;
  }

  Forces
  computeForces(const NodalDisplacements & displacements) const {
    
    Strain strain = computeBarStrain(displacements);
    Stress stress = _materialModel->computeStress(strain);
    
    Vector x0,x1;
    x0 = _X0 + displacements[0];
    x1 = _X1 + displacements[1];
    
    Vector forceVector = (x1-x0)/((x1-x0).norm());
    
    Forces forces;
    forces[0]= -forceVector * stress(0) * _properties._area;
    forces[1]=  forceVector * stress(0) * _properties._area;
    
    // Return                               
    return forces;
  }
  
  StiffnessMatrix
  computeStiffnessMatrix(const NodalDisplacements & displacements) const {
    
    StiffnessMatrix stiffnessMatrix;
    stiffnessMatrix.fill(0.);
      
    // TODO: this method computes the nodal forces for given nodal displacements
    Vector undeformedLength 		    = _X1-_X0;
    double undeformedLengthMagnitude= undeformedLength.norm();
    
    Vector deformedLength 			    = (_X1+displacements[1])-(_X0+displacements[0]);
    double deformedLengthMagnitude  = deformedLength.norm();
    Vector deformedLengthNormalized = deformedLength/deformedLength.norm(); 
    
    // Evaluate element strain, stress, tangent
    Strain strain 					= computeBarStrain(displacements);
    Stress stressElement 			= (*_materialModel).computeStress(strain);
    TangentMatrix tangentElement 	= (*_materialModel).computeTangentMatrix(strain);
    
    NodeStiffnessMatrix nodeStiffnessMatrix = 	
    _properties._area*(  tangentElement(0,0)/undeformedLengthMagnitude
              -stressElement(0,0) /deformedLengthMagnitude  )
      *deformedLengthNormalized*deformedLengthNormalized.transpose()
    +_properties._area*	 stressElement(0,0) /deformedLengthMagnitude
      *Matrix<double,DegreesOfFreedom,DegreesOfFreedom>::Identity();
    
    // TODO: this method computes the element stiffness matrix.
    for (unsigned int a = 0; a < NumberOfNodes; a++){
      for (unsigned int b = 0; b < NumberOfNodes; b++){
        for (unsigned int i=0; i < DegreesOfFreedom; i++){
          for (unsigned int j=0; j < DegreesOfFreedom; j++){
            stiffnessMatrix(a*DegreesOfFreedom+i,
                    b*DegreesOfFreedom+j)
              = (-1+2*(a==b))*nodeStiffnessMatrix(i,j);
          }
        }
      }
    }
	
    return stiffnessMatrix;
  }

  Strain
  computeBarStrain(const NodalDisplacements & displacements) const {
    
    Strain strain = Strain::Zero();
    
    double lengthDeformed = ( _X0+displacements[0] - (_X1+displacements[1]) ).norm();
    strain(0)             = (lengthDeformed-_lengthUndeformed)/_lengthUndeformed;
    
    return strain;
  }
  
  Stress
  computeBarStress(const NodalDisplacements & displacements) const {

    Stress stress;
    Strain strain = computeBarStrain(displacements);

    stress = (*_materialModel).computeStress(strain);

    return stress;
  }

  array<size_t, 2>
  getNodeIds() const {
    return _nodeIds;
  }

  array<Strain, QuadPoints>
  computeStrainsAtGaussPoints(const NodalDisplacements & displacements) const {
    array<Strain,1> strainsAtGaussPoints;
    strainsAtGaussPoints[0] = computeBarStrain(displacements);
    return strainsAtGaussPoints;
  }
  
  array<Stress, QuadPoints>
  computeStressesAtGaussPoints(const NodalDisplacements & displacements) const {
    array<Stress,1> stressesAtGaussPoints;
    stressesAtGaussPoints[0] = computeBarStress(displacements);
    return stressesAtGaussPoints;
  }
  
private:
  
  // Private members    
  const MaterialModel * _materialModel;
  
  
};

}

#endif //ELEMENT_TWO_NODE_BAR
