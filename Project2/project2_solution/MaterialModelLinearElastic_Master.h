// -*- C++ -*-
#ifndef MATERIAL_MODEL_LINEAR_ELASTIC
#define MATERIAL_MODEL_LINEAR_ELASTIC

#include "/src/Definitions.h"
#include "PJ2Utilities.h"

namespace MaterialModels {

class MaterialModelLinearElastic {

public:
    
  // TODO: Given that DisplacementGradient, Stress, Strain and TangentMatrix are in Voigt-notation,
  //       determine the number of rows and columns of each of those matrices
  // REMINDER: The second template parameter for Matrix is the number or rows, the third template
  //           parameter denotes the number of columns  
  typedef Matrix<double, 9, 1>  DisplacementGradient;
  typedef Matrix<double, 9, 1>  Stress;
  typedef Matrix<double, 9, 1>  Strain;
  typedef Matrix<double, 3, 3>  StandardStrainMatrix;
  typedef Matrix<double, 9, 9>  TangentMatrix;

  MaterialModelLinearElastic(const double youngsModulus, const double poissonsRatio):
    _lambda(poissonsRatio*youngsModulus/((1.0+poissonsRatio)*(1.0-2.0*poissonsRatio))),
    _mu    (youngsModulus/(2.0*(1.0+_mu))){
    }
  
  Strain
  computeStrain(const DisplacementGradient & displacementGradient) const {
    
    // TODO: Evaluate the strain vector (here: epsilon) based on the displacement gradient.
    // NOTE: The functionality to convert from standard to Voigt notation can be found
    //       in the Utilities namespace (see PJ2Utilities.h) (it might come in handy here...)
    const StandardStrainMatrix displacementGradientMatrix
      = Utilities::convertTensorFromVoigtToStandard<3>(displacementGradient);
    const StandardStrainMatrix strainMatrix
      = 0.5*(displacementGradientMatrix+displacementGradientMatrix.transpose());
    const Strain strain
      = Utilities::convertTensorFromStandardToVoigt<3>(strainMatrix);
      
    // Return
    return strain;
  };
  
  double
  computeEnergy(const DisplacementGradient & displacementGradient) const {
    
    // TODO: Evaluate strain in vector-form from displacementGradient using the
    //       computeStrain function you've defined before:
    Strain strainVector = computeStrain(displacementGradient);
    StandardStrainMatrix epsilon = Utilities::convertTensorFromVoigtToStandard<3>(strainVector);
    
    // TODO: If necessary, based on strain matrix, evaluate useful scalars such as the
    //       trace of the strain
    double traceOfStrain = epsilon.trace();
    
    // TODO: Evaluate the energy density
    double energyDensity = 0.5*_lambda*pow(traceOfStrain,2) + _mu*strainVector.transpose()*strainVector;
    
    // Return
    return energyDensity;
  };

  
  Stress
  computeStress(const DisplacementGradient & displacementGradient) const {
    
    // TODO: Evaluate strain in vector-form from displacementGradient using the
    //       computeStrain function you've defined before:
    Strain strainVector = computeStrain(displacementGradient);
    StandardStrainMatrix epsilon = Utilities::convertTensorFromVoigtToStandard<3>(strainVector);
    
    // TODO: If necessary, based on strain matrix, evaluate useful scalars such as the
    //       trace of the strain
    double traceOfStrain = epsilon.trace();
    
    // TODO: Evaluate the 2nd order stress tensor in Voigt-form (!)
    Stress stress
      = Utilities::convertTensorFromStandardToVoigt<3>(_lambda*traceOfStrain*StandardStrainMatrix::Identity())+ 2.0*_mu*strainVector;
    
    // Return
    return  stress;
  };
      
  TangentMatrix
  computeTangentMatrix(const DisplacementGradient & displacementGradient) const {
    
    // TODO: We keep our lives simple and define a 4th order tensor as an array first.
    //       1st: Set the right dimensions (i.e. replace the 1's by the correct dimension)
    //       2nd: Set all entries of tangentMatrixAsArray. You may have to use nested for-loops.
    array<array<array<array<double,3>,3>,3>,3> tangentMatrixAsArray;
    tangentMatrixAsArray[0][0][0][0] = 1.;
    
    
    for(unsigned int i = 0; i < 3; i++) {
      for(unsigned int j = 0; j < 3; j++){
        for(unsigned int k = 0; k < 3; k++){
          for(unsigned int l = 0; l < 3; l++){
            tangentMatrixAsArray[i][j][k][l] =
              _lambda * (i==j) * (k==l) + _mu *((i==k) * (j==l) + (i==l) * (j==k));
          }
        }
      }
    }
    
    
    // TODO: Simply replace the 1 again by the right dimension
    // We'll take care of the rest of the conversion. The functionality to convert between standard
    // and Voigt is found in the Utilities namespace in PJ2Utilities.h.
    TangentMatrix tangentMatrix
      = Utilities::convertFourthOrderTensorFromStandardToVoigt<3>(tangentMatrixAsArray);
    
    ignoreUnusedVariables(displacementGradient);
    
    return tangentMatrix;
    
  };
 


private:
  const double  _lambda;
  const double  _mu;
};


} // MaterialModels
#endif // MATERIAL_MODEL_LINEAR_ELASTIC
