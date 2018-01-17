// -*- C++ -*-
#ifndef MATERIAL_MODEL_HYPERELASTIC
#define MATERIAL_MODEL_HYPERELASTIC

#include "/src/Definitions.h"
#include "PJ2Utilities.h"

namespace MaterialModels {

class MaterialModelHyperelastic {

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

  // Constructor
  MaterialModelHyperelastic(const double lambda, const double mu):
    // TODO: Set the Lame parameters _lambda and _mu based on the input
    // CAUTION: _lambda and _mu are endowed with the 'const'-qualifier and hence can only
    //          be defined in this way as part of the constructor.
    _lambda(lambda),
    _mu    (mu    ){}

  
  Strain
  computeStrain(const DisplacementGradient & displacementGradient) const {
    
    // Finite Deformations: "strainVector" corresponds to F, i.e. I+\grad u
    Strain strainVector
      = displacementGradient
        + Utilities::convertTensorFromStandardToVoigt<3>(StandardStrainMatrix::Identity());
    return strainVector;
    
  };
  
  double
  computeEnergy(const DisplacementGradient & displacementGradient) const {

    // Evaluate F based on displacementGradient using computeStrain
    Strain strainVector = computeStrain(displacementGradient);
    StandardStrainMatrix F
      = Utilities::convertTensorFromVoigtToStandard<3>(strainVector);
      
    // TODO: If necessary, evaluate useful scalars such as the I1 invariant
    //       and the volumetric expansion J based on F.
    StandardStrainMatrix C
      = F.transpose()*F;
    double I1 = C.trace();
    double J  = F.determinant();
    
    // TODO: Evaluate the energy density
    double energyDensity = 0.5*_mu*(I1-3.0)+0.5*_lambda*pow(log(J),2) - _mu *log(J);
    
    // Return
    return energyDensity;
  };

  
  Stress
  computeStress(const DisplacementGradient & displacementGradient) const {
    
    // Evaluate F based on displacementGradient using computeStrain
    Strain strainVector = computeStrain(displacementGradient);
    StandardStrainMatrix F
      = Utilities::convertTensorFromVoigtToStandard<3>(strainVector);
      
    // TODO: If necessary, evaluate useful scalars such as the I1 invariant
    //       and the volumetric expansion J or tensors such as the right Cauchy-Green
    //       tensor C or the inverse-transpose of F
    StandardStrainMatrix FInverseTranspose
      = (F.inverse()).transpose();
    Strain strainVectorInvTrans
      = Utilities::convertTensorFromStandardToVoigt<3>(FInverseTranspose);
    double J  = F.determinant();
    
   // TODO: Based on the possibly-useful scalars and tensors, you've obtained above,
    //       determine the stress tensor.
    Stress stress = _mu * strainVector
                    + _lambda * log(J) * strainVectorInvTrans
                    - _mu * strainVectorInvTrans;
    
    // Return
    return stress;
  };

  // method computeTangentMatrix, which receives a DisplacementGradient vector and returns a TangentMatrix C:

  TangentMatrix
  computeTangentMatrix(const DisplacementGradient & displacementGradient) const {
    
    // Evaluate F based on displacementGradient using computeStrain
    Strain strainVector = computeStrain(displacementGradient);
    StandardStrainMatrix F
      = Utilities::convertTensorFromVoigtToStandard<3>(strainVector);
      
    // TODO: If necessary, evaluate useful scalars such as the I1 invariant
    //       and the volumetric expansion J or tensors such as the right Cauchy-Green
    //       tensor C or the inverse-transpose of F
    StandardStrainMatrix FInvTrans
      = (F.inverse()).transpose();
    double J  = F.determinant();
    
    
    // TODO: We keep our lives simple and define a 4th order tensor as an array first.
    //       1st: Set the right dimensions (i.e. replace the 1's by the correct dimension)
    //       2nd: Set all entries of tangentMatrixAsArray. You may have to use nested for-loops.
    array<array<array<array<double,3>,3>,3>,3> tangentMatrixAsArray;
    for(unsigned int i = 0; i < 3; i++) {
      for(unsigned int j = 0; j < 3; j++){
        for(unsigned int k = 0; k < 3; k++){
          for(unsigned int l = 0; l < 3; l++){
            tangentMatrixAsArray[i][j][k][l]
              = _mu*(i==k)*(j==l)
                +_lambda*( FInvTrans(k,l)*FInvTrans(i,j) - log(J)*FInvTrans(k,j)*FInvTrans(i,l) )
                + _mu*FInvTrans(k,j)*FInvTrans(i,l);
          }
        }
      }
    }
    
    // TODO: Simply replace the 1 again by the right dimension
    // We'll take care of the rest of the conversion. The functionality to convert between standard
    // and Voigt is found in the Utilities namespace in PJ2Utilities.h.
    TangentMatrix tangentMatrix
      = Utilities::convertFourthOrderTensorFromStandardToVoigt<3>(tangentMatrixAsArray);
    
    // Return
    return tangentMatrix;
  };

  private:
  const double _lambda;
  const double _mu;

  };
}
#endif // MATERIAL_MODEL_HYPERELASTIC
