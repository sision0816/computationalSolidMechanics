// -*- C++ -*-
#ifndef MATERIALMODEL_DERIVATIVETEST_H
#define MATERIALMODEL_DERIVATIVETEST_H

#include "/src/Definitions.h"
#include "PJ2Utilities.h"
 
namespace MaterialModels {

  template <class MaterialModel>
  bool
  testMaterialModelDerivatives(const MaterialModel & model) {

    // Typedefs are created through copying from the existing ones in the MaterialModel
    typedef typename MaterialModel::DisplacementGradient DisplacementGradient;
    typedef typename MaterialModel::Stress               Stress;
    typedef typename MaterialModel::TangentMatrix        TangentMatrix;

    DisplacementGradient displacementGradient = 0.1*DisplacementGradient::Random();
    unsigned int numberOfDisplacementGradientComponents = displacementGradient.size();
    
    // Exact energy, stress, tangent
    double        energy  = model.computeEnergy(displacementGradient);
    Stress        stress  = model.computeStress(displacementGradient);
    TangentMatrix tangent = model.computeTangentMatrix(displacementGradient);

    cout << "[MMT] Materials::testMaterialModelDerivatives started: " << endl;
    
    double perturbation = 1.0e-6;
    double tolerance    = 1.0e-4;

    

    // ===================================
    // === TEST : STRESS - ANA vs. NUM ===
    // ===================================
     
    Stress numericalStresses = Stress::Zero();
    for (unsigned int index_I = 0; index_I < numberOfDisplacementGradientComponents; index_I++ ){

      // Evaluate perturbated energy
      
      // Add perturbation
      displacementGradient(index_I,0)
        += perturbation; 
        
      // Evaluate energy of perturbed displacement state
      double energyPertubated
        = model.computeEnergy(displacementGradient);
        
      // Remove perturbation
      displacementGradient(index_I,0)
        -= perturbation;

      // Evaluate numerical stress approximation
      numericalStresses(index_I,0) = (energyPertubated-energy)/perturbation;

      printf("[MMT] Stress (Ana): %8.5f | Stress (Num): %8.5f\n",
        stress           (index_I,0),
        numericalStresses(index_I,0));
    }

    // We evaluate the Frobenius norm
    double errorStresses = (stress - numericalStresses).norm()/stress.norm();
    cout << "[MMT] Error of method computeStress = " << errorStresses << endl;



    // ====================================
    // === TEST : TANGENT - ANA vs. NUM ===
    // ====================================
    
    TangentMatrix numericalTangent = TangentMatrix::Zero();
    for (unsigned int index_I = 0; index_I < numberOfDisplacementGradientComponents; index_I++){

      // Evaluate perturbated energy
      displacementGradient(index_I,0)
        += perturbation; 
      Stress stressPertubated
        = model.computeStress(displacementGradient);
      displacementGradient(index_I,0)
        -= perturbation;

      // Evaluate numerical stress approximation
      for (unsigned int index_J = 0; index_J < numberOfDisplacementGradientComponents; index_J++){

        numericalTangent(index_J,index_I)
          = (stressPertubated(index_J,0)-stress(index_J,0))/perturbation;
        printf("[MMT] (%d,%d) Tangent (Ana): %8.5f | Tangent (Num): %8.5f\n",
          index_I,
          index_J,
          tangent(index_J,index_I),
          numericalTangent(index_J, index_I) );	
      }

    }

    // Error between analytical and numerical tangent matrices
    double errorTangent = (tangent - numericalTangent).norm()/tangent.norm();
    cout << "[MMT] Error of method computeTangentMatrix = " << errorTangent << endl;

    
    
    // =============================
    // === FINAL OUTPUT & RETURN ===
    // =============================
    
    if (errorStresses >= tolerance || errorTangent >= tolerance) {
      cout << "[MMT] ================================================" << endl;
      cout << "[MMT] Warning: material model derivatives test failed." << endl << endl;
      cout << "[MMT] ================================================" << endl;
      return false;
    }
    else {
      cout << endl << "[MMT] =======================================" << endl;
      cout << "[MMT] Material model derivatives test passed." << endl;
      cout << "[MMT] =======================================" << endl << endl;
      return true;
    }

    return true;
  }



}

#endif  // MATERIALMODEL_DERIVATIVETEST_H
