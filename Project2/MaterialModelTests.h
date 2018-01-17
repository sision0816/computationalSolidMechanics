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

    DisplacementGradient displacementGradient
      = DisplacementGradient::Random();
    unsigned int numberOfDisplacementGradientComponents
      = displacementGradient.size();

    //TODO: compute the exact energy, stresses and tangent matrix using the input object "model"
    double 	      energy  = model.computeEnergy(displacementGradient);
    Stress 	      stress = model.computeStress(displacementGradient);
    TangentMatrix tangent = model.computeTangentMatrix(displacementGradient);
    
    cout << "[MMT] Materials::testMaterialModelDerivatives started: " << endl;
    
    //TODO: define the perturbation and tolerance values
    double perturbation = 0.000001;
    double tolerance    = 0.0001;

    
    
    // ===================================
    // === TEST : STRESS - ANA vs. NUM ===
    // ===================================
    
    Stress numericalStresses = Stress::Zero();    
    //TODO: recompute stresses from the material model energy by taking numerical derivatives
    for (	unsigned int index_I = 0; index_I < numberOfDisplacementGradientComponents; index_I++ ){
		

      // Evaluate perturbated energy
      // Add perturbation
      displacementGradient(index_I,0) += perturbation;

      // Evaluate energy of perturbed displacement state
      double energyOfPertubedDisplacementState = model.computeEnergy(displacementGradient);
	  
      // Remove perturbation
      displacementGradient(index_I,0) -= perturbation;

      // Evaluate numerical stress approximation
	  
      numericalStresses(index_I,0) = (energyOfPertubedDisplacementState - energy) / perturbation;
	  

      printf("[MMT] Stress (Ana): %8.5f | Stress (Num): %8.5f\n",
              stress           (index_I,0),
              numericalStresses(index_I,0));
    }

    // TODO: compute Frobenius norm
    double errorStresses = sqrt((numericalStresses - stress).dot(numericalStresses - stress)) / sqrt(stress.dot(stress));
    cout << "[MMT] Error of method computeStress = " << errorStresses << endl;


    
    // ====================================
    // === TEST : TANGENT - ANA vs. NUM ===
    // ====================================
    
    TangentMatrix numericalTangent = TangentMatrix::Zero();
    //TODO: recompute tangent matrix from the material model stresses by taking numerical derivatives
    for (	unsigned int index_I = 0; index_I < numberOfDisplacementGradientComponents; index_I++ ){

      // Add perturbation
      displacementGradient(index_I,0) += perturbation;

      // Evaluate perturbated stress
	  Stress stressOfPertubedDisplacementState = Stress::Zero();
      stressOfPertubedDisplacementState = model.computeStress(displacementGradient);

      // Remove perturbation
      displacementGradient(index_I,0) -= perturbation;

      // Evaluate numerical Tangent approximation
	  
      for (	unsigned int index_J = 0;
        index_J < numberOfDisplacementGradientComponents;
        index_J++){

        numericalTangent(index_J,index_I) = (stressOfPertubedDisplacementState(index_J,0) - stress(index_J,0)) / perturbation;
		

        printf("[MMT] (%d,%d) Tangent (Ana): %8.5f | Tangent (Num): %8.5f\n",
                index_I,
                index_J,
                tangent         (index_J,index_I),
                numericalTangent(index_J,index_I) );	
      }

    }

    //TODO: compute error between analytical and numerical tangent matrices
    //double errorTangent = (tangent - numericalTangent).norm() / tangent.norm();
	//double errorTangent = (tangent - numericalTangent).determinant() / tangent.determinant();
    Matrix<double, 9, 9> tangentMinusNumericalTangent = tangent - numericalTangent;
    double tangentMinusNumericalTangentDotTangentMinusNumericalTangent = 0.0;
    double tangentDotTangent = 0.0;
    for (int i = 0; i<9; i++){
      for (int j=0; j<9; j++){
        tangentMinusNumericalTangentDotTangentMinusNumericalTangent += tangentMinusNumericalTangent(i,j)*tangentMinusNumericalTangent(i,j);
        tangentDotTangent += tangent(i,j)*tangent(i,j);
      }
    }
    double normOfTangentMinusNumericalTangent = sqrt(tangentMinusNumericalTangentDotTangentMinusNumericalTangent);
    double normOfTangent =sqrt(tangentDotTangent);
    double errorTangent = normOfTangentMinusNumericalTangent / normOfTangent;
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
	/*if (errorStresses >= tolerance ) {
      cout << "[MMT] ================================================" << endl;
      cout << "[MMT] Warning: material model stress derivatives test failed." << endl << endl;
      cout << "[MMT] ================================================" << endl;
      return false;
    }
	if (errorTangent >= tolerance ) {
      cout << "[MMT] ================================================" << endl;
      cout << "[MMT] Warning: material model Tangent derivatives test failed." << endl << endl;
      cout << "[MMT] ================================================" << endl;
      return false;
    }*/
    else {
      cout << endl << "[MMT] =======================================" << endl;
      cout << "[MMT] Material model derivatives test passed." << endl;
      cout << "[MMT] =======================================" << endl << endl;
      return true;
    }

    ignoreUnusedVariables(energy,perturbation);
    
    return true;
  }



}

#endif  // MATERIALMODEL_DERIVATIVETEST_H
