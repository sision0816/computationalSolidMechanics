#include "/src/Definitions.h"
#include "MaterialModelHyperelasticity.h"
#include "MaterialModelLinearElastic.h"
#include "MaterialModelTests.h"
#include "PJ2Utilities.h"

// TODO: The typedef for the linear elastic model has been graciously provided here
//       Note, that for the finite deformation material model, you will have to change
//       this line
typedef MaterialModels::MaterialModelLinearElastic LinearElasticModel;
typedef MaterialModels::MaterialModelHyperelastic HyperElasticMaterialModel;

int main() {

  printf("Testing MaterialModelLinearElastic in %1u dimensions \n\n", 3);

  // Define constants that your material model will use, for example:
  const double lambda = 1.33;
  const double mu = 1.67;
  const double youngsModulus = 4.0;
  const double poissonsRatio = 0.23; 

  // Create a material model:
  LinearElasticModel linearElasticModel(lambda,mu);
  HyperElasticMaterialModel hyperElasticModel(youngsModulus,poissonsRatio);

  // We feed your newly created materialModel to the material model test
  //function defined in the namespace MaterialModels to confirm that it is correct
  MaterialModels::testMaterialModelDerivatives<LinearElasticModel>(linearElasticModel);
  MaterialModels::testMaterialModelDerivatives<HyperElasticMaterialModel>(hyperElasticModel);
  
  //REMINDER: Functions from a namespace are accessed via NameOfNamespace::FunctionName<class>(object)


  return 0;
}
