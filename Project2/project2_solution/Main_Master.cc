#include "/src/Definitions.h"
#include "MaterialModelHyperelasticity_Master.h"
#include "MaterialModelLinearElastic_Master.h"
#include "MaterialModelTests_Master.h"
#include "PJ2Utilities.h"

// TODO: The typedef for the linear elastic model has been graciously provided here
//       Note, that for the finite deformation material model, you will have to change
//       this line
typedef MaterialModels::MaterialModelLinearElastic MaterialModel;

int main() {

  printf("Testing MaterialModelLinearElastic in %1u dimensions \n\n", 3);

  // Define constants that your material model will use, for example:
  /* const double lambda = 1.33;
  const double mu     = 1.67; */
  const double E  = 1.00;
  const double nu = 0.30;

  // Create a material model:
  MaterialModel model(E,nu);
  // MaterialModel model(lambda,mu);

  // We feed your newly created materialModel to the material model test
  // function defined in the namespace MaterialModels to confirm that it is correct
  MaterialModels::testMaterialModelDerivatives<MaterialModel>(model); 

  //REMINDER: Functions from a namespace are accessed via NameOfNamespace::FunctionName<class>(object)

  return 0;
}
