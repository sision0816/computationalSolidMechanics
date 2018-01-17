#ifndef SOLVER_CLASS_H
#define SOLVER_CLASS_H
#include "Definitions.h"
#include "Utilities.h"
#include "Assembler.h"
#include <Eigen/LU>
  
template <class Assembler0,class Assembler1>
class Solver {

  public:

  typedef typename Assembler0::ElementVector ElementVector;

  // "Special-Case"-Constructor
  Solver(Assembler0 & assembler0):
      _assembler0(assembler0),
      _assembler1(Assembler1()){
  }
  
  // "General-Case"-Constructor
  Solver(Assembler0 & assembler0, Assembler1 & assembler1):
      _assembler0(assembler0),_assembler1(assembler1){
  }
  
  // Newton-Raphon-Solver
  vector<ElementVector>
  computeNewtonRaphsonSolution(const vector<EssentialBoundaryCondition> & essentialBCs          ,
                               const vector<ElementVector> &              initialGuess          ,
                               const unsigned int                         maxIterations = 1000  ,
                               const double                               tolerance     = 1e-8  ,
                               const bool                                 verbose       = false ){

    if (verbose == true) {
      printf("Newton Raphson solver trying to achieve a tolerance of %e in %u "
           "maximum iterations\n", tolerance, maxIterations);
    }

    const unsigned int DegreesOfFreedom = Assembler0::DegreesOfFreedom;
    const unsigned int numberOfDOFs     = _assembler0.getNumberOfNodes() *DegreesOfFreedom;

    
    
    
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%            Start with an initial guess         %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    
    // Todo: keeping in mind, that initialGuess is an array of displacements at all nodes,
    //       load it onto solution, which is the ginormous vector of length
    //       DegreesOfFreedom * _assembler0.getNumberOfNodes()
    VectorXd solution(numberOfDOFs);
    for (unsigned int nodeIndex = 0; nodeIndex < _assembler0.getNumberOfNodes(); nodeIndex++) {
      for (unsigned int dofIndex = 0; dofIndex < DegreesOfFreedom; dofIndex++) {
        solution(nodeIndex * DegreesOfFreedom + dofIndex) = initialGuess[nodeIndex](dofIndex); // ...
      }
    }
   
   
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
      
      // NOTE: Much of the following structure will be used at subsequent parts of the code,
      //       where essential boundary conditions are "read out" and "implemented". So the
      //       next line of code can for example pretty much be copied into every single one
      //       of these parts of the code. You'll see! :)
      
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];

      // you can delete the following three lines - their sole purpose is to illustrate how
      // you read from essential boundary conditions
      const size_t someNodeID     = bc._nodeId;
      const size_t someCoordinate = bc._coordinate;
      const double someConstraint = bc._constraint;
      
      ignoreUnusedVariables(someNodeID,someCoordinate,someConstraint);
      
      // ToDo:
      // set the right entries of solution so as to impose the right boundary constraints
      const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
      solution(dofIndex)    = bc._constraint;
      
    }

    // make nodal-wise displacements
    vector<ElementVector> nodalDisplacements =
      Utilities::distributeGlobalVectorToLocalVectors<Assembler0>(solution);
    
    
    
    
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%             Compute global force vector        %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    
    // ToDo: use the nodal displacements as well as the assembleForceVector functionality you
    //       implemented in your assembler, to obtain the global force.
    // Note: Don't forget, that you have two assemblers...
    VectorXd globalForces = _assembler0.assembleForceVector (nodalDisplacements) + _assembler1.assembleForceVector (nodalDisplacements);

    
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
     
     // ToDo: same spiel as above - call every boundary condition, read out the equivalent 
     //       dofIndex and set the corresponding entry in glob
      //
      //
      // alForces equal to 0.0
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      globalForces(bc._nodeId * DegreesOfFreedom + bc._coordinate) = 0.0;
     
    }

    
    
    
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%             Setup of iterative loop            %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    
    // ToDo: evaluate residue based on the norm of globalForces
    double residue                  = globalForces.norm();
    unsigned int numberOfIterations = 0;
    
    while(residue > tolerance && numberOfIterations < maxIterations) {

    
      // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
      // %%%%%%%%%%%%%%%%%%%%        Compute global stiffness matrix         %%%%%%%%%%%%%%%%%%%%
      // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
      
      // ToDo: use the nodal displacements as well as the assembleStiffnessMatrix functionality you
      //       implemented in your assembler, to obtain the tangentMatrix.
      // Note: Don't forget, that you have two assemblers...
      MatrixXd tangentMatrix = _assembler0.assembleStiffnessMatrix(nodalDisplacements) + _assembler1.assembleStiffnessMatrix(nodalDisplacements);
      
      //   apply boundary conditions to stiffness matrix
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
     
       // ToDo: same spiel as above - call every boundary condition, read out the equivalent 
       //       dofIndex and now ...
       //       ... first zero the entire dofIndexth row (Note: tangentMatrix.row(3).fill(1.0);
       //           fills the 4th row with 1.0s ... you may want to use this functionality :) )
       //       ... set the dofIndexth entry of the dofIndexth row equal to 1
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        tangentMatrix.row(bc._nodeId * DegreesOfFreedom + bc._coordinate).fill(0.0);
        tangentMatrix(bc._nodeId * DegreesOfFreedom + bc._coordinate, bc._nodeId * DegreesOfFreedom + bc._coordinate) = 1.0;
      }

      
      
      
      // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
      // %%%%%%%%%%%%%%%%%%%%   TangentMatrix check (nothing to do here)     %%%%%%%%%%%%%%%%%%%%
      // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
      
      vector<size_t> zeroRows;
      for (size_t row = 0; row < size_t(tangentMatrix.rows()); ++row) {
        if (tangentMatrix.row(row).norm() < 1e-3){
          zeroRows.push_back(row);
        }
      }
      if (zeroRows.size() > 0) {
        printf("Error: Newton Raphson solver found that on iteration %u, "
               "the tangent matrix has %zu rows full of zeros.  The indices are: \n",
               numberOfIterations, zeroRows.size());
        for (size_t i = 0; i < zeroRows.size(); ++i) {
          printf("%zu, ", zeroRows[i]);
        }
        printf("\nSolver stops.\n");
        exit(1);
      }

      
      
      
      // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
      // %%%%%%%%%%%%%%%%%%%%        Perform Newton-Raphson interation       %%%%%%%%%%%%%%%%%%%%
      // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
      
      // update solution (nothing to be done here, but please do read this line carefully :) )
      solution -= tangentMatrix.lu().solve(globalForces);
      
      // reassemble displacements from solution
      nodalDisplacements =
        Utilities::distributeGlobalVectorToLocalVectors<Assembler0>(solution);

      // ToDo: use the nodal displacements as well as the assembleForceVector functionality you
      //       implemented in your assembler, to obtain the global force.
      // Note: Don't forget, that you have two assemblers..
      globalForces = _assembler0.assembleForceVector (nodalDisplacements) + _assembler1.assembleForceVector (nodalDisplacements);

      // zero out entries of global forces that are not considered
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
        
        // ToDo: same spiel as above - call every boundary condition, read out the equivalent 
        //       dofIndex and set the corresponding entry in globalForces equal to 0.0

        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        globalForces(bc._nodeId * DegreesOfFreedom + bc._coordinate) = 0.0;
        
      }
      
      // ToDo: evaluate residue based on the norm of globalForces
      residue = globalForces.norm();
      if (verbose == true) {
        printf("Newton Raphson iteration %4u, residue = %8.3e\n",
               numberOfIterations, residue);
      }

      // ToDo:
      // Increase numberOfIterations by one
      
      numberOfIterations += 1;
      
    }
    
    
    
    
    // %%%%%%%%%%%%%%%%%%%%                          %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%        FINAL STEPS       %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%                          %%%%%%%%%%%%%%%%%%%%
    
    if (numberOfIterations == maxIterations) {
      printf("Error: Newton Raphson solver could not converge in %u "
           "iterations.\n", maxIterations);
      exit(1);
    }

    // finally, we transform the global vector back into local element vectors
    vector<ElementVector> localVectors =
      Utilities::distributeGlobalVectorToLocalVectors<Assembler0>(solution);

    return localVectors;
    
  }


  private:
  Assembler0 _assembler0;
  Assembler1 _assembler1;

};

#endif  // SOLVER_CLASS_H
