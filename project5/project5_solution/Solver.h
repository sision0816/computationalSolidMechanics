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
    
    VectorXd solution(numberOfDOFs);
    for (unsigned int nodeIndex = 0; nodeIndex < _assembler0.getNumberOfNodes(); nodeIndex++) {
      for (unsigned int dofIndex = 0; dofIndex < DegreesOfFreedom; dofIndex++) {
        solution(nodeIndex * DegreesOfFreedom + dofIndex) =
          initialGuess[nodeIndex](dofIndex);
      }
    }
   
   
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
      
                                                                                            
                                                                                           
                                                                                            
                                                         
      
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      
      const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
      solution(dofIndex)    = bc._constraint;
      
    }

    // make nodal-wise displacements
    vector<ElementVector> nodalDisplacements =
      Utilities::distributeGlobalVectorToLocalVectors<Assembler0>(solution);
    
    
    
    
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%             Compute global force vector        %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    
    // use the nodal displacements as well as the assembleForceVector functionality you
    // implemented in your assembler, to obtain the global force
    VectorXd globalForces =
       _assembler0.assembleForceVector(nodalDisplacements)
      +_assembler1.assembleForceVector(nodalDisplacements);

    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
     
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      const size_t dofIndex   = bc._nodeId * DegreesOfFreedom + bc._coordinate;
      globalForces(dofIndex)  = 0;
     
    }
    
    
    
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%             Setup of iterative loop            %%%%%%%%%%%%%%%%%%%%
    // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
    
    // evaluate residue based on the norm of globalForces
    double residue                  = globalForces.norm();
    unsigned int numberOfIterations = 0;
    
    while(residue > tolerance && numberOfIterations < maxIterations) {

    
      // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
      // %%%%%%%%%%%%%%%%%%%%        Compute global stiffness matrix         %%%%%%%%%%%%%%%%%%%%
      // %%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%
      
      // use the nodal displacements as well as the assembleStiffnessMatrix functionality you
      // implemented in your assembler, to obtain the tangentMatrix
      MatrixXd tangentMatrix =
         _assembler0.assembleStiffnessMatrix(nodalDisplacements)
        +_assembler1.assembleStiffnessMatrix(nodalDisplacements);
      
      //   apply boundary conditions to stiffness matrix
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
     
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        
        const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;

        tangentMatrix.row(dofIndex).fill(0);

        tangentMatrix(dofIndex, dofIndex) = 1;
       
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

      // use the nodal displacements as well as the assembleForceVector functionality you
      // implemented in your assembler, to obtain the global force
      globalForces =
         _assembler0.assembleForceVector(nodalDisplacements)
        +_assembler1.assembleForceVector(nodalDisplacements);

      // zero out entries of global forces that are not considered
       for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
                
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        const size_t dofIndex   = bc._nodeId * DegreesOfFreedom + bc._coordinate;
        globalForces(dofIndex)  = 0;
        
      }
      
      // Tevaluate residue based on the norm of globalForces
      residue = globalForces.norm();
      if (verbose == true) {
        printf("Newton Raphson iteration %4u, residue = %8.3e\n",
               numberOfIterations, residue);
      }

      // Increase numberOfIterations by one
      numberOfIterations++;
      
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
