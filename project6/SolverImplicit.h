#ifndef SOLVER_CLASS_H
#define SOLVER_CLASS_H
#include "Definitions.h"
#include "Utilities.h"
#include "Assembler.h"
#include <Eigen/LU>
  
template <class Assembler0,class Assembler1, class Assembler2>
class SolverImplicitDynamics {

  public:

  typedef typename Assembler0::ElementVector ElementVector;

  // "Special-Case"-Constructor
  SolverImplicitDynamics(Assembler0 & assembler0,
                         const double timestep    ,
                         const double dampingAlpha,
                         const double dampingBeta ,
                         const double newmarkBeta ,
                         const double newmarkGamma):
      _assembler0(assembler0  ),
      _assembler1(Assembler1(assembler0.getNumberOfNodes())),
      _assembler2(Assembler2(assembler0.getNumberOfNodes())),
      _dampingAlpha(dampingAlpha),
      _dampingBeta (dampingBeta ),
      _newmarkBeta (newmarkBeta ),
      _newmarkGamma(newmarkGamma){
  }
  
  // "Special-Case"-Constructor
  SolverImplicitDynamics(Assembler0 & assembler0,
                         Assembler1 & assembler1,
                         const double timestep    ,
                         const double dampingAlpha,
                         const double dampingBeta ,
                         const double newmarkBeta ,
                         const double newmarkGamma):
      _assembler0(assembler0),
      _assembler1(assembler1),
      _assembler2(Assembler2(assembler0.getNumberOfNodes())), 
      _timestep    (timestep    ), 
      _dampingAlpha(dampingAlpha),
      _dampingBeta (dampingBeta ),
      _newmarkBeta (newmarkBeta ),
      _newmarkGamma(newmarkGamma){
  }
  
  // "General-Case"-Constructor
  SolverImplicitDynamics(Assembler0 & assembler0,
                         Assembler1 & assembler1,
                         Assembler2 & assembler2,
                         const double timestep    ,
                         const double dampingAlpha,
                         const double dampingBeta ,
                         const double newmarkBeta ,
                         const double newmarkGamma):
      _assembler0(assembler0),
      _assembler1(assembler1),
      _assembler2(assembler2),
      _timestep    (timestep    ),
      _dampingAlpha(dampingAlpha),
      _dampingBeta (dampingBeta ),
      _newmarkBeta (newmarkBeta ),
      _newmarkGamma(newmarkGamma){
  }
  
  
  
  void
  computeNewmarkUpdate(const vector<EssentialBoundaryCondition> & essentialBCs,
                       vector<ElementVector> &       currentNodalDisplacement ,
                       vector<ElementVector> &       currentNodalVelocity     ,
                       vector<ElementVector> &       currentNodalAcceleration ,
                       const unsigned int            maxIterations = 1000,
                       const double                  tolerance = 1e-4    ,
                       const bool                    verbose = true      ) {
                         
    // Solving for the updated displacements using Newton-Raphson iterations
    if (verbose) {
      printf("Implicit Dynamics solver trying to achieve a tolerance of %e in %u "
             "maximum iterations\n", tolerance, maxIterations);
    }

    // Some parameters
    size_t DegreesOfFreedom = Assembler0::DegreesOfFreedom;
    size_t numberOfDOFs     = currentNodalDisplacement.size()*DegreesOfFreedom;

    // TODO: create three VectorXd's currentDisplacement, currentVelocity and currentAcceleration
    //       these should remain unchanged so you can even define them as "const"
    
    VectorXd currentDisplacement(numberOfDOFs);
	for (unsigned int nodeIndex = 0; nodeIndex < currentNodalDisplacement.size(); nodeIndex++) {
      for (unsigned int dofIndex = 0; dofIndex < DegreesOfFreedom; dofIndex++) {
        currentDisplacement(nodeIndex * DegreesOfFreedom + dofIndex) = currentNodalDisplacement[nodeIndex](dofIndex); // ...
      }
    }
    VectorXd currentVelocity(numberOfDOFs);
	for (unsigned int nodeIndex = 0; nodeIndex < currentNodalDisplacement.size(); nodeIndex++) {
      for (unsigned int dofIndex = 0; dofIndex < DegreesOfFreedom; dofIndex++) {
        currentVelocity(nodeIndex * DegreesOfFreedom + dofIndex) = currentNodalVelocity[nodeIndex](dofIndex); // ...
      }
    }
	VectorXd currentAcceleration(numberOfDOFs);
    for (unsigned int nodeIndex = 0; nodeIndex < currentNodalDisplacement.size(); nodeIndex++) {
      for (unsigned int dofIndex = 0; dofIndex < DegreesOfFreedom; dofIndex++) {
        currentAcceleration(nodeIndex * DegreesOfFreedom + dofIndex) = currentNodalAcceleration[nodeIndex](dofIndex); // ...
      }
    }
    // ...
    
    // TODO: further define a VectorXd newDisplacement, which is the solution, we will
    //       iteratively (try to) improve
    
    VectorXd newDisplacement(numberOfDOFs);
	newDisplacement.fill(0);
    
    // TODO: Boundary conditions I - Solution
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
       
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
      currentDisplacement(dofIndex)    = bc._constraint;
      
    }

    // TODO: newDisplacement in vector<ElementVector> form - needed in this form to evaluate
    //       stiffness, forces, etc. etc.
    vector<ElementVector> nodalNewDisplacements
        	=  Utilities::distributeGlobalVectorToLocalVectors<Assembler0>(currentDisplacement);
    
    // TODO: Evaluate the consistent mass matrix and damping matrix
    Eigen::MatrixXd consistentMassMatrix(numberOfDofs, numberOfDofs);
	consistentMassMatrix = _assembler0.assembleConsistentMassMatrix() + _assembler1.assembleConsistentMassMatrix() + _assembler2.assembleConsistentMassMatrix();
	Eigen::MatrixXd stiffnessMatrix(numberOfDofs, numberOfDofs);
	stiffnessMatrix = _assembler0.assembleStiffnessMatrix(nodalNewDisplacements) + _assembler1.assembleStiffnessMatrix(nodalNewDisplacements) + _assembler2.assembleStiffnessMatrix(nodalNewDisplacements);
	Eigen::MatrixXd dampingMatrix(numberOfDofs, numberOfDofs);
	dampingMatrix = dampingAlpha * consistentMassMatrix + dampingBeta * stiffnessMatrix;
    
    //TODO: evaluate the effective force vector
    VectorXd effectiveForceVector;
    VectorXd globalForceVector;// ...
    VectorXd rVector;
    globalForceVector = _assembler0.assembleForceVector(nodalNewDisplacements) + _assembler1.assembleForceVector(nodalNewDisplacements) + _assembler2.assembleForceVector(nodalNewDisplacements);
    rVector = consistentMassMatrix.solve((1/(_newmarkBeta*_timestep*_timestep) * currentNodalDisplacement) + 1/(_newmarkBeta*_timestep) * currentNodalVelocity + (1/(2*_newmarkBeta) -1)*currentNodalAcceleration)
              + dampingMatrix.solve(_newmarkGamma/(_newmarkBeta*_timestep)*currentNodalDisplacement + (_newmarkGamma/_newmarkBeta -1)*currentNodalVelocity + _timestep*(_newmarkGamma/(2*_newmarkBeta)-1)*currentNodalAcceleration);
    effectiveForceVector = ((1/(_newmarkBeta*_timestep*_timestep)*consistentMassMatrix)+((_newmarkGamma/(_newmarkBeta*_timestep))*dampingMatrix)).slove(nodalNewDisplacements)
                           + globalForceVector
                           - rVector;

    // TODO: Boundary conditions II - Force
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
      
      // ...
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      effectiveForceVector(bc._nodeId * DegreesOfFreedom + bc._coordinate) = 0.0;
      
    }

    // TODO: Evaluate the residual based on the effectiveForceVector incl. BCs
    double residue = 0.0;
    residue = effectiveForceVector.norm();
    if (verbose == true){
      printf("Initial residue = %9.3e\n", residue);
    }

    // While the residue > tolerance compute Newton-Raphson iterations
    unsigned int numberOfIterations = 0;
    
    MatrixXd effectiveTangentMatrix;
    
    while( (residue > tolerance) && (numberOfIterations < maxIterations) ) {
    
      // TODO: Set the efficient tangentMatrix
      // ...
      effectiveTangentMatrix = stiffnessMatrix + (consistentMassMatrix/(_newmarkBeta*_timestep*_timestep) + _newmarkGamma * dampingMatrix/(_newmarkBeta*_timestep));
      // Boundary condition III - Effective Tangent Matrix
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
        
        // ...
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        effectiveTangentMatrix.row(bc._nodeId * DegreesOfFreedom + bc._coordinate).fill(0.0);
        effectiveTangentMatrix(bc._nodeId * DegreesOfFreedom + bc._coordinate, bc._nodeId * DegreesOfFreedom + bc._coordinate) = 1.0;
      }
      
      // Update newDisplacement using the Newmark method update rule
      newDisplacement -= effectiveTangentMatrix.lu().solve(effectiveForceVector);
      
      // TODO :Boundary conditions IV - Solution
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
        
        // ...
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
        newDisplacement(dofIndex)    = bc._constraint;
      }

      // TODO: again convert newDisplacement into nodal form (i.e. update nodalNewDisplacements)
      
      // ...
      vector<ElementVector> nodalNewDisplacements
              =  Utilities::distributeGlobalVectorToLocalVectors<Assembler0>(newDisplacement);
        
      // TODO: Evaluate the new damping matrix
      // ...
      Eigen::MatrixXd consistentMassMatrix(numberOfDofs, numberOfDofs);
      consistentMassMatrix = _assembler0.assembleConsistentMassMatrix() + _assembler1.assembleConsistentMassMatrix() + _assembler2.assembleConsistentMassMatrix();
      Eigen::MatrixXd stiffnessMatrix(numberOfDofs, numberOfDofs);
      stiffnessMatrix = _assembler0.assembleStiffnessMatrix(nodalNewDisplacements) + _assembler1.assembleStiffnessMatrix(nodalNewDisplacements) + _assembler2.assembleStiffnessMatrix(nodalNewDisplacements);
      Eigen::MatrixXd dampingMatrix(numberOfDofs, numberOfDofs);
      dampingMatrix = dampingAlpha * consistentMassMatrix + dampingBeta * stiffnessMatrix;

      // TODO: evaluate the new effective force vector
      // ...
      VectorXd effectiveForceVector; // ...
      VectorXd globalForceVector;
      globalForceVector = _assembler0.assembleForceVector(nodalNewDisplacements) + _assembler1.assembleForceVector(nodalNewDisplacements) + _assembler2.assembleForceVector(nodalNewDisplacements);
      effectiveForceVector = ((1/(_newmarkBeta*_timestep*_timestep)*consistentMassMatrix)+((_newmarkGamma/(_newmarkBeta*_timestep))*dampingMatrix)).slove(nodalNewDisplacements)
                             + globalForceVector
                             - rVector;


      // TODO: Boundary conditions V - Force
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) { 
        
        // ...
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        effectiveForceVector(bc._nodeId * DegreesOfFreedom + bc._coordinate) = 0.0;
        
      }

      // TODO: evaluate the residual based on the norm of effectiveForceVector and divide it by the numberOfDOFs
      residue = effectiveForceVector.norm()/numberOfDOFs; // ...

      if (verbose == true) {
        printf("Newton Raphson iteration %4u, residue = %8.3e\n", numberOfIterations, residue);
      }

      numberOfIterations++;

    }
     
    // Error check
    if (numberOfIterations == maxIterations) {
      throwException("Newton Raphson solver could not converge "
                     "in %u iterations.\nTolerance: %e \nResidue: %e",
                     maxIterations,tolerance,residue);
    }
    
    // TODO: Update the states, i.e. save the new displacement, velocity and acceleration onto
    //       currentNodalDisplacement, currentNodalVelocity and currentNodalAcceleration
    
    // ...
    tempNodalVelocity = currentNodalVelocity;
    tempNodalAcceleration = currentNodalAcceleration;
    tempNodalDisplacement = currentNodalDisplacement;
    currentNodalAcceleration = 1/(_newmarkBeta*_timestep*_timestep)(nodalNewDisplacements - tempNodalDisplacement - _timestep*tempNodalVelocity)  - ((1-2*_newmarkBeta)/(2*_newmarkBeta))*tempNodalAcceleration;
    currentNodalVelocity = (1 - _newmarkGamma/_newmarkBeta)*tempNodalVelocity + (_newmarkGamma/(_newmarkBeta*_timestep))*(nodalNewDisplacements - tempNodalDisplacement) - _timestep*(_newmarkGamma/(2*_newmarkBeta) -1)*tempNodalAcceleration;
    currentNodalDisplacement = nodalNewDisplacements;
    
    // TODO: Delete the following when you're done
    //ignoreUnusedVariables(numberOfDOFs);

  }
  
  
  
  
  private:
  
  Assembler0 _assembler0;
  Assembler1 _assembler1;
  Assembler2 _assembler2;
  
  const double _timestep    ;
  const double _dampingAlpha;
  const double _dampingBeta ;
  const double _newmarkBeta ;
  const double _newmarkGamma;
  

};

#endif  // SOLVER_CLASS_H
