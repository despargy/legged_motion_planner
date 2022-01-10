#include "motion_planner/MotionNLP.hpp"
#include "motion_planner/Agent.hpp"

#include "IpIpoptApplication.hpp"
#include <Eigen/Dense>
#include <iostream>

#include "ros/ros.h"



int main(int argc, char** argv)
{

  // 3D world
  int WORLD_DIM = 3;
  //define legs
  int NUM_LEGS = 4;
  // define boxes
  int BOXES = 10; // #TODO
  // time to achieve goal
  double HORIZON_TIME = 6.0;
  //define goal for CoM
  Eigen::Vector3d TARGET_XYZ;
  TARGET_XYZ = Eigen::Vector3d(10.0,5.0,1.0);

  Agent *p_agent = new Agent();

  p_agent->setDim(WORLD_DIM);
  p_agent->setLegs(NUM_LEGS);
  p_agent->setOptimizationBoxes(BOXES);
  p_agent->setHorizonTime(HORIZON_TIME);
  // while
  p_agent->setTargetCoM_xyz(TARGET_XYZ);
  p_agent->setDesAvVel();
  p_agent->storeStartingCoM_PosVel(); //TODO
  //
  // // 2 cycles, 5 phase duration for each foot, 4 itnegration step = 20 phase duration
  // //TODO
  p_agent->defineCoMConstraint(); //TODO
  // p_agent->defineFeetConstraint(); //TODO
  // p_agent->defineContactForcesConstraint(); //TODO
  // p_agent->definePhaseConstraint(); //TODO

  std::cout<<"pass"<<std::endl;

  /* set problem - IPOPT throw IpoptApplication */

  // Create a new instance of our nlp
  // Problem for this Agent, pass pointer

  SmartPtr<MotionNLP> mynlp  = new MotionNLP(p_agent);
  std::cout<<mynlp->p_agent->getLegs()<<std::endl;

  // Create an instance of the IpoptApplication
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  /* define parameters */
  app->RethrowNonIpoptException(true);

  app->Options()->SetNumericValue("tol", 1e-6);
  app->Options()->SetNumericValue("mu_init", 1e-2);
  // app->Options()->SetNumericValue("max_iter", 3000);
  // app->Options()->SetNumericValue("print_level", 5);

  app->Options()->SetStringValue("nlp_scaling_method", "none");
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");
  app->Options()->SetStringValue("hessian_approximation", "limited-memory");

  // // #TODO what they are doing????????????????
  // ////////////////////////////////////////////// TO CHANGE ///////////////
  // //
  // // mynlp->addCoMOrientConstraint();
  // // mynlp->addCoMPosConstraint();
  // ////////////////////////////////////////////// END ///////////////

  // Initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if( status != Solve_Succeeded )
  {
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return (int) status;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);

  if (status == Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
  }
  else {
    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
  }



  return 0;
}
