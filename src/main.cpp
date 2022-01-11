#include "legged_motion_planner/CentroidalNLP.hpp"
#include "legged_motion_planner/Agent.hpp"

#include "IpIpoptApplication.hpp"
#include <Eigen/Dense>
#include <iostream>
#include "legged_motion_planner/class_csv.h"
#include "ros/ros.h"



int main(int argc, char** argv)
{
  bool SOLVE_IPOPT = false;
  // read data
  vector<labels> desired = pass_csv("/home/atom/master_ws/src/legged_motion_planner/records/DesiredTrajectory3.csv");

  // 3D world
  int WORLD_DIM = 3;
  //define legs
  int NUM_LEGS = 2;
  double IND_STEP_TIME = 2.0 ; // double sup. -> 0.5sec + single sup. -> 1.5sec
  // dt
  double dt = 0.01; // in sec
  // frequency
  double freq = 1/dt;
  // points to solved
  int POINTS = desired.size(); // all data - same as desired
  std::cout<<POINTS<<std::endl;

  // time to achieve goal
  double HORIZON_TIME = POINTS*dt;
  // define boxes
  int BOXES = (int)HORIZON_TIME/(IND_STEP_TIME*NUM_LEGS); // total time div. (time for each leg * #k_legs) -

  //define goal for CoM
  Eigen::Vector3d TARGET_XYZ;
  TARGET_XYZ = desired[POINTS-1].CoM_positions;

  Agent *p_agent = new Agent();

  p_agent->setDim(WORLD_DIM);
  p_agent->setLegs(NUM_LEGS);
  p_agent->setDt(dt);
  p_agent->setOptimizationBoxes(BOXES);
  p_agent->setHorizonTime(HORIZON_TIME);
  // while
  p_agent->setTargetCoM_xyz(TARGET_XYZ);
  // p_agent->setDesAvVel();
  p_agent->storeStartingCoM_PosVel(); //TODO
  //
  // // 2 cycles, 5 phase duration for each foot, 4 itnegration step = 20 phase duration
  // //TODO
  p_agent->defineCoMConstraint(); //TODO
  // p_agent->defineFeetConstraint(); //TODO
  // p_agent->defineContactForcesConstraint(); //TODO
  // p_agent->definePhaseConstraint(); //TODO

  std::cout<<"pass"<<std::endl;

  if(SOLVE_IPOPT){

    /* set problem - IPOPT throw IpoptApplication */

    // Create a new instance of our nlp
    // Problem for this Agent, pass pointer

    SmartPtr<CentroidalNLP> mynlp  = new CentroidalNLP(p_agent);

    mynlp->dt = dt;
    mynlp->n_points = POINTS;

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

  }

  return 0;
}
