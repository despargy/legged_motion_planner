#include "legged_motion_planner/CentroidalNLP.hpp"
#include "legged_motion_planner/Agent.hpp"

#include "IpIpoptApplication.hpp"
#include <Eigen/Dense>
#include <iostream>
#include "legged_motion_planner/class_csv.h"
#include "ros/ros.h"



int main(int argc, char** argv)
{
  bool SOLVE_IPOPT = true;

  // 3D world
  int WORLD_DIM = 3;
  //define legs
  int NUM_LEGS = 2;
  // dt
  double dt = 0.01; // in sec
  // frequency
  double freq = 1/dt;
  // read data
  // vector<labels> desired = pass_csv("/home/despargy/master_g15_ws/src/legged_motion_planner/records/DesiredTrajectory01.csv", NUM_LEGS);
  vector<labels> desired = pass_csv("/home/atom/master_ws/src/legged_motion_planner/records/DesiredTrajectory01.csv", NUM_LEGS);

  // points to solved
  int POINTS = desired.size(); // all data - same as desired
  std::cout<<POINTS<<std::endl;

  double IND_STEP_TIME = 2.0 ; // double sup. -> 0.5sec + single sup. -> 1.5sec

  // time to achieve goal
  double HORIZON_TIME = POINTS*dt;

  //define goal for CoM
  Eigen::Vector3d TARGET_XYZ;
  TARGET_XYZ = desired[POINTS-1].CoM_positions;

  Agent *p_agent = new Agent();

  p_agent->setDim(WORLD_DIM);
  p_agent->setLegs(NUM_LEGS);
  p_agent->setDt(dt);

  p_agent->setHorizonTime(HORIZON_TIME);

  p_agent->setTargetCoM_xyz(TARGET_XYZ);
  // p_agent->storeStartingCoM_PosVel(); //TODO


  std::cout<<"pass"<<std::endl;

  if(SOLVE_IPOPT){

    /* set problem - IPOPT throw IpoptApplication */

    // Create a new instance of our nlp
    // Problem for this Agent, pass pointer

    SmartPtr<CentroidalNLP> mynlp  = new CentroidalNLP(p_agent);

    mynlp->dt = dt;
    mynlp->n_points = POINTS-1; // USE THE FIRST AS INIT POINT - NOT SOLVE ANYTHING FOR THAT
    mynlp->desired = desired;
    bool init_state = mynlp->init_before_start();

    std::cout<<init_state<<std::endl;

    // Create an instance of the IpoptApplication
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    /* define parameters */
    app->RethrowNonIpoptException(true);

    app->Options()->SetNumericValue("tol", 1e-6);
    app->Options()->SetNumericValue("mu_init", 1e-2);
    app->Options()->SetNumericValue("max_iter", 3000);
    app->Options()->SetNumericValue("print_level", 5);

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
