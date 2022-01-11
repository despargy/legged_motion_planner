#include "legged_motion_planner/Agent.hpp"

Agent::Agent()
{
  printf("This is from 'AGENT' constructor\n");
  // init variables
  variables_per_box = 0;
  constraints_per_box = 0;
  nnz_jac_g_per_box = 0 ;

  next_index = 0;

 }

Agent::~Agent()
{ }

void Agent::setDim(int dim_){

  world_dim = dim_;
}

void Agent::setLegs(int legs_){

  k_legs = legs_;
}

void Agent::setOptimizationBoxes(int boxes_){

  boxes = boxes_;
}

void Agent::setDt(double dt_){

  dt = dt_;
}

void Agent::setHorizonTime(int time_h_){

  horizon_time = time_h_;
}

void Agent::setTargetCoM_xyz(const Eigen::Vector3d& target_xyz_){

  target_xyz = target_xyz_;
}

void Agent::storeStartingCoM_PosVel(){

  starting_com_pos.resize(3);
  starting_com_pos.setZero();
  // starting_com_pos = Eigen::Vector3d(0.0,0.0,0.0);
  // starting_com_pos << 1.0, 2.0, 3.0;

  // std::cout<<starting_com_pos<<std::endl;


  starting_com_vel.resize(3);
  starting_com_vel.setZero();

  starting_com_orient.resize(3);
  starting_com_orient.setZero();

  starting_com_ang_vel.resize(3);
  starting_com_ang_vel.setZero();
}

// void Agent::setDesAvVel(){
//
//   desired_com_lin_vel(0) = target_xyz(0)/((double)boxes);
//   desired_com_lin_vel(1) = target_xyz(1)/((double)boxes);
//   desired_com_lin_vel(2) = target_xyz(2)/((double)boxes);
//
// }

int Agent::getLegs(){

  return k_legs;
}

void Agent::defineCoMConstraint(){

  // position
  variables_per_box += 6; //(xyz)
  constraints_per_box += 3;
  nnz_jac_g_per_box += 9;

  CoM_pos_index = next_index;
  next_index += 6*boxes; //TODO -1 sta boxes?

  // //orientation
  // variables_per_box += 6; //(xyz, rpy)
  // constraints_per_box += 3;
  // nnz_jac_g_per_box += 9;
  //
  // // CoM_ori_index = next_index;
  // // next_index += variables_per_box*boxes;
}


void Agent::defineFeetConstraint(){

  // k_legs
  variables_per_box += k_legs;
  constraints_per_box += k_legs;
  nnz_jac_g_per_box += 9;

  // Feet_index = next_index;
  // next_index += variables_per_box*boxes;
}

void Agent::defineContactForcesConstraint(){

  // k_legs
  variables_per_box += k_legs;
  constraints_per_box += k_legs;
  nnz_jac_g_per_box += 9;

  // ContactForces_index = next_index;
  // next_index += variables_per_box*boxes;
}

void Agent::definePhaseConstraint(){

  // k_legs
  variables_per_box += k_legs;
  constraints_per_box += k_legs;
  nnz_jac_g_per_box += 9;

  // Phase_index = next_index;
  // next_index += variables_per_box*boxes;
}
