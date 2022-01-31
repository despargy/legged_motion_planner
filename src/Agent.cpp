#include "legged_motion_planner/Agent.hpp"

Agent::Agent()
{
  printf("This is from 'AGENT' constructor\n");

  g = 9.8; //m/s^2
  mass = 72.0; //kg
  height = 1.20952; //m

   wl =1; wk=2; wf=2; ws =2; //TODO fix weights
   mu = 0.8;

 }

Agent::~Agent()
{ }

void Agent::setDim(int dim_){

  world_dim = dim_;
}

void Agent::setLegs(int legs_){

  k_legs = legs_;
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

  starting_com_vel.resize(3);
  starting_com_vel.setZero();

  starting_com_orient.resize(3);
  starting_com_orient.setZero();

  starting_com_ang_vel.resize(3);
  starting_com_ang_vel.setZero();
}


int Agent::getLegs(){

  return k_legs;
}
