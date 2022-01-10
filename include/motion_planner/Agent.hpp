#ifndef __AGENT_HPP__
#define __AGENT_HPP__

#include <Eigen/Dense>
#include <iostream>
#include <math.h>

class Agent
{
public:

  int world_dim, k_legs, boxes,
    variables_per_box, constraints_per_box, nnz_jac_g_per_box;
  double mass;
  double horizon_time;

  int next_index, CoM_pos_index;

  Eigen::Vector3d desired_com_lin_vel;

  Eigen::VectorXd target_xyz;
  Eigen::VectorXd starting_com_pos, starting_com_orient,
                      starting_com_vel, starting_com_ang_vel;

  /** default constructor */
  Agent();

  /** default destructor */
  virtual ~Agent();

  void setDim(int dim_);
  void setLegs(int legs_);
  void setOptimizationBoxes(int boxes_);
  void setHorizonTime(int time_h_);
  void setTargetCoM_xyz(const Eigen::Vector3d& target_xyz_);
  void setDesAvVel();
  void storeStartingCoM_PosVel();
  void defineCoMConstraint();
  void defineFeetConstraint();
  void defineContactForcesConstraint();
  void definePhaseConstraint();

  int getLegs();
  ////////////////////////////////////////////// TO CHANGE ///////////////
  //
  // void addCoMPosConstraint();
  // void addCoMOrientConstraint();
  ////////////////////////////////////////////// END ///////////////



};

#endif
