#ifndef __AGENT_HPP__
#define __AGENT_HPP__

#include <Eigen/Dense>
#include <iostream>
#include <math.h>

class Agent
{
public:
  // by boxes we get the one cycle of leg locomotion
  //(ex. in bipedal:move forward  left + right foot - which is two individual steps)
  // (ex. in quaduapedal 4 individual steps is 1 cycle - means 1 box)
  int world_dim, k_legs, contact_legs;
  double mass, g, height, mu;
  double horizon_time, dt;
  double wl, wk, wf, ws;

  Eigen::VectorXd target_xyz;
  Eigen::VectorXd starting_com_pos, starting_com_orient,
                      starting_com_vel, starting_com_ang_vel;

  /** default constructor */
  Agent();

  /** default destructor */
  virtual ~Agent();

  void setDim(int dim_);
  void setLegs(int legs_, int legs_cont_);
  void setDt(double dt_);

  void setHorizonTime(int time_h_);
  void setTargetCoM_xyz(const Eigen::Vector3d& target_xyz_);

  void storeStartingCoM_PosVel();

  int getLegs();




};

#endif
