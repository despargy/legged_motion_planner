#ifndef __LABELS_H__
#define __LABELS_H__

#include <Eigen/Dense>

using namespace std;

struct labels{

  Eigen::Vector3d CoM_positions;
  Eigen::Vector3d CoM_linear_velocities;
  Eigen::Vector3d CoM_linear_accelerations;

  vector<Eigen::Vector3d> foot_position_v;

  Eigen::Vector3d ZMP_positions;

  void set_foot_size(int k_legs){
    foot_position_v.resize(k_legs);
  };

};
#endif
