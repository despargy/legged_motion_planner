

#ifndef __CENTROIDALNNLP_HPP__
#define __CENTROIDALNLP_HPP__

#include "IpTNLP.hpp"
#include <math.h>
// #include <Eigen/Dense>
#include <iostream>
#include "legged_motion_planner/Agent.hpp"
#include <legged_motion_planner/labels.h>

using namespace Ipopt;

class CentroidalNLP: public TNLP
{


public:

  Agent* p_agent;
  int n_points, n_per_dt, g_per_dt, v_per_dt, nnz_jac_g_per_dt;
  double dt;
  int index_CoM_pos, index_CoM_lin, index_CoM_ang, index_Forces, index_Step, index_ends;
  int index_conF_eq_1, index_conF_eq_2, index_conF_eq_3, index_conF_eq_4, index_conF_eq_5, index_conF_eq_6, index_conF_eq_7;
  int index_partial_eq_1, index_partial_eq_2, index_partial_eq_3, index_partial_eq_4, index_partial_eq_5, index_partial_eq_6, index_partial_eq_7;
  vector<labels> desired;
  /** default constructor */
  CentroidalNLP();
  CentroidalNLP(Agent* p_agent);
  bool init_before_start();

  /** default destructor */
  virtual ~CentroidalNLP();

     /**@name Overloaded from TNLP */
     //@{
     /** Method to return some info about the nlp */
     virtual bool get_nlp_info(
        Index&          n,
        Index&          m,
        Index&          nnz_jac_g,
        Index&          nnz_h_lag,
        IndexStyleEnum& index_style
     );

     /** Method to return the bounds for my problem */
     virtual bool get_bounds_info(
        Index   n,
        Number* x_l,
        Number* x_u,
        Index   m,
        Number* g_l,
        Number* g_u
     );

     /** Method to return the starting point for the algorithm */
     virtual bool get_starting_point(
        Index   n,
        bool    init_x,
        Number* x,
        bool    init_z,
        Number* z_L,
        Number* z_U,
        Index   m,
        bool    init_lambda,
        Number* lambda
     );

     /** Method to return the objective value */
     virtual bool eval_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number&       obj_value
     );

     /** Method to return the gradient of the objective */
     virtual bool eval_grad_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number*       grad_f
     );

     /** Method to return the constraint residuals */
     virtual bool eval_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Number*       g
     );

     /** Method to return:
      *   1) The structure of the Jacobian (if "values" is NULL)
      *   2) The values of the Jacobian (if "values" is not NULL)
      */
     virtual bool eval_jac_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Index         nele_jac,
        Index*        iRow,
        Index*        jCol,
        Number*       values
     );

     /** Method to return:
      *   1) The structure of the Hessian of the Lagrangian (if "values" is NULL)
      *   2) The values of the Hessian of the Lagrangian (if "values" is not NULL)
      */
     virtual bool eval_h(
        Index         n,
        const Number* x,
        bool          new_x,
        Number        obj_factor,
        Index         m,
        const Number* lambda,
        bool          new_lambda,
        Index         nele_hess,
        Index*        iRow,
        Index*        jCol,
        Number*       values
     );

     /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
     virtual void finalize_solution(
        SolverReturn               status,
        Index                      n,
        const Number*              x,
        const Number*              z_L,
        const Number*              z_U,
        Index                      m,
        const Number*              g,
        const Number*              lambda,
        Number                     obj_value,
        const IpoptData*           ip_data,
        IpoptCalculatedQuantities* ip_cq
     );
     //@}
     double getQuadraticCost(const Eigen::Vector3d & final_state,
                                                         const Eigen::Vector3d & target,
                                                         Eigen::Vector3d & gradient_f);

     Eigen::Vector3d crossProduct(Eigen::Vector3d vect_A, Eigen::Vector3d vect_B);

  private:
     /**@name Methods to block default compiler methods.
      *
      * The compiler automatically generates the following three methods.
      *  Since the default compiler implementation is generally not what
      *  you want (for all but the most simple classes), we usually
      *  put the declarations of these methods in the private section
      *  and never implement them. This prevents the compiler from
      *  implementing an incorrect "default" behavior without us
      *  knowing. (See Scott Meyers book, "Effective C++")
      */
     //@{
     CentroidalNLP(
        const CentroidalNLP&
     );

     CentroidalNLP& operator=(
        const CentroidalNLP&
     );

  };

  #endif
