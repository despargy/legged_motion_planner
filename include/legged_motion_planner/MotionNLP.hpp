

#ifndef __MOTIONNLP_HPP__
#define __MOTIONNLP_HPP__

#include "IpTNLP.hpp"
#include <math.h>
#include <Eigen/Dense>
#include <iostream>

#include "motion_planner/Agent.hpp"
#include <motion_planner/NumericalIntegrationConstraints.hpp>

using namespace Ipopt;

class MotionNLP: public TNLP
{


public:

  Agent* p_agent;

  /** default constructor */
  MotionNLP();
  MotionNLP(Agent* p_agent);

  /** default destructor */
  virtual ~MotionNLP();

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
     MotionNLP(
        const MotionNLP&
     );

     MotionNLP& operator=(
        const MotionNLP&
     );
     //@}
     NumericalIntegrationConstraints integration_constr;

  };

  #endif
