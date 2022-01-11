#include "legged_motion_planner/CentroidalNLP.hpp"

#include <cassert>

using namespace Ipopt;

/* Constructor. */
CentroidalNLP::CentroidalNLP()
{ printf("This is from 'MOTION NLP' constructor\n"); }

CentroidalNLP::CentroidalNLP(Agent* p_agent_)
{
  printf("This is from 'MOTION NLP' constructor\n");
  p_agent = p_agent_;

 }

CentroidalNLP::~CentroidalNLP()
{ }

/*
  CoM position - 3var = 3D
  CoM orientation - 3var = 3D
  feet position - 3*k_legs = 3D*k_legs
  contact forces - 3*k_legs = 3D*k_legs
  phase duration -
*/

bool CentroidalNLP::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
   // The problem described in CentroidalNLP.hpp has 2 variables, x, & u,
   n = p_agent->variables_per_box*(p_agent->boxes); //TODO without -1?

   // equality constraint,
   m = p_agent->constraints_per_box*(p_agent->boxes-1);

   // 2 nonzeros in the jacobian (one for x1, and one for x2),
   nnz_jac_g = p_agent->nnz_jac_g_per_box*(p_agent->boxes-1);

   // and 2 nonzeros in the hessian of the lagrangian
   // (one in the hessian of the objective for x2,
   //  and one in the hessian of the constraints for x1)
   nnz_h_lag = 0;

   // We use the standard fortran index style for row/col entries
   index_style = C_STYLE;//FORTRAN_STYLE;

   return true;
}

bool CentroidalNLP::get_bounds_info(
   Index   n,
   Number* x_l,
   Number* x_u,
   Index   m,
   Number* g_l,
   Number* g_u
)
{

   // CoM Position
   x_l[p_agent->CoM_pos_index + 0] = p_agent->starting_com_pos(0); //x
   x_u[p_agent->CoM_pos_index +0] = p_agent->starting_com_pos(0);

   x_l[p_agent->CoM_pos_index +1] = p_agent->starting_com_pos(1); //y
   x_u[p_agent->CoM_pos_index +1] = p_agent->starting_com_pos(1);

   x_l[p_agent->CoM_pos_index +2] = p_agent->starting_com_pos(2); //z
   x_u[p_agent->CoM_pos_index +2] = p_agent->starting_com_pos(2);

   // CoM Velocity SOS mazi ta pira
   x_l[p_agent->CoM_pos_index + 3] = p_agent->starting_com_vel(0); //x
   x_u[p_agent->CoM_pos_index +3] = p_agent->starting_com_vel(0);

   x_l[p_agent->CoM_pos_index +4] = p_agent->starting_com_vel(1); //y
   x_u[p_agent->CoM_pos_index +4] = p_agent->starting_com_vel(1);

   x_l[p_agent->CoM_pos_index +5] = p_agent->starting_com_vel(2); //z
   x_u[p_agent->CoM_pos_index +5] = p_agent->starting_com_vel(2);

   for(int b = 1; b<p_agent->boxes; b++){
     //CoM Position in the world frame:

     for(int j = p_agent->CoM_pos_index; j<p_agent->variables_per_box; j++){
       x_l[b * p_agent->variables_per_box + j] = -100.0; x_u[b * p_agent->variables_per_box + j] = +100.0;
     }
   }

   for(int k =  p_agent->CoM_pos_index; k< (p_agent->boxes-1); k++){
     // CoM possition has 3 constraints_per_box
     for(int c= 0; c<3; c++){
       g_l[k * p_agent->constraints_per_box + c] = 0.0; g_u[k * p_agent->constraints_per_box + c] = 0.0;
     }
   }
   //THE SAME FOR OTHER VARS
   // // CoM Velocity
   // x_l[3] = p_agent->starting_com_vel(0);
   // x_u[3] = p_agent->starting_com_vel(0);
   //
   // x_l[4] = p_agent->starting_com_vel(1);
   // x_u[4] = p_agent->starting_com_vel(1);
   //
   // x_l[5] = p_agent->starting_com_vel(2);
   // x_u[5] = p_agent->starting_com_vel(2);

   // // x2 has no upper or lower bound, so we set them to
   // // a large negative and a large positive number.
   // // The value that is interpretted as -/+infinity can be
   // // set in the options, but it defaults to -/+1e19
   // x_l[1] = -1.0e19;
   // x_u[1] = +1.0e19;

   // we have one equality constraint, so we set the bounds on this constraint
   // to be equal (and zero).
   // g_l[0] = g_u[0] = 0.0;

   return true;
}

bool CentroidalNLP::get_starting_point(
   Index   n,
   bool    init_x,
   Number* x,
   bool    init_z,
   Number* z_L,
   Number* z_U,
   Index   m,
   bool    init_lambda,
   Number* lambda
)
{
   // Here, we assume we only have starting values for x, if you code
   // your own NLP, you can provide starting values for the others if
   // you wish.
   assert(init_x == true);
   assert(init_z == false);
   assert(init_lambda == false);

   //TODO why for each box??????
   for(int k = 0; k<p_agent->boxes; k++){

     for(int j = 0; j<p_agent->variables_per_box; j++){
       x[k * p_agent->variables_per_box + j] = 0.0;
     }

   }

   return true;
}

bool CentroidalNLP::eval_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number&       obj_value
)
{
  // return the value of the objective function
  const Eigen::Map<const Eigen::VectorXd> decision_vars(x, n);

  Eigen::Vector3d final_com_state, com_lin_vel, grad_f;
  final_com_state(0) = decision_vars((p_agent->boxes-1)*6+0);
  final_com_state(1) = decision_vars((p_agent->boxes-1)*6+1);
  final_com_state(2) = decision_vars((p_agent->boxes-1)*6+2);
  double cost_temp = 0.0;
  for(int k = 0; k<p_agent->boxes; k++){
    com_lin_vel(0) = decision_vars(k * p_agent->variables_per_box + 3);
    com_lin_vel(1) = decision_vars(k * p_agent->variables_per_box + 4);
    com_lin_vel(2) = decision_vars(k * p_agent->variables_per_box + 5);
    cost_temp += getQuadraticCost(com_lin_vel, p_agent->desired_com_lin_vel, grad_f);
  }
  obj_value = cost_temp + getQuadraticCost(final_com_state, p_agent->target_xyz, grad_f);

  return true;
}

bool CentroidalNLP::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
   // return the gradient of the objective function grad_{x} f(x)
   const Eigen::Map<const Eigen::VectorXd> decision_vars(x, n);
   Eigen::Vector3d final_com_state;
   final_com_state(0) = decision_vars(n-6);
   final_com_state(1) = decision_vars(n-5);
   final_com_state(2) = decision_vars(n-4);

   Eigen::Vector3d com_cost_gradient;
   getQuadraticCost(final_com_state, p_agent->target_xyz, com_cost_gradient);

   for(int k = 0; k<p_agent->boxes; k++){
     Eigen::Vector3d com_lin_vel, com_lin_vel_cost_gradient;
     com_lin_vel(0) = decision_vars(k * p_agent->variables_per_box + 3);
     com_lin_vel(1) = decision_vars(k * p_agent->variables_per_box + 4);
     com_lin_vel(2) = decision_vars(k * p_agent->variables_per_box + 5);
     getQuadraticCost(com_lin_vel, p_agent->desired_com_lin_vel, com_lin_vel_cost_gradient);

     if(k == p_agent->boxes-1){
         grad_f[k * p_agent->variables_per_box + 0] = com_cost_gradient(0);
         grad_f[k * p_agent->variables_per_box + 1] = com_cost_gradient(1);
         grad_f[k * p_agent->variables_per_box + 2] = com_cost_gradient(2);
         grad_f[k * p_agent->variables_per_box + 3] = com_lin_vel_cost_gradient(0);
         grad_f[k * p_agent->variables_per_box + 4] = com_lin_vel_cost_gradient(1);
         grad_f[k * p_agent->variables_per_box + 5] = com_lin_vel_cost_gradient(2);
      }else{
         grad_f[k * p_agent->variables_per_box + 0] = 0.0;
         grad_f[k * p_agent->variables_per_box + 1] = 0.0;
         grad_f[k * p_agent->variables_per_box + 2] = 0.0;
         grad_f[k * p_agent->variables_per_box + 3] = com_lin_vel_cost_gradient(0);
         grad_f[k * p_agent->variables_per_box + 4] = com_lin_vel_cost_gradient(1);
         grad_f[k * p_agent->variables_per_box + 5] = com_lin_vel_cost_gradient(2);
     }
   }
   return true;
}

bool CentroidalNLP::eval_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Number*       g
)
{
   // return the value of the constraints: g(x)
   // return the value of the constraints: g(x)
 const Eigen::Map<const Eigen::VectorXd> decision_vars(x, n);
 const Eigen::Map<const Eigen::VectorXd> constraints(g, m);

 Eigen::Vector3d next_com_pos, com_pos, com_lin_vel, constr_violation;
 Eigen::VectorXd delta_t(p_agent->boxes);
 delta_t.setZero();

 for(int k = 0; k<p_agent->boxes-1; k++){
   double kd = (double)k;
   double incr = p_agent->horizon_time / (double)p_agent->boxes;
   delta_t(k) = incr;
 }

 for(int k = 0; k<p_agent->boxes-1; k++){
   com_pos(0) = decision_vars(k * p_agent->variables_per_box + 0);
   com_pos(1) = decision_vars(k * p_agent->variables_per_box + 1);
   com_pos(2) = decision_vars(k * p_agent->variables_per_box + 2);
   com_lin_vel(0) = decision_vars(k * p_agent->variables_per_box + 3);
   com_lin_vel(1) = decision_vars(k * p_agent->variables_per_box + 4);
   com_lin_vel(2) = decision_vars(k * p_agent->variables_per_box + 5);
   next_com_pos(0) = decision_vars((k+1) * p_agent->variables_per_box + 0);
   next_com_pos(1) = decision_vars((k+1) * p_agent->variables_per_box + 1);
   next_com_pos(2) = decision_vars((k+1) * p_agent->variables_per_box + 2);
   integration_constr.firstOrderMethod(delta_t(k), com_pos, com_lin_vel, next_com_pos, constr_violation);
   g[k * p_agent->constraints_per_box + 0] = constr_violation(0);
   g[k * p_agent->constraints_per_box + 1] = constr_violation(1);
   g[k * p_agent->constraints_per_box + 2] = constr_violation(2);
 }

  return true;
}

bool CentroidalNLP::eval_jac_g(
   Index         n,
   const Number* x,
   bool          new_x,
   Index         m,
   Index         nele_jac,
   Index*        iRow,
   Index*        jCol,
   Number*       values
)
{

  Eigen::VectorXd delta_t(p_agent->boxes);
  delta_t.setZero();
  for(int k = 0; k < p_agent->boxes-1; k++){
    double kd = (double)k;
    double incr = p_agent->horizon_time / (double)p_agent->boxes;
    delta_t(k) = kd * incr + 1.0;
  }
  //std::cout<<delta_t.transpose()<<std::endl;

    for(int k = 0; k<p_agent->boxes-1; k++){
      if (values == NULL) {
        // return the structure of the Jacobian

        // this particular Jacobian is dense
        iRow[k*p_agent->nnz_jac_g_per_box + 0] = 0 + k * 3; jCol[k*p_agent->nnz_jac_g_per_box + 0] = 0 + k * 6;
        iRow[k*p_agent->nnz_jac_g_per_box + 1] = 0 + k * 3; jCol[k*p_agent->nnz_jac_g_per_box + 1] = 3 + k * 6;
        iRow[k*p_agent->nnz_jac_g_per_box + 2] = 0 + k * 3; jCol[k*p_agent->nnz_jac_g_per_box + 2] = 6 + k * 6;

        iRow[k*p_agent->nnz_jac_g_per_box + 3] = 1 + k * 3; jCol[k*p_agent->nnz_jac_g_per_box + 3] = 1 + k * 6;
        iRow[k*p_agent->nnz_jac_g_per_box + 4] = 1 + k * 3; jCol[k*p_agent->nnz_jac_g_per_box + 4] = 4 + k * 6;
        iRow[k*p_agent->nnz_jac_g_per_box + 5] = 1 + k * 3; jCol[k*p_agent->nnz_jac_g_per_box + 5] = 7 + k * 6;

        iRow[k*p_agent->nnz_jac_g_per_box + 6] = 2 + k * 3; jCol[k*p_agent->nnz_jac_g_per_box + 6] = 2 + k * 6;
        iRow[k*p_agent->nnz_jac_g_per_box + 7] = 2 + k * 3; jCol[k*p_agent->nnz_jac_g_per_box + 7] = 5 + k * 6;
        iRow[k*p_agent->nnz_jac_g_per_box + 8] = 2 + k * 3; jCol[k*p_agent->nnz_jac_g_per_box + 8] = 8 + k * 6;
      }
      else {
        // return the values of the Jacobian of the constraints

        values[k*p_agent->nnz_jac_g_per_box + 0] = -1.0;
        values[k*p_agent->nnz_jac_g_per_box + 1] = -delta_t(k);
        values[k*p_agent->nnz_jac_g_per_box + 2] = 1.0;

        values[k*p_agent->nnz_jac_g_per_box + 3] = -1.0;
        values[k*p_agent->nnz_jac_g_per_box + 4] = -delta_t(k);
        values[k*p_agent->nnz_jac_g_per_box + 5] = 1.0;

        values[k*p_agent->nnz_jac_g_per_box + 6] = -1.0;
        values[k*p_agent->nnz_jac_g_per_box + 7] = -delta_t(k);
        values[k*p_agent->nnz_jac_g_per_box + 8] = 1.0;
      }
    }

   return true;
}

bool CentroidalNLP::eval_h(
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
)
{

   return false;
}

void CentroidalNLP::finalize_solution(
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
)
{
   // here is where we would store the solution to variables, or write to a file, etc
   // so we could use the solution. Since the solution is displayed to the console,
   // we currently do nothing here.
   printf("\n\nObjective value\n");
   printf("f(x*) = %e\n", obj_value);
}

double CentroidalNLP::getQuadraticCost(const Eigen::Vector3d & final_state,
                                                    const Eigen::Vector3d & target,
                                                    Eigen::Vector3d & gradient_f)
{

  gradient_f(0) = 2.0*(final_state(0)-target(0));
  gradient_f(1) = 2.0*(final_state(1)-target(1));
  gradient_f(2) = 2.0*(final_state(2)-target(2));

  return pow((final_state(0)-target(0)),2.0) +
        pow((final_state(1)-target(1)),2.0) +
        pow((final_state(2)-target(2)),2.0);
}
