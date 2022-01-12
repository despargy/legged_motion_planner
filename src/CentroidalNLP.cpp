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

bool CentroidalNLP::init_before_start(
){
  // don't foret to init indexs to help you map x variables with CoM, l, k, Forces, Steps
  index_CoM_pos = 0; // it has 3 so the next index is 0+3
  index_CoM_lin = index_CoM_pos + 3; // it has 3 so the next index is 3+3
  index_CoM_ang = index_CoM_lin + 3; // it has 3 so the next index is 6+3
  index_Forces = index_CoM_ang + 3; // it has 3 but MANY legs(+3*p_agent->k_legs)
  index_Step = index_Forces + 3*p_agent->k_legs;
  index_ends = index_Step + 3*p_agent->k_legs;

  index_conF_eq_1 = 0;
  index_conF_eq_2 = index_conF_eq_1 + p_agent->k_legs;
  index_conF_eq_3 = index_conF_eq_2 + p_agent->k_legs;
  index_conF_eq_4 = index_conF_eq_3 + p_agent->k_legs;
  index_conF_eq_5 = index_conF_eq_4 + p_agent->k_legs;
  index_conF_eq_6 = index_conF_eq_5 + 3;
  index_conF_eq_7 = index_conF_eq_5 + 3;

}


bool CentroidalNLP::get_nlp_info(
   Index&          n,
   Index&          m,
   Index&          nnz_jac_g,
   Index&          nnz_h_lag,
   IndexStyleEnum& index_style
)
{
   // The problem described in CentroidalNLP.hpp has for each dt:
   // CoM -> 3, l ->3 , k -> 3, Forces*legs -> 3*legs, Stepposition*legs -> 3*legs
   // index_CoM_pos, index_CoM_lin, index_CoM_ang, index_Forces, index_Step;

   n_per_dt = index_ends;
   n = n_per_dt*n_points; //TODO without -1?

   // equality constraints and inequality constraints
   // Fx*legs +- -> 2, Fy*legs +- -> 2, l->3, k->3, CoM->3 BUT for each dt (multi. with dt)
   m = 13*n_points; //TODO

   // 2 nonzeros in the jacobian (one for x1, and one for x2),
   nnz_jac_g = 9*n_points; //TODO;

   // and 2 nonzeros in the hessian of the lagrangian
   // (one in the hessian of the objective for x2,
   //  and one in the hessian of the constraints for x1)
   nnz_h_lag = 0; //TODO;

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

  // set all as inf, lets say by default
  for (int i=0; i < n; i++)
  {
    x_l[i] = -1.0e19;
    x_u[i] = +1.0e19;
  }

  //CoM z bounds ground - height
  for(int pp =0; pp<n_points; pp++)
  {
    x_l[pp+index_CoM_pos+2] = 0.0 ; x_u[pp + index_CoM_pos+2] = p_agent->height;
  }

  // cons bounds

  // for each point
  for(int pp =0; pp<n_points; pp++)
  {
    // eq1 0<= μFz - Fx
    for(int ll=0; ll<p_agent->k_legs; ll++)
    {
      g_l[pp + index_conF_eq_1 + ll ] = 0.0; g_u[pp + index_conF_eq_1 + ll ] = +1.0e19;
    }

    // eq2 μFz + Fx >= 0
    for(int ll=0; ll<p_agent->k_legs; ll++)
    {
      g_l[pp + index_conF_eq_2 + ll ] = 0.0; g_u[pp + index_conF_eq_2 + ll ] = +1.0e19;
    }

    // eq3 0<= μFz - Fy
    for(int ll=0; ll<p_agent->k_legs; ll++)
    {
      g_l[pp + index_conF_eq_3 + ll ] = 0.0; g_u[pp + index_conF_eq_3 + ll ] = +1.0e19;
    }

    // eq4 μFz + Fy >= 0
    for(int ll=0; ll<p_agent->k_legs; ll++)
    {
      g_l[pp + index_conF_eq_4 + ll ] = 0.0; g_u[pp + index_conF_eq_4 + ll ] = +1.0e19;
    }

    // eq5 but x y z l
    g_l[pp + index_conF_eq_5 + 0 ] = 0.0; g_u[pp + index_conF_eq_5 + 0 ] = 0.0;
    g_l[pp + index_conF_eq_5 + 1 ] = 0.0; g_u[pp + index_conF_eq_5 + 1 ] = 0.0;
    g_l[pp + index_conF_eq_5 + 2 ] = 0.0; g_u[pp + index_conF_eq_5 + 2 ] = 0.0;

    // eq6 but x y z k
    g_l[pp + index_conF_eq_6 + 0 ] = 0.0; g_u[pp + index_conF_eq_6 + 0 ] = 0.0;
    g_l[pp + index_conF_eq_6 + 1 ] = 0.0; g_u[pp + index_conF_eq_6 + 1 ] = 0.0;
    g_l[pp + index_conF_eq_6 + 2 ] = 0.0; g_u[pp + index_conF_eq_6 + 2 ] = 0.0;

    // eq5 but x y z c
    g_l[pp + index_conF_eq_7 + 0 ] = 0.0; g_u[pp + index_conF_eq_7 + 0 ] = 0.0;
    g_l[pp + index_conF_eq_7 + 1 ] = 0.0; g_u[pp + index_conF_eq_7 + 1 ] = 0.0;
    g_l[pp + index_conF_eq_7 + 2 ] = 0.0; g_u[pp + index_conF_eq_7 + 2 ] = 0.0;

  }

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
   // CoM Position
   x[index_CoM_pos+0] = desired[0].CoM_positions(0);
   x[index_CoM_pos+1] = desired[0].CoM_positions(1);
   x[index_CoM_pos+2] = desired[0].CoM_positions(2);

   for(int i=index_CoM_lin; i < index_CoM_ang; i++) //3
   {
     // CoM Lin
     x[i] = 0;
   }
   for(int i=index_CoM_ang; i < index_Forces; i++) //3
   {
     // CoM Ang
     x[i] = 0;
   }

   //init Forces, for each leg
   for(int ll=0; ll < p_agent->k_legs; ll++) //k_legs //TODO
   {
     // Init Foot Step Posiotion
     x[index_Forces +ll+0] = 0; //x
     x[index_Forces +ll+1] = 0; //y
     //split grav into each force foots
     x[index_Forces +ll+2] = (p_agent->mass * p_agent->g) / p_agent->k_legs; // z
   }

   // Steps, for each leg
   for(int ll=0; ll < p_agent->k_legs; ll++) //k_legs //TODO
   {
     // Init Foot Step Posiotion
     x[index_Step +ll+0] = desired[0].foot_position_v[ll](0); //x
     x[index_Step +ll+1] = desired[0].foot_position_v[ll](1); //y
     x[index_Step +ll+2] = desired[0].foot_position_v[ll](2); //z
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
  // double wl = 1, wk=0.8, wf=0.4, ws =0.4;
  //
  // for(int pp =0; pp < n_points; pp++)
  // {
  //   wl*sum2_l + wk*sum2_k + wf*sum2_Flegs + ws*sum2_Slegs
  //
  // }
  // obj_value = wl*
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
