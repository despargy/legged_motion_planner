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

  n_per_dt = index_ends;

  index_conF_eq_1 = 0;
  index_conF_eq_2 = index_conF_eq_1 + p_agent->k_legs;
  index_conF_eq_3 = index_conF_eq_2 + p_agent->k_legs;
  index_conF_eq_4 = index_conF_eq_3 + p_agent->k_legs;
  index_conF_eq_5 = index_conF_eq_4 + p_agent->k_legs;
  index_conF_eq_6 = index_conF_eq_5 + 3;
  index_conF_eq_7 = index_conF_eq_6 + 3;

  g_per_dt = index_conF_eq_7 + 3; // g_per_dt has to be 17 in case of 2 legs

  index_partial_eq_1 = 0;
  index_partial_eq_2 = index_partial_eq_1 + p_agent->k_legs*2;
  index_partial_eq_3 = index_partial_eq_2 + p_agent->k_legs*2;
  index_partial_eq_4 = index_partial_eq_3 + p_agent->k_legs*2;
  index_partial_eq_5 = index_partial_eq_4 + p_agent->k_legs*2;
  index_partial_eq_6 = index_partial_eq_5 + 3*(1+p_agent->k_legs*1);
  index_partial_eq_7 = index_partial_eq_6 + 3*(2+p_agent->k_legs*4);


  // std::cout<<"n_per_dt = "<< n_per_dt<<std::endl;
  // std::cout<<"g_per_dt = "<< g_per_dt<<std::endl;
  // std::cout <<"index_partial_eq_7 + 5"<< index_partial_eq_7 + 6 << '\n';

  return true;

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

   n = n_per_dt*n_points; // n_points is without the init first
   std::cout<<"n = "<< n<<std::endl;

   // equality constraints and inequality constraints
   // Fx*legs +- -> 2*legs, Fy*legs +- -> 2*legs, l->3, k->3, CoM->3, -> 1 BUT for each dt (multi. with dt)
   m = ( 2*p_agent->k_legs + 2*p_agent->k_legs + 3 + 3 + 3)*n_points; //TODO check if equals ind7 * n_points
   std::cout<<"m = "<< m<<std::endl;

   // 2 nonzeros in the jacobian (one for x1, and one for x2),
   // sum index_partial_eq_7 + 2*3
   nnz_jac_g_per_dt = index_partial_eq_7 + 2*3;// p_agent->k_legs*8 + 3*(1 + p_agent->k_legs) + 3*(2 + 4*p_agent->k_legs) + 6 ;
   nnz_jac_g = nnz_jac_g_per_dt*n_points; //TODO;

   std::cout<<"nnz_jac_g_per_dt = "<< nnz_jac_g_per_dt<<"\n"<<std::endl;

   std::cout<<"nnz_jac_g = "<< nnz_jac_g<<"\n"<<std::endl;

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
  // assert(n == 21);
  // assert(m == 17);
  // std::cout<<"here get_bounds_info"<<std::endl;

  // set all as inf, lets say by default
  for (int i=0; i < n; i++)
  {
    x_l[i] = -1.0e19;
    x_u[i] = +1.0e19;
  }


  //CoM z bounds ground - height
  for(int pp=0; pp<n_points; pp++)
  {
    x_l[pp*n_per_dt + index_CoM_pos + 2] = p_agent->height; x_u[pp*n_per_dt + index_CoM_pos + 2] = p_agent->height ;
  }

  // constraints bounds

  // for each point
  for(int pp=0; pp<n_points; pp++)
  {
    // std::cout<<"pp*g_per_dt = "<<pp*g_per_dt<<"\n"<<std::endl;

    for(int ll=0; ll<p_agent->k_legs; ll++)
    {
      // eq1 0<= μFz - Fx
      g_l[pp*g_per_dt + index_conF_eq_1 + ll ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_1 + ll ] = +1.0e19;

      // eq2 μFz + Fx >= 0
      g_l[pp*g_per_dt + index_conF_eq_2 + ll ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_2 + ll ] = +1.0e19;

      // eq3 0<= μFz - Fy
      g_l[pp*g_per_dt + index_conF_eq_3 + ll ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_3 + ll ] = +1.0e19;

      // eq4 μFz + Fy >= 0
      g_l[pp*g_per_dt + index_conF_eq_4 + ll ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_4 + ll ] = +1.0e19;

    }

    // eq5 but x y z l
    g_l[pp*g_per_dt + index_conF_eq_5 + 0 ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_5 + 0 ] = 0.0;
    g_l[pp*g_per_dt + index_conF_eq_5 + 1 ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_5 + 1 ] = 0.0;
    g_l[pp*g_per_dt + index_conF_eq_5 + 2 ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_5 + 2 ] = 0.0;

    // eq6 but x y z k
    g_l[pp*g_per_dt + index_conF_eq_6 + 0 ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_6 + 0 ] = 0.0;
    g_l[pp*g_per_dt + index_conF_eq_6 + 1 ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_6 + 1 ] = 0.0;
    g_l[pp*g_per_dt + index_conF_eq_6 + 2 ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_6 + 2 ] = 0.0;

    // eq7 but x y z c
    g_l[pp*g_per_dt + index_conF_eq_7 + 0 ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_7 + 0 ] = 0.0;
    g_l[pp*g_per_dt + index_conF_eq_7 + 1 ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_7 + 1 ] = 0.0;
    g_l[pp*g_per_dt + index_conF_eq_7 + 2 ] = 0.0; g_u[pp*g_per_dt + index_conF_eq_7 + 2 ] = 0.0;

  }
  // for(int ll=0; ll<m; ll++)
  // {
  //     std::cout<< g_l[ll]<<" "<< g_u[ll]<<std::endl;
  // }
  // std::cout<<"end get_bounds_info"<<std::endl;

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
  // std::cout<<"here get_starting_point"<<std::endl;

   // Here, we assume we only have starting values for x, if you code
   // your own NLP, you can provide starting values for the others if
   // you wish.
   assert(init_x == true);
   assert(init_z == false);
   assert(init_lambda == false);

   // CoM Position
   x[index_CoM_pos+0] = desired[0].CoM_positions(0);
   x[index_CoM_pos+1] = desired[0].CoM_positions(1);
   x[index_CoM_pos+2] = desired[0].CoM_positions(2);

   // // CoM Lin
   x[index_CoM_lin + 0] = desired[0].CoM_linear_velocities(0);
   x[index_CoM_lin + 1] = desired[0].CoM_linear_velocities(1);
   x[index_CoM_lin + 2] = desired[0].CoM_linear_velocities(2);

   for(int i=index_CoM_ang; i < index_Forces; i++) //3
   {
     // CoM Ang
     x[i] = 0;
   }

   //init Forces, for each leg
   for(int ll=0; ll < p_agent->k_legs; ll++) //k_legs //TODO
   {
     // Init Foot Step Posiotion

    //multi 3* since each legs has x y z
     x[index_Forces + 3*ll + 0] = 0; //x
     x[index_Forces + 3*ll + 1] = 0; //y
     //split grav into each force foots
     x[index_Forces + 3*ll + 2] = (p_agent->mass * p_agent->g) / p_agent->k_legs; // z
   }

   // Steps, for each leg
   for(int ll=0; ll < p_agent->k_legs; ll++) //k_legs //TODO
   {
     // Init Foot Step Posiotion
     x[index_Step + 3*ll + 0] = desired[0].foot_position_v[ll](0); //x
     x[index_Step + 3*ll + 1] = desired[0].foot_position_v[ll](1); //y
     x[index_Step + 3*ll + 2] = desired[0].foot_position_v[ll](2); //z
   }
   // std::cout<<"here"<<std::endl;

   for(int pp=1; pp<n_points; pp++)
   {
     // std::cout<<"pp*g_per_dt = "<<pp*n_per_dt<<"\n"<<std::endl;

        // CoM Position
        x[pp*n_per_dt + index_CoM_pos+0] = x[(pp-1)*n_per_dt + index_CoM_pos+0];
        x[pp*n_per_dt + index_CoM_pos+1] = x[(pp-1)*n_per_dt + index_CoM_pos+1];
        x[pp*n_per_dt + index_CoM_pos+2] = x[(pp-1)*n_per_dt + index_CoM_pos+2];

        // // CoM Lin
        x[pp*n_per_dt + index_CoM_lin + 0] = x[(pp-1)*n_per_dt + index_CoM_lin + 0];
        x[pp*n_per_dt + index_CoM_lin + 1] = x[(pp-1)*n_per_dt + index_CoM_lin + 1];
        x[pp*n_per_dt + index_CoM_lin + 2] = x[(pp-1)*n_per_dt + index_CoM_lin + 2];

        for(int i=index_CoM_ang; i < index_Forces; i++) //3
        {
          // CoM Ang
          x[pp*n_per_dt + i] = x[(pp-1)*n_per_dt + i];
        }

        //init Forces, for each leg
        for(int ll=0; ll < p_agent->k_legs; ll++) //k_legs //TODO
        {
          // Init Foot Step Posiotion

         //multi 3* since each legs has x y z
          x[pp*n_per_dt + index_Forces + 3*ll + 0] = x[(pp-1)*n_per_dt + index_Forces + 3*ll + 0]; //x
          x[pp*n_per_dt + index_Forces + 3*ll + 1] = x[(pp-1)*n_per_dt + index_Forces + 3*ll + 1]; //y
          //split grav into each force foots
          x[pp*n_per_dt + index_Forces + 3*ll + 2] = x[(pp-1)*n_per_dt + index_Forces + 3*ll + 2]; // z
        }

        // Steps, for each leg
        for(int ll=0; ll < p_agent->k_legs; ll++) //k_legs //TODO
        {
          // Init Foot Step Posiotion
          x[pp*n_per_dt + index_Step + 3*ll + 0] = x[(pp-1)*n_per_dt + index_Step + 3*ll + 0]; //x
          x[pp*n_per_dt + index_Step + 3*ll + 1] = x[(pp-1)*n_per_dt + index_Step + 3*ll + 1]; //y
          x[pp*n_per_dt + index_Step + 3*ll + 2] = x[(pp-1)*n_per_dt + index_Step + 3*ll + 2]; //z
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
  // std::cout<<"here eval_f"<<std::endl;

  double sum2_l=0, sum2_k=0, sum2_Flegs=0, sum2_Slegs=0;

  for(int pp =0; pp < n_points; pp++) // Σ from 1 to N
  {
    // sum linears
    sum2_l += pow( x[pp*n_per_dt + index_CoM_lin + 0]-desired[pp+1].CoM_linear_velocities(0), 2); //x linear
    sum2_l += pow( x[pp*n_per_dt + index_CoM_lin + 1]-desired[pp+1].CoM_linear_velocities(1), 2); //y linear
    sum2_l += pow( x[pp*n_per_dt + index_CoM_lin + 2]-desired[pp+1].CoM_linear_velocities(2), 2); //z linear

    // sum angulars
    sum2_k += pow( x[pp*n_per_dt + index_CoM_ang + 0]-0, 2); //x angular
    sum2_k += pow( x[pp*n_per_dt + index_CoM_ang + 1]-0, 2); //y angular
    sum2_k += pow( x[pp*n_per_dt + index_CoM_ang + 2]-0, 2); //z angular

    for(int ll=0; ll < p_agent->k_legs; ll++) //k_legs //TODO
    {
      sum2_Flegs += pow(x[pp*n_per_dt + index_Forces + 3*ll + 0], 2) ; //x
      sum2_Flegs += pow(x[pp*n_per_dt + index_Forces + 3*ll + 1], 2) ; //y
      sum2_Flegs += pow(x[pp*n_per_dt + index_Forces + 3*ll + 2], 2) ; //z

      sum2_Slegs += pow(x[pp*n_per_dt + index_Step + 3*ll + 0]-desired[pp+1].foot_position_v[ll](0) ,2); //x
      sum2_Slegs += pow(x[pp*n_per_dt + index_Step + 3*ll + 1]-desired[pp+1].foot_position_v[ll](1) ,2); //y
      sum2_Slegs += pow(x[pp*n_per_dt + index_Step + 3*ll + 2]-desired[pp+1].foot_position_v[ll](2) ,2); //z

    }

  }

  obj_value = 1/2*(p_agent->wl*sum2_l + p_agent->wk*sum2_k + p_agent->wf*sum2_Flegs + p_agent->ws*sum2_Slegs);

  return true;
}

bool CentroidalNLP::eval_grad_f(
   Index         n,
   const Number* x,
   bool          new_x,
   Number*       grad_f
)
{
  // std::cout<<"here eval_grad_f"<<std::endl;

   // return the gradient of the objective function grad_{x} f(x)
   for(int pp =0; pp < n_points; pp++)
   {
     // std::cout<<"pp*g_per_dt = "<<pp*n_per_dt<<"\n"<<std::endl;

      grad_f[pp*n_per_dt + index_CoM_pos + 0] = 0;//  CoM x
      grad_f[pp*n_per_dt + index_CoM_pos + 1] = 0;//  CoM y
      grad_f[pp*n_per_dt + index_CoM_pos + 2] = 0;//  CoM z

      grad_f[pp*n_per_dt + index_CoM_lin + 0] = p_agent->wl*(x[pp*n_per_dt + index_CoM_lin + 0] - desired[pp+1].CoM_linear_velocities(0));//  linear x
      grad_f[pp*n_per_dt + index_CoM_lin + 1] = p_agent->wl*(x[pp*n_per_dt + index_CoM_lin + 1] - desired[pp+1].CoM_linear_velocities(1));//  linear y
      grad_f[pp*n_per_dt + index_CoM_lin + 2] = p_agent->wl*(x[pp*n_per_dt + index_CoM_lin + 2] - desired[pp+1].CoM_linear_velocities(2));//  linear z

      grad_f[pp*n_per_dt + index_CoM_ang + 0] = p_agent->wk*(x[pp*n_per_dt + index_CoM_ang + 0] - 0);//  angular x
      grad_f[pp*n_per_dt + index_CoM_ang + 1] = p_agent->wk*(x[pp*n_per_dt + index_CoM_ang + 1] - 0);//  angular y
      grad_f[pp*n_per_dt + index_CoM_ang + 2] = p_agent->wk*(x[pp*n_per_dt + index_CoM_ang + 2] - 0);//  angular z

      for(int ll=0; ll < p_agent->k_legs; ll++) //k_legs //TODO
      {
        grad_f[pp*n_per_dt + index_Forces + 3*ll + 0] = p_agent->wf*(x[pp*n_per_dt + index_Forces + 3*ll + 0]);//  i leg force x
        grad_f[pp*n_per_dt + index_Forces + 3*ll + 1] = p_agent->wf*(x[pp*n_per_dt + index_Forces + 3*ll + 1]);//  i leg force y
        grad_f[pp*n_per_dt + index_Forces + 3*ll + 2] = p_agent->wf*(x[pp*n_per_dt + index_Forces + 3*ll + 2]);//  i leg force z

        grad_f[pp*n_per_dt + index_Step + 3*ll + 0] = p_agent->ws*(x[pp*n_per_dt + index_Step + 3*ll + 0] - desired[pp+1].foot_position_v[ll](0));//  i leg step x
        grad_f[pp*n_per_dt + index_Step + 3*ll + 1] = p_agent->ws*(x[pp*n_per_dt + index_Step + 3*ll + 1] - desired[pp+1].foot_position_v[ll](1));//  i leg step y
        grad_f[pp*n_per_dt + index_Step + 3*ll + 2] = p_agent->ws*(x[pp*n_per_dt + index_Step + 3*ll + 2] - desired[pp+1].foot_position_v[ll](2));//  i leg step z
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
  // std::cout<<"here eval_g"<<std::endl;

   // return the value of the constraints: g(x)
   // return the value of the constraints: g(x)
   for(int pp =0; pp < n_points; pp++)
   {

     for(int ll=0; ll<p_agent->k_legs; ll++)
     {
       // eq1 0<= μFz - Fx
       g[pp*g_per_dt + index_conF_eq_1 + ll] = p_agent->mu*x[pp*n_per_dt+ index_Forces + 3*ll + 2] - x[pp*n_per_dt + index_Forces + 3*ll + 0];

       // eq2 μFz + Fx >= 0
       g[pp*g_per_dt + index_conF_eq_2 + ll] = p_agent->mu*x[pp*n_per_dt + index_Forces + 3*ll + 2] + x[pp*n_per_dt + index_Forces + 3*ll + 0];

       // eq3 0<= μFz - Fy
       g[pp*g_per_dt + index_conF_eq_3 + ll] = p_agent->mu*x[pp*n_per_dt + index_Forces + 3*ll + 2] - x[pp*n_per_dt + index_Forces + 3*ll + 1];

       // eq4 μFz + Fy >= 0
       g[pp*g_per_dt + index_conF_eq_4 + ll] = p_agent->mu*x[pp*n_per_dt + index_Forces + 3*ll + 2] + x[pp*n_per_dt + index_Forces + 3*ll + 1];

     }

       double sumf_x = 0, sumf_y = 0, sumf_z = 0;
       for(int ll=0; ll<p_agent->k_legs; ll++)
       {
         sumf_x += x[ pp*n_per_dt + index_Forces + 3*ll + 0];
         sumf_y += x[ pp*n_per_dt + index_Forces + 3*ll + 1];
         sumf_z += x[ pp*n_per_dt + index_Forces + 3*ll + 2];
       }

       g[pp*g_per_dt + index_conF_eq_5 + 0] = x[(pp+1)*n_per_dt + index_CoM_lin + 0] - (sumf_x + p_agent->mass*p_agent->g)*dt -  x[(pp)*n_per_dt + index_CoM_lin + 0]; //x
       g[pp*g_per_dt + index_conF_eq_5 + 1] = x[(pp+1)*n_per_dt + index_CoM_lin + 1] - (sumf_y + p_agent->mass*p_agent->g)*dt -  x[(pp)*n_per_dt + index_CoM_lin + 1]; //y
       g[pp*g_per_dt + index_conF_eq_5 + 2] = x[(pp+1)*n_per_dt + index_CoM_lin + 2] - (sumf_z + p_agent->mass*p_agent->g)*dt -  x[(pp)*n_per_dt + index_CoM_lin + 2]; //z

       double sum_cross_x = 0, sum_cross_y = 0, sum_cross_z = 0;
       for(int ll=0; ll<p_agent->k_legs; ll++)
       {
         // cross_P(0) = vect_A(0) * vect_B(0) - vect_A(2) * vect_B(1);
         sum_cross_x += x[(pp)*n_per_dt + index_CoM_pos + 0] - ( x[(pp)*n_per_dt + index_Step + 3*ll + 0]*x[(pp)*n_per_dt + index_Forces + 3*ll + 0]
                                                                                -  x[(pp)*n_per_dt + index_Step + 3*ll + 2]*x[(pp)*n_per_dt + index_Forces + 3*ll + 1]);

         // cross_P(1) = vect_A(2) * vect_B(0) - vect_A(0) * vect_B(2);
         sum_cross_y += x[(pp)*n_per_dt + index_CoM_pos + 1] - ( x[(pp)*n_per_dt + index_Step + 3*ll + 2]*x[(pp)*n_per_dt + index_Forces + 3*ll + 0]
                                                                                -  x[(pp)*n_per_dt + index_Step + 3*ll + 0]*x[(pp)*n_per_dt + index_Forces + 3*ll + 2]);

         // cross_P(2) = vect_A(0) * vect_B(1) - vect_A(1) * vect_B(0);
         sum_cross_z += x[(pp)*n_per_dt + index_CoM_pos + 2] - ( x[(pp)*n_per_dt + index_Step + 3*ll + 0]*x[(pp)*n_per_dt + index_Forces + 3*ll + 1]
                                                                                -  x[(pp)*n_per_dt + index_Step + 3*ll + 1]*x[(pp)*n_per_dt + index_Forces + 3*ll + 0]);

       }
       g[pp*g_per_dt + index_conF_eq_6 + 0] = x[(pp+1)*n_per_dt + index_CoM_ang + 0] - sum_cross_x*dt -  x[(pp)*n_per_dt + index_CoM_ang + 0]; //x
       g[pp*g_per_dt + index_conF_eq_6 + 1] = x[(pp+1)*n_per_dt + index_CoM_ang + 1] - sum_cross_y*dt -  x[(pp)*n_per_dt + index_CoM_ang + 1]; //y
       g[pp*g_per_dt + index_conF_eq_6 + 2] = x[(pp+1)*n_per_dt + index_CoM_ang + 2] - sum_cross_z*dt -  x[(pp)*n_per_dt + index_CoM_ang + 2]; //z

       g[pp*g_per_dt + index_conF_eq_7 + 0] = x[(pp+1)*n_per_dt + index_CoM_pos + 0] -  x[(pp)*n_per_dt + index_CoM_pos + 0] - 1/p_agent->mass*x[(pp)*n_per_dt + index_CoM_lin + 0]*dt; //x
       g[pp*g_per_dt + index_conF_eq_7 + 1] = x[(pp+1)*n_per_dt + index_CoM_pos + 1] -  x[(pp)*n_per_dt + index_CoM_pos + 1] - 1/p_agent->mass*x[(pp)*n_per_dt + index_CoM_lin + 1]*dt; //y
       g[pp*g_per_dt + index_conF_eq_7 + 2] = x[(pp+1)*n_per_dt + index_CoM_pos + 2] -  x[(pp)*n_per_dt + index_CoM_pos + 2] - 1/p_agent->mass*x[(pp)*n_per_dt + index_CoM_lin + 2]*dt; //z

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
  // std::cout<<"here eval_jac_g"<<std::endl;

      if (values == NULL) {
        // std::cout<<"here when NULL"<<std::endl;

        for(int pp =0; pp < n_points; pp++)
        {
          // return the structure of the Jacobian

          // this particular Jacobian is sparse
          // https://coin-or.github.io/Ipopt/IMPL.html

          // eq 1. (0<= μFz - Fx)/delta **i -> 2 nonzeros for each leg
          for(int ll=0; ll<p_agent->k_legs; ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_1 + 2*ll + 0] = pp*g_per_dt + index_conF_eq_1 + ll; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_1 + 2*ll + 0] = pp*n_per_dt + index_Forces + 3*ll + 0; // has /dFi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_1 + 2*ll + 1] = pp*g_per_dt + index_conF_eq_1 + ll; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_1 + 2*ll + 1] = pp*n_per_dt + index_Forces + 3*ll + 2; // has /dFi_z

            // std::cout << pp*nnz_jac_g_per_dt + index_partial_eq_1 + 2*ll + 0 << '\n';
            // std::cout << pp*nnz_jac_g_per_dt + index_partial_eq_1 + 2*ll + 1 << '\n';
          }

          // eq 2. (μFz + Fx >= 0)/delta **i -> 2 nonzeros for each leg
          for(int ll=0; ll<p_agent->k_legs; ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_2 + 2*ll + 0] = pp*g_per_dt + index_conF_eq_2 + ll; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_2 + 2*ll + 0] = pp*n_per_dt + index_Forces + 3*ll + 0; // has /dFi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_2 + 2*ll + 1] = pp*g_per_dt + index_conF_eq_2 + ll; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_2 + 2*ll + 1] = pp*n_per_dt + index_Forces + 3*ll + 2; // has /dFi_z
            // std::cout << pp*nnz_jac_g_per_dt + index_partial_eq_2 + 2*ll + 0 << '\n';
            // std::cout << pp*nnz_jac_g_per_dt + index_partial_eq_2 + 2*ll + 1 << '\n';
          }


          // eq 3. (0<= μFz - Fy)/delta **i -> 2 nonzeros for each leg
          for(int ll=0; ll<p_agent->k_legs; ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_3 + 2*ll + 0] = pp*g_per_dt + index_conF_eq_3 + ll; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_3 + 2*ll + 0] = pp*n_per_dt + index_Forces + 3*ll + 1; // has /dFi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_3 + 2*ll + 1] = pp*g_per_dt + index_conF_eq_3 + ll; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_3 + 2*ll + 1] = pp*n_per_dt + index_Forces + 3*ll + 2; // has /dFi_z

            // std::cout << pp*nnz_jac_g_per_dt + index_partial_eq_3 + 2*ll + 0 << '\n';
            // std::cout << pp*nnz_jac_g_per_dt + index_partial_eq_3 + 2*ll + 1 << '\n';
          }
          // eq 4. (μFz + Fy >= 0)/delta **i -> 2 nonzeros for each leg
          for(int ll=0; ll<p_agent->k_legs; ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_4 + 2*ll + 0] = pp*g_per_dt + index_conF_eq_4 + ll; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_4 + 2*ll + 0] = pp*n_per_dt + index_Forces + 3*ll + 1; // has /dFi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_4 + 2*ll + 1] = pp*g_per_dt + index_conF_eq_4 + ll; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_4 + 2*ll + 1] = pp*n_per_dt + index_Forces + 3*ll + 2; // has /dFi_z

            // std::cout << pp*nnz_jac_g_per_dt + index_partial_eq_4 + 2*ll + 0 << '\n';
            // std::cout << pp*nnz_jac_g_per_dt + index_partial_eq_4 + 2*ll + 1 << '\n';

          }

          //eq 5

          // l_x
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_5 + 0] = pp*g_per_dt + index_conF_eq_5 + 0; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_5 + 0] = pp*n_per_dt + index_CoM_lin + 0; // has /dl_x
          // std::cout <<  pp*nnz_jac_g_per_dt + index_partial_eq_5 + 0 << '\n';

          for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+inner_ll)] = pp*g_per_dt + index_conF_eq_5 + 0;
            jCol[pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+inner_ll)] = pp*n_per_dt + index_Forces + inner_ll*3 + 0; // has /dl_x

            // std::cout <<  pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+inner_ll) << '\n';

          }

          // l_y
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+p_agent->k_legs)] = pp*g_per_dt + index_conF_eq_5 + 1; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+p_agent->k_legs)] = pp*n_per_dt + index_CoM_lin + 1; // has /dl_x
          // std::cout << pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+p_agent->k_legs) << '\n';

          for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+p_agent->k_legs) + (1+inner_ll)] = pp*g_per_dt + index_conF_eq_5 + 1;
            jCol[pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+p_agent->k_legs) + (1+inner_ll)] = pp*n_per_dt + index_Forces + inner_ll*3 + 1; // has /dl_x

            // std::cout <<  pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+p_agent->k_legs) + (1+inner_ll) << '\n';

          }


          // l_z
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_5 + 2*(1+p_agent->k_legs)] = pp*g_per_dt + index_conF_eq_5 + 2;
          jCol[pp*nnz_jac_g_per_dt + index_partial_eq_5 + 2*(1+p_agent->k_legs)] = pp*n_per_dt + index_CoM_lin + 2; // has /dl_x

          // std::cout <<  pp*nnz_jac_g_per_dt + index_partial_eq_5 + 2*(1+p_agent->k_legs) << '\n';

          for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_5 + 2*(1+p_agent->k_legs) + (1+inner_ll)] = pp*g_per_dt + index_conF_eq_5 + 2;
            jCol[pp*nnz_jac_g_per_dt + index_partial_eq_5 + 2*(1+p_agent->k_legs) + (1+inner_ll)] = pp*n_per_dt + index_Forces + inner_ll*3 + 2; // has /dl_x

            // std::cout <<  pp*nnz_jac_g_per_dt + index_partial_eq_5 + 2*(1+p_agent->k_legs) + (1+inner_ll) << '\n';

          }

          //eq 6

          // k_x
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 0] = pp*g_per_dt + index_conF_eq_6 + 0; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 0] = pp*n_per_dt + index_CoM_pos + 0; // has /dc_x
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 1] = pp*g_per_dt + index_conF_eq_6 + 0; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 1] = pp*n_per_dt + index_CoM_ang + 0; // has /dk_x


          for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2 + 2*inner_ll + 0] = pp*g_per_dt + index_conF_eq_6 + 0; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2 + 2*inner_ll + 0] = pp*n_per_dt + index_Forces + 3*inner_ll + 0; // has /dfi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2 + 2*inner_ll + 1] = pp*g_per_dt + index_conF_eq_6 + 0; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2 + 2*inner_ll + 1] = pp*n_per_dt + index_Forces + 3*inner_ll + 1; // has /dfi_y

          }

          for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2+2*p_agent->k_legs) + 2*inner_ll + 0] = pp*g_per_dt + index_conF_eq_6 + 0; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2+2*p_agent->k_legs) + 2*inner_ll + 0] = pp*n_per_dt + index_Step + 3*inner_ll + 0; // has /dsi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2+2*p_agent->k_legs) + 2*inner_ll + 1] = pp*g_per_dt + index_conF_eq_6 + 0; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2+2*p_agent->k_legs) + 2*inner_ll + 1] = pp*n_per_dt + index_Step + 3*inner_ll + 2; // has /dsi_z
          }

          // k_y
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 0] = pp*g_per_dt + index_conF_eq_6 + 1; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 0] = pp*n_per_dt + index_CoM_pos + 1; // has /dc_y
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 1] = pp*g_per_dt + index_conF_eq_6 + 1; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 1] = pp*n_per_dt + index_CoM_ang + 1; // has /dk_y


          for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 0] = pp*g_per_dt + index_conF_eq_6 + 1; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 0] = pp*n_per_dt + index_Forces + 3*inner_ll + 0; // has /dfi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 1] = pp*g_per_dt + index_conF_eq_6 + 1; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 1] = pp*n_per_dt + index_Forces + 3*inner_ll + 2; // has /dfi_y

          }

          for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 0] = pp*g_per_dt + index_conF_eq_6 + 1; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 0] = pp*n_per_dt + index_Step + 3*inner_ll + 0; // has /dsi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 1] = pp*g_per_dt + index_conF_eq_6 + 1; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 1] = pp*n_per_dt + index_Step + 3*inner_ll + 2; // has /dsi_z
          }

          // k_z
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 0] = pp*g_per_dt + index_conF_eq_6 + 2; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 0] = pp*n_per_dt + index_CoM_pos + 2; // has /dc_z
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 1] = pp*g_per_dt + index_conF_eq_6 + 2; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 1] = pp*n_per_dt + index_CoM_ang + 2; // has /dk_z

          for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 0] = pp*g_per_dt + index_conF_eq_6 + 2; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 0] = pp*n_per_dt + index_Forces + 3*inner_ll + 0; // has /dfi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 1] = pp*g_per_dt + index_conF_eq_6 + 2; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 1] = pp*n_per_dt + index_Forces + 3*inner_ll + 1; // has /dfi_y

          }
          for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
          {
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 0] = pp*g_per_dt + index_conF_eq_6 + 2; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 0] = pp*n_per_dt + index_Step + 3*inner_ll + 0; // has /dsi_x
            iRow[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 1] = pp*g_per_dt + index_conF_eq_6 + 2; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 1] = pp*n_per_dt + index_Step + 3*inner_ll + 2; // has /dsi_z

          }
          //eq 7

          // c_x
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 0] = pp*g_per_dt + index_conF_eq_7 + 0; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 0] = pp*n_per_dt + index_CoM_pos + 0; // has /dc_x
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 1] = pp*g_per_dt + index_conF_eq_7 + 0; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 1] = pp*n_per_dt + index_CoM_lin + 0; // has /dl_x

          // c_y
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 2] = pp*g_per_dt + index_conF_eq_7 + 1; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 2] = pp*n_per_dt + index_CoM_pos + 1; // has /dc_y
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 3] = pp*g_per_dt + index_conF_eq_7 + 1; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 3] = pp*n_per_dt + index_CoM_lin + 1; // has /dl_y


          // c_z
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 4] = pp*g_per_dt + index_conF_eq_7 + 2; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 4] = pp*n_per_dt + index_CoM_pos + 2; // has /dc_z
          iRow[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 5] = pp*g_per_dt + index_conF_eq_7 + 2; jCol[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 5] = pp*n_per_dt + index_CoM_lin + 2; // has /dl_z


        }
        // std::cout<<"end when NULL"<<std::endl;

        // for(int l=0; l<n_points*nnz_jac_g_per_dt; l++)
        // {
        //   std::cout << iRow[l] <<" "<<jCol[l]<< '\n';
        //   std::cout << "l= "<<l << '\n';
        //
        // }
      }
      else
      {
        // std::cout<<"here else not NULL"<<std::endl;

        for(int pp =0; pp < n_points; pp++)
        {
          // return the values of the Jacobian of the constraints
  /////////////////////////////////////////
            //eq 1
            for(int ll=0; ll<p_agent->k_legs; ll++)
            {


              for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
              {

                if( inner_ll == ll)
                {
                  values[pp*nnz_jac_g_per_dt + index_partial_eq_1 + ll*2 + 0] = -1; //dfi_x
                  values[pp*nnz_jac_g_per_dt + index_partial_eq_1 + ll*2 + 1] = p_agent->mu; //dfi_z
                }

              }


            }


  //////////////////////////////////////////////
            //eq 2
            for(int ll=0; ll<p_agent->k_legs; ll++)
            {

              for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
              {

                if( inner_ll == ll)
                {
                  values[pp*nnz_jac_g_per_dt + index_partial_eq_2 + ll*2 + 0] = 1; //dfi_x
                  values[pp*nnz_jac_g_per_dt + index_partial_eq_2 + ll*2 + 1] = p_agent->mu; //dfi_z
                }

              }


            }


  ////////////////////////////////////
            //eq 3
            for(int ll=0; ll<p_agent->k_legs; ll++)
            {

              for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
              {

                if( inner_ll == ll)
                {
                  values[pp*nnz_jac_g_per_dt + index_partial_eq_3 + ll*2 + 0] = -1; //dfi_y
                  values[pp*nnz_jac_g_per_dt + index_partial_eq_3 + ll*2 + 1] = p_agent->mu; //dfi_z
                }

              }

            }

  ///////////////////////
            //eq 4
            for(int ll=0; ll<p_agent->k_legs; ll++)
            {

              for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
              {

                if( inner_ll == ll)
                {
                  values[pp*nnz_jac_g_per_dt + index_partial_eq_4 + ll*2 + 0] = 1; //dfi_y //ll*n_per_dt +
                  values[pp*nnz_jac_g_per_dt + index_partial_eq_4 + ll*2 + 1] = p_agent->mu; //dfi_z //ll*n_per_dt +
                }

              }

            }

  ////////////////////
            //eq 5

            // l_x
            values[pp*nnz_jac_g_per_dt + index_partial_eq_5 + 0] = -1; //dl_x

            for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
            {
              values[pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+inner_ll) ] = -dt; //dfi_x
            }

            // l_y
            values[pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+p_agent->k_legs) ] = -1; //dl_y


            for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
            {
              values[pp*nnz_jac_g_per_dt + index_partial_eq_5 + (1+p_agent->k_legs) + (1+inner_ll)] = -dt; //dfi_y
            }

            // l_z
            values[pp*nnz_jac_g_per_dt + index_partial_eq_5 + 2*(1+p_agent->k_legs)] = -1; //dl_z


            for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
            {
              values[pp*nnz_jac_g_per_dt + index_partial_eq_5 + 2*(1+p_agent->k_legs) + (1+inner_ll)] = -dt; //dfi_z
            }

  ////////////////////
            //eq 6

            // k_x
            values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 0] = -2*dt; //dc_x
            values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 1] = -1; //dk_x


            for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
            {
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2 + 2*inner_ll + 0] =  x[pp*n_per_dt + index_Step + 3*inner_ll + 0]*dt; //dfi_x
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2 + 2*inner_ll + 1] = -x[pp*n_per_dt + index_Step + 3*inner_ll + 2]*dt; //dfi_y
            }

            for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
            {
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2+2*p_agent->k_legs) + 2*inner_ll + 0] =  x[pp*n_per_dt + index_Forces + 3*inner_ll + 0]*dt; //dsi_x
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2+2*p_agent->k_legs) + 2*inner_ll + 1] = -x[pp*n_per_dt + index_Forces + 3*inner_ll + 1]*dt; //dsi_z
            }

            // k_y
            values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 0] = -2*dt; //dc_y
            values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 1] = -1; //dk_y

            for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
            {
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 0] =  x[pp*n_per_dt + index_Step + 3*inner_ll + 2]*dt; //dfi_x
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 1] = -x[pp*n_per_dt + index_Step + 3*inner_ll + 0]*dt; //dfi_z
            }

            for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
            {
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 0] = -x[pp*n_per_dt + index_Forces + 3*inner_ll + 2]*dt; //dsi_x
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + (2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 1] =  x[pp*n_per_dt + index_Forces + 3*inner_ll + 0]*dt; //dsi_z
            }

            // k_z
            values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 0] = -2*dt; //dc_z
            values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 1] = -1; //dk_z

            for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
            {
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 0] = -x[pp*n_per_dt + index_Step + 3*inner_ll + 1]*dt; //dfi_x
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + 2 + 2*inner_ll + 1] =  x[pp*n_per_dt + index_Step + 3*inner_ll + 0]*dt; //dfi_y
            }
            for(int inner_ll=0; inner_ll<p_agent->k_legs; inner_ll++)
            {
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 0] =  x[pp*n_per_dt + index_Forces + 3*inner_ll + 1]*dt; //dsi_x
              values[pp*nnz_jac_g_per_dt + index_partial_eq_6 + 2*(2 + p_agent->k_legs*4) + (2+2*p_agent->k_legs) + 2*inner_ll + 1] = -x[pp*n_per_dt + index_Forces + 3*inner_ll + 0]*dt; //dsi_y
            }
  ////////////////////
            //eq 7

            // c_x
            values[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 0] = -1; //dc_x
            values[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 1] = -dt/p_agent->mass; //dl_x


            // c_y
            values[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 2] = -1; //dc_y
            values[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 3] = -dt/p_agent->mass; //dl_y

            // c_z
            values[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 4] = -1; //dc_z
            values[pp*nnz_jac_g_per_dt + index_partial_eq_7 + 5] = -dt/p_agent->mass; //dl_z

        }

        // for(int l=0; l<n_points*nnz_jac_g_per_dt; l++)
        // {
        //   std::cout << values[l]<< '\n';
        //   std::cout << "l= "<<l << '\n';
        //
        // }
    }
    std::cout<<"end eval_jac_g"<<std::endl;

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

   std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
   for( Index i = 0; i < n; i++ )
   {
     std::cout << "x[" << i << "] = " << x[i] << std::endl;
   }

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

// Function to find
// cross product of two vector array.
Eigen::Vector3d CentroidalNLP::crossProduct(Eigen::Vector3d vect_A, Eigen::Vector3d vect_B)

{
    Eigen::Vector3d cross_P;
    cross_P(0) = vect_A(0) * vect_B(0) - vect_A(2) * vect_B(1);
    cross_P(1) = vect_A(2) * vect_B(0) - vect_A(0) * vect_B(2);
    cross_P(2) = vect_A(0) * vect_B(1) - vect_A(1) * vect_B(0);

    return cross_P;
}
