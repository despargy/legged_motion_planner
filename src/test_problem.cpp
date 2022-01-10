#include <Eigen/Dense>
#include <iostream>
#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

using namespace Ipopt;

int main()
{
  Index n=-1; /* number of variables */
  Index m=-1; /* number of constraints */
  Number* x_L = NULL; /* lower bounds on x */
  Number* x_U = NULL; /* upper bounds on x */
  Number* g_L = NULL; /* lower bounds on g */
  Number* g_U = NULL; /* upper bounds on g */
  IpoptProblem nlp = NULL; /* IpoptProblem */
  enum ApplicationReturnStatus status; /* Solve return code */
  Number* x = NULL; /* starting point and solution vector */
  Number* mult_x_L = NULL; /* lower bound multipliers at the solution */
  Number* mult_x_U = NULL; /* upper bound multipliers at the solution */
  Number obj; /* objective value */
  Index i; /* generic counter */
  /* set the number of variables and allocate space for the bounds */
  n=4;
  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  /* set the values for the variable bounds */
  for (i=0; i<n; i++) {
  x_L[i] = 1.0;
  x_U[i] = 5.0;
  }
  /* set the number of constraints and allocate space for the bounds */
  m=2;
  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);
  /* set the values of the constraint bounds */
  g_L[0] = 25; g_U[0] = 2e19;
  g_L[1] = 40; g_U[1] = 40;
  /* create the IpoptProblem */
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, 8, 10, 0,
  &eval_f, &eval_g, &eval_grad_f,
  &eval_jac_g, &eval_h);
  /* We can free the memory now - the values for the bounds have been
  copied internally in CreateIpoptProblem */
  free(x_L);
  free(x_U);
  free(g_L);
  free(g_U);
  /* set some options */
  AddIpoptNumOption(nlp, "tol", 1e-9);
  AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
  /* allocate space for the initial point and set the values */
  x = (Number*)malloc(sizeof(Number)*n);
  x[0] = 1.0;
  x[1] = 5.0;
  x[2] = 5.0;
  x[3] = 1.0;
  /* allocate space to store the bound multipliers at the solution *
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);
  /* solve the problem */
  status = IpoptSolve(nlp, x, NULL, &obj, NULL, mult_x_L, mult_x_U, NULL);
  if (status == Solve_Succeeded) {
  printf("\n\nSolution of the primal variables, x\n");
  for (i=0; i<n; i++)
  printf("x[%d] = %e\n", i, x[i]);
  printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
  for (i=0; i<n; i++)
  printf("z_L[%d] = %e\n", i, mult_x_L[i]);
  for (i=0; i<n; i++)
  printf("z_U[%d] = %e\n", i, mult_x_U[i]);
  printf("\n\nObjective value\nf(x*) = %e\n", obj);
  }
  /* free allocated memory */
  FreeIpoptProblem(nlp);
  free(x);
  free(mult_x_L);
  free(mult_x_U);
  return 0;
}
