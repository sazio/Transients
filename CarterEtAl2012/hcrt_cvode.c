/*
 This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Header files with a description of contents used here */

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, functions, and macros */

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.*/

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
//#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NEQ   7                /* number of equations  */
#define ATOL  RCONST(1.0e-14)   /* vector absolute tolerance components */
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.4)      /* first output time      */
#define TMULT RCONST(0.01)     /* output time factor     */


/* Functions Called by the Solver */

static int f(double t, N_Vector y, N_Vector ydot, void *user_data);

/* Private functions to output results */

static void PrintOutput(double t, double y1, double y2);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

void set_initial_conditions(N_Vector );

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */
#define PARAMETERS 4
#define TIMEMAX   1
#define TOLERANCE 2
#define IDC       3

double Idc;

int main(int argc, char **argv)
{
  double t, tout;
  N_Vector y=NULL;
  void *cvode_mem=NULL;
  int flag;
  double tmax;

  if (argc!=PARAMETERS) {
    puts("Max time in seconds");
    puts("Relative tolerance");
    puts("Current");
    return -1;
  }

  Idc=atof(argv[IDC]);
  tmax=atof(argv[TIMEMAX])*1000;

  /* Create serial vector of length NEQ for I.C. */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  /* Initialize y */
  set_initial_conditions(y);

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  
  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  flag = CVodeSStolerances(cvode_mem, atof(argv[TOLERANCE]), ATOL);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);
  /* Set the Jacobian routine to internal estimation */
  flag=CVDlsSetDenseJacFn(cvode_mem, NULL);
  
  // flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  // if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */

  tout = T1;
  while(1) {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    PrintOutput(t, Ith(y,1), Ith(y,7));

    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS) {
      tout += TMULT;
      //printf("%lg %lg\n",tout,TMULT);
    }

    if (tout>tmax) break;
  }

  /* Free y vector */
  N_VDestroy_Serial(y);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  return(0);
}


/*
 * f routine. Compute function f(t,y). 
 */
#define CA  10.0 /* pico F */
#define CS  10.0 /* pico F */
#define gl  1.6 /* nano S */
#define gAS 65.0  /* nano S */ 
#define gNa 260.0   /* nano S */
#define gKd 80.0   /* nano S */
#define gCa 0.88   /* nano S */
#define MU  1.5    /* Calcium dynamics dissipation*/

#define ABS( X ) (((X)>0.0)?(X):-(X))
#define Cagam(X,Y,Z) (1.0/(1.0+exp(((X)-(Y))/(Z)))) 

static int f(double t, N_Vector y, N_Vector ydot, void *user_data)
{
  double Va=Ith(y,1);
  double Vs=Ith(y,2);
  double Vt=-50.14905;
  double v;
  double a;
  double INa, IKd, ICa;
  
  /* INa current [Traub & Miles, 1991] */
  INa=gNa*Ith(y,3)*Ith(y,3)*Ith(y,4)*(Va-50);

  /* INa activation Ith(y,3)*/
  v=18.0+Vt-Va;
  if (ABS(v)<0.0001) {
    a=1.28-0.16*v;
  } else 
    a=0.32*v/(exp(v/4.0)-1);    
  
  Ith(ydot,3)=a*(1-Ith(y,3));

  v=-40.0-Vt+Va;
  if (ABS(v)<0.0001) {
    a=1.40-0.14*v;
  } else 
    a=0.28*v/(exp(v/5.0)-1.0);
  Ith(ydot,3)-=a*(Ith(y,3));

  /* INa inactivation Ith(y,4) */
  v=17+Vt-Va;
  Ith(ydot,4)=0.128*exp(v/18.0)*(1-Ith(y,4));
  v=40.0+Vt-Va;
  Ith(ydot,4)-=4.0*Ith(y,4)/(1.0+exp(v/5.0));

  /* IKd delayed rectifier current */
  IKd=gKd*Ith(y,5)*(Va+60);
  
  /* IKd activation Ith(y,5) */
  v=35.0+Vt-Va;
  if (ABS(v)<0.0001) {
    a=0.08-0.008*v;
  } else 
    a=0.016*v/(exp(v/5.0)-1.0);
  v=20.0+Vt-Va;
  Ith(ydot,5)=a*(1-Ith(y,5))-0.25*exp(v/40)*Ith(y,5);
  
  /* ICa current  */
  if (ABS(Vs)<0.001) {
    a=0.5*(- 24.42002442+Vs);
  } else {
    a=Vs/(1-exp(2*Vs/ 24.42002442));
  }
  ICa=gCa*Ith(y,6)*Ith(y,6)*Ith(y,6)*a;
  /* Calcium activation Ith(y,6) */
  Ith(ydot,6)=0.1*(Cagam(-Vs,39.1,2)-Ith(y,6));
  //printf("%lg ydot %lf\n",ICa,Ith(ydot,6));exit(0);
  /* Axon membrane potential */

  /* [Ca] Ith(y,7)*/
  Ith(ydot,7)=0.001*(-0.35*ICa-MU*MU*Ith(y,7)+0.04*MU*MU);

  Ith(ydot,1) = (1.0/CA)*(-gl*(Va+45.0 /*resting potential*/)-gAS*(Va-Vs)
			  -INa-IKd);
  
  //printf("ydo1=%lf \n",Ith(ydot,1));
  //exit(0);
  /* Soma membrane pontential */
  Ith(ydot,2) = (1.0/CS)*(-gl*(Vs+45.0 /*resting potential*/)-gAS*(Vs-Va)
			  +ICa
			 +Idc);
  //printf("ydo1=%lg (%lf,%lf,%lf,%lf)\n",Ith(ydot,1),gl*(Va+45.0 /*resting potential*/),gAS*(Va-Vs),INa,IKd);
  //exit(0);
  return(0);      
}


/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(double t, double y1, double y2)
{

  printf("%0.4e %14.6e  %14.6e\n", t, y1, y2*100);

  return;
}

void set_initial_conditions(N_Vector y) {
  int i;
  for(i=3;i<=NEQ;++i)   Ith(y,i)=0.1;
  
  Ith(y,1) = -45.0;
  Ith(y,2) = -46.0;
  Ith(y,7) = 0.04;

}

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
