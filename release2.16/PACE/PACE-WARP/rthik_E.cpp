#include "mex.h"        // to compile in Matlab
#include <iostream> 	// to check outputs here and there
#include <cstdlib>    	// to use srand (if you uncomment the relevant line)
#include <string>   	// to read in the kernel name
#include <vector>       // to construct the index list
#include <map>	    	// to map kernels to integers for the switch
#include <Eigen/Core>   // to use basic Eigen structures
using namespace Eigen;
using namespace std; 

//Commented lines contain various checks on integrity

//Date of last editing: Saturday 8th of September 2012

float LinearInterpolation (VectorXf X , VectorXf Y, float X_PointOfInterest){
//Produce Y_point_of_interest given X,Y and target X_point_of_interest
//X : vector containing the X variables of the interpolant 
//Y : vector containing the Y variables of the interpolant 
//PointOfInterest : Point of X to estimate the new point of Y
    
  float   xk, xkp1, yk, ykp1;  //Points adjecent to the point of interpolation
  if ( X.size() != Y.size() ){cout << "Problem with vector sizes" << endl; return(-1);}
//cout <<  " X(0): " <<  X(0) <<" X(Y.size()-1): " <<X(Y.size()-1)   <<   " Point of interest: " << X_PointOfInterest<< endl;
  if ( X_PointOfInterest < X(0) || X_PointOfInterest > X(Y.size()-1) ){cout << "You interpolate out of the curve boundaries" << endl; return(-1);}
//Find the points right before and right after the point of interest
  for (int i=1; i<X.size() ; i++){
    if (X(i)>= X_PointOfInterest){
      xkp1 = X(i);
      xk = X(i-1);
      ykp1 = Y(i);
      yk = Y(i-1);
      break;}
  }
//point-slope form for a line formula
  float t = (X_PointOfInterest -xk)/(xkp1 -xk);
  float yPOI = (1-t) * yk + t * ykp1;  // estimate point of interest
// cout << "(" << xk << ",  " << X_PointOfInterest << " , " << xkp1 << ") & (" << yk << ", " << yPOI  << ", " << ykp1 << ")"<< endl;
  return (yPOI);   
} 

VectorXf fnvalspapi(VectorXf X, VectorXf Y, VectorXf X_target){
  //evaluate Y_target for X_target given X and Y
    
	int N =  X_target.size();
	VectorXf rr(N); 
	for (int i=0; i<N ; i++){  
      rr(i) =  LinearInterpolation(X, Y, X_target(i)) ;  }
	return(rr);
}

int MonotonicityCheck(VectorXf X){
   // evaluate vector monotonicity of vector X
   int N =  X.size(); 
   for (int i=0; i<(N-1) ; i++){  
     if ( X(i) >X (1+i) ) return (-1);   }
   return (1);
}

float rttemp1( VectorXf t_reg,  VectorXf curvei, VectorXf curvek,int nknots, float lambda, VectorXf initial){
  //Calculate the cost of the given warping
  // t_reg  : time grid of y_reg
  // curvei : query curve
  // curvek : reference curve
  // nknots : number of knots
  // lambda : time distortion penalty parameter
  // initial: position of knots on the simplex 
    
  int N = t_reg.size();
  VectorXf struct_ = VectorXf::LinSpaced( nknots+2, t_reg(0)  , t_reg.maxCoeff() ); 
  VectorXf hik(N);
  VectorXf Q(2+ initial.size()) ;           //Solution with the placement of the knots on the simplex
  Q(0) = t_reg(0) ; Q(1+ initial.size()) = t_reg.maxCoeff()  ; Q.segment(1, initial.size()) = initial;
  hik =  fnvalspapi( struct_, Q , t_reg);   // compute the new internal time-scale  
  
  //cout << "hik: " << hik.transpose() << endl;
  //cout << "Monotonicity Checked on Hik: " << MonotonicityCheck(hik) << endl;
  //if(  MonotonicityCheck(hik) == -1) { cout <<" Q.transpose()  is :"<< Q.transpose()  << endl;}
  
  return ( (fnvalspapi(t_reg,curvei,hik)-fnvalspapi(t_reg,curvek ,t_reg)).array().pow(2).sum() + lambda * (hik - t_reg).array().pow(2).sum() );
}

VectorXf NewSolution( VectorXf x0 , float Step, int point, VectorXf t_reg){ 
    //generate new solution on the simplex defined in [t_reg(0)-x0-t_reg(N-1)] using a displacement of size Step
    // x0   : initial solution
    // Step : displacement size
    // point: knot to perturb
    // t_reg: time_grid
    
    VectorXf InitConf (2+ x0.size());
    InitConf(0) = t_reg(0); InitConf( x0.size()+1) = t_reg.maxCoeff()  ; InitConf.segment(1, x0.size()) = x0; 
    float LowBou =  InitConf(point-1);
    float UppBou =  InitConf(point+1);
    float New_State = (UppBou - LowBou) * Step + LowBou; 
    InitConf(point) = New_State; 
    //if (  MonotonicityCheck( InitConf.segment(1, x0.size()) ) == -1) { 
    //    cout << "We generated a unacceptable solutiion" << endl << "Initial seed was : " << x0.transpose() 
    //         << endl << " and we produced :"<<InitConf.segment(1, x0.size()).transpose() << endl;}
    return (InitConf.segment(1, x0.size()) ) ;
  } 
    
//====
typedef struct { float Val; VectorXf Mapping;} rthink_Output; //define structure that holds the results.
//====

rthink_Output rthik_SA(VectorXf t_reg,  VectorXf curvei, VectorXf curvek,int nknots, float lambda){
  // Random Search solver for the optimization problem of pairwise warping
  // t_reg : time_grid
  // curvei: query curve
  // curvek: reference curve 
  // nknots : number of knots
  // lambda : time distortion penalty parameter
    
  int k=0; float OldSol, bk; 
  VectorXf xk(nknots);   
  bk = (t_reg.maxCoeff() - t_reg(0))/(1+ float(nknots ));   //Distance between adjacent knots and edges-knots
  xk  = VectorXf::LinSpaced(nknots , t_reg(0) + bk  ,  - bk + t_reg.maxCoeff()  ); //Initial candidate solution with equispaced knots 
  VectorXf xn(nknots);   VectorXf help(2+nknots); 
 
  float NewSol; 
  OldSol =  rttemp1(t_reg, curvei, curvek, nknots , lambda, xk);  //Cost of initial solution
  int z= 99*nknots;                                               //Number of random search to do (proportional to the # of knots)
  //srand(1);                                                     //Fix the seed to have reproducable behaviour
  VectorXf Steps(z); Steps = (ArrayXf::Random(z)+1.)/2.;          //Generate possible random pertubations magnitude
  VectorXi Posit(z); for (int u=0; u <z; u++ ) Posit(u) =1+ rand()%(nknots+0); //Generate list of positions to purturb
 
  k=0;
  while((OldSol > .0001) && (k<z)) {
    xn = NewSolution(xk,Steps(k), Posit(k), t_reg );              //Get a new solution
    NewSol = rttemp1(t_reg, curvei, curvek, nknots, lambda, xn);  //Cost of new solution
    if (  (NewSol < OldSol)  ) {                                  //If it's better than the old one, use it.
      OldSol= NewSol; xk= xn; 
    }
  k++;
}

  VectorXf x3 =  VectorXf::LinSpaced(2+nknots, t_reg(0), t_reg( t_reg.size()-1)); 
  help(0) =  t_reg(0) ; help(nknots+1) = t_reg( t_reg.size()-1) ; help.segment(1,nknots) = xk;  

  rthink_Output G;
  G.Val = OldSol;
  G.Mapping = fnvalspapi(  x3 , help  , t_reg);
  return G ;
}

//Use with syntax :  [hik(indd,:),ftemp]=rthik_E( curvei.coefs,curvek.coefs,t_reg,nknots,lambda); 
#define IS_REAL_1D_FULL_DOUBLE(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))

void mexFunction (int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{    
    #define p_curvei	prhs[0]  
	#define p_curvek	prhs[1]  
	#define p_treg		prhs[2]  
	#define p_knots		prhs[3]  
	#define p_lambda	prhs[4]  
               
    #define p_hik         plhs[0]
    #define p_f           plhs[1]        
     
    if(nrhs != 5)       // Check the number of input arguments  
    mexErrMsgTxt("Wrong number of input arguments.");            
    else if(!IS_REAL_1D_FULL_DOUBLE(p_curvei))   //  Check lw
     mexErrMsgTxt("curvei needs be a real 2D full double array."); 
    else if(!IS_REAL_1D_FULL_DOUBLE(p_curvek))   //  Check lw
     mexErrMsgTxt("curvek needs be a real 2D full double array."); 
    else if(!IS_REAL_1D_FULL_DOUBLE(p_treg))   //  Check lw
     mexErrMsgTxt("t_reg needs be a real 2D full double array.");   
   
   const int LengthOfVector = mxGetN(p_treg); 
   
   //Read in the variables
   Map<VectorXd> CurvI(mxGetPr(p_curvei) , LengthOfVector); 
   Map<VectorXd> CurvK(mxGetPr(p_curvek) , LengthOfVector); 
   Map<VectorXd> T_Reg(mxGetPr(p_treg)   , LengthOfVector); 
        
   int knots = *mxGetPr(p_knots);
   double lamdba = *mxGetPr(p_lambda);
   srand(12);    
   //Past as floats cause really nobody cares about the small digits there.
   rthink_Output Res = rthik_SA(T_Reg.cast<float>(),CurvI.cast<float>(),CurvK.cast<float>(), knots, float(lamdba));
   
  //Declare the outputs 
  p_hik = mxCreateDoubleMatrix(1, LengthOfVector, mxREAL);
  p_f   = mxCreateDoubleMatrix(1, 1, mxREAL); 
  
  double *hik = mxGetPr(p_hik );
  double *f = mxGetPr(p_f );
  for (int u=0; u< LengthOfVector; u++) { hik[u] = Res.Mapping(u);}
  f[0] = Res.Val;
 
}
