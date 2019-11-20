/*bessel_mod.c: A fast bessel function calculator, useful for evaluating
 *              bessel functions J(n, z) and J'(n, z) for large n. Adapted
 *              with permission from Leung, Gammie, and Noble (2011).
 *
 * References:
 * [1] Pandya, Zhang, Chandra, Gammie (Accepted to ApJ 2016)
 * [2] Leung, Gammie, and Noble (2011)
 * [3] Xiao (2006)
 */

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>

/* other C header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#define C_pi       (3.1415926535897932384       ) /* pi */
#define BESSEL_EPSILON_ORDER (16) 


// Parameter for lines that deliminate between various approximations to J_n(x)
//  ( values found empirically for n = 100..1e7 )
#define SLOPE1 (-6.627757624078600696e-01)
#define SLOPE2 (-6.656260931106707801e-01)
#define SLOPE3 (-6.543033585865805080e-01)
#define INTERCEPT1 (1.063380408975875602e+00)
#define INTERCEPT2 (2.563324856985127465e-01)
#define INTERCEPT3 (2.720161927383055733e-01)
#define N_JN 	(30.)

double exp_factor(  double f_factor, double f_exp )  ;

/******************************************************************************************/
/******************************************************************************************
   my_Bessel_dJ():
   ----------------

  //#if FLAG_JNprime_EQ == JNprime_EQ1
       -- returns the derivative of J_n(x) based on recurrence relation:

            J_n'(x) = -J_{n+1}(x) + J_n(x)*(n/x)

       -- I use this relation so we do not have to check for n=0 case. 
  //#elif FLAG_JNprime_EQ == JNprime_EQ2
       -- returns the derivative of J_n(x) based on recurrence relation:

            2 J_n'(x) = J_{n-1}(x) + J_{n+1}(x)
  //#endif

******************************************************************************************/
double my_Bessel_dJ( double n, double x) 
{
  double jnp1;
  double my_Bessel_J( double n, double x ) ;
  double bessel_func;
  
  bessel_func  = my_Bessel_J(   n,   x );
  jnp1 = my_Bessel_J( (n+1), x );
#if FLAG_JNprime_EQ == JNprime_EQ1
  // **** problem: how about if n is between 0 and 1?
  if(x == 0.) {
    if(n >= 2.) return(0.); /* J_n(0) = 0 for n >= 1, then the recursive relation gives
                             * a zero derivative for n >= 2 */
    if(n == 0.) return(-jnp1); /* d(J_0(z))/dz = -J_1(z) */
    return( (n*bessel_func)/(x+DBL_MIN) - jnp1 );
  }
  return( n*(bessel_func)/x - jnp1 );
#elif FLAG_JNprime_EQ == JNprime_EQ2
  double jnm1 = my_Bessel_J( (n-1), x );

  // **** problem: how about if n is between 0 and 1?
  if(x == 0.) {
    if(n >= 2.) return(0.); /* J_n(0) = 0 for n >= 1, then the recursive relation gives
                             * a zero derivative for n >= 2 */
    if(n == 0.) return(-jnp1); /* d(J_0(z))/dz = -J_1(z) */
  }

  /*  2 d(J_n(z))/dz = J_{n-1}(z) - J_{n+1)(z) */
  return( 0.5*(jnm1-jnp1) );
#endif
}

/******************************************************************************************/
/******************************************************************************************
   my_Bessel_J():
   ----------------
       -- requires x > 0, and  n > 0; 
       -- only recommend for  n > 100  ,  use glibc function jn() otherwise ;
       -- returns approximate value of the Bessel function "J_n(x)" ;
       -- approximations based on those given in Chishtie et al. 2005
******************************************************************************************/
double my_Bessel_J( double n, double x ) 
{ 
  double fn, y, logn;
  double BesselJ_Debye_Eps_Exp(  double n, double x ) ;
  double BesselJ_Meissel_First(  double n, double x ) ;
  double BesselJ_Meissel_Second( double n, double x ) ;
  double BesselJ_Asympt1(        double n, double x);
  double BesselJ_bigx(           double n , double x );

//#if FLAG_N_JN == JN_C_LIB
  /* At least for GNU C Library, jn(n,x) requires n to be integer. */
 if( n < N_JN ) { 
    return(gsl_sf_bessel_Jn((int)n,x));
  }
//#elif FLAG_N_JN == JN_C_LIB_interpolate
//  double dn;
//  int ni;

  /* At least for GNU C Library, jn(n,x) requires n to be integer. Therefore linear 
   *  interpolation is used to calculate J_n. */
  /* If we just return jn(n,x), the bessel function at small n is a step function with 
   *  jumps at every integer. The integrand of n integration becomes spiky for n < N_JN */
//  if( n < N_JN ) { 
//    ni = (int)n;
//    dn = n - ni;

//    return( (1.-dn)*jn(ni,x) + dn*jn(ni+1,x) );
//  }
//#elif FLAG_N_JN == JN_GSL
//  if( n < N_JN ) { 
//    return( gsl_sf_bessel_Jnu(n,x) );
//  }
//#elif FLAG_N_JN == JN_Chishtie
//#else
//  #error "FLAG_N_JN is not set correctly"
//#endif

  fn = n;
  logn = log10(fn);
  
  if( x < fn )  {
    y = log10( (fn - x)/fn );

    if( x!=fn ) { 
      // May be useful to speed things up :    
//      if( y > (SLOPE1*logn+INTERCEPT1)  ) { 
//	return( 0. );
//      }
      if( y > (SLOPE2*logn+INTERCEPT2)  ) { 
	return( BesselJ_Meissel_First( n , x ) );
      }
    }
    return( BesselJ_Debye_Eps_Exp( n , x ) );
  }
  else { 
    y = log10( (x - fn)/x );

    if( (x!=fn) && (y > (SLOPE3*logn+INTERCEPT3))  ) { 
      return( BesselJ_Meissel_Second( n , x ) );
    } 
    return( BesselJ_Debye_Eps_Exp( n , x ) );
  }


}

/******************************************************************************************/
/******************************************************************************************
  BesselJ_Debye_Eps_Exp():
  -----------------------
        -- the x -> n  limits from either side;
        -- valid for a very narrow region around  x=n;
        -- expansion w.r.t. variable "eps",  eps = x - n 
      
        -- because of the order of the eps used, this routine is limited to x < 1e55;

 ******************************************************************************************/

double BesselJ_Debye_Eps_Exp( double n , double x ) 
{ 

  static int  first_time = 1 ;
  double z, ez;
  double t3, t4, t10, t38, t44, t70, t93, t107, t114, t146,t149;
  static double At[BESSEL_EPSILON_ORDER];
  void  set_At( double At[BESSEL_EPSILON_ORDER] );

  if( first_time ) { 
    set_At(At);
    first_time = 0;
  }

  if( x > 1.e55 ) { 
    return( 0. ); 
  }

  ez = x - n;
  z = pow( x, (1./3.) );

  t3 = z * z;
  //t4 =  t3 * t3;   // z^4 = x^(4/3)
  t4 =  x * z ;      // z^4 = x^(4/3)
  t10 = t4 *  t4;    // z^8 = x^(8/3)
  t38 = 810485676000000 * At[4] * t3;
  //t44 = t3 * z;
  t44 = x;
  t70 = At[7] * t3;
  t93 = At[10] * t3;
  t107 = At[13] * t3;
  t114 =  ez * ez;
  t146 = t10 * t10;
  t149 = ((-5403237840000 * At[6] + (69470200800000 * At[4] + 19451656224000000 * At[0] * t4) * t3) * t10 * z + (-401283384 * At[15] + (4707059994 * At[13] + (3012121710 * At[12] + (-36011689560 * At[10] + (8027667648000 * At[9] + (-8027667648000 * At[7] + (-1296777081600000 * At[3] + 19451656224000000 * At[1] * t3) * t4) * t3) * z) * t3) * z) * t3 + ((-41423013450 * At[12] + (484040056500 * At[10] + (67540473000000 * At[6] - t38) * t4) * t3) * t44 + (1748257220 * At[15] + (-19964735910 * At[13] + (-2594411820000 * At[9] + (29331862560000 * At[7] + 3241942704000000 * At[3] * t4) * t3) * t4) * t3 + ((78248884350 * At[12] + (-860873013000 * At[10] + (-94556662200000 * At[6] + t38) * t4) * t3) * t44 + (-1938419560 * At[15] + (20997160275 * At[13] + (2283511230000 * At[9] - 21612951360000 * t70) * t4) * t3 + ((-47153256150 * At[12] + (459918459000 * At[10] + 27016189200000 * At[6] * t4) * t3) * t44 + (849093050 * At[15] + (-8397889500 * At[13] + (-643242600000 * At[9] + 3859455600000 * t70) * t4) * t3 + ((11448186750 * At[12] - 88445857500 * t93) * t44 + (-173573400 * At[15] + (1474097625 * At[13] + 53603550000 * At[9] * t4) * t3 + ((-1161410250 * At[12] + 5360355000 * t93) * t44 + (-113704500 * t107 + 17481100 * At[15] + (40608750 * At[12] * t44 + (-833000 * At[15] + 3123750 * t107 + 14875 * At[15] * t114) * ez) * ez) * ez) * ez) * ez) * ez) * ez) * ez) * ez) * ez) * ez) * ez) * ez) / ( C_pi * t146 * 0.58354968672000000e17 );

  return( t149 ) ; 


}


/******************************************************************************************/
/******************************************************************************************
 set_At():
 --------
    -- utility routine to set constant array used by   BesselJ_Debye_Eps_Exp():
 ******************************************************************************************/
void  set_At( double At[BESSEL_EPSILON_ORDER] )
{
  int m ; 
  double M; 

  for( m = 0 ; m < BESSEL_EPSILON_ORDER ; m++) { 
    M = (m+1)/3.;
    At[m] = sin(C_pi*M) * pow(6., M) * exp(lgamma(M)) ; 
  }
  return; 
}

/******************************************************************************************/
/******************************************************************************************
  BesselJ_Meissel_First():
  -----------------------
   -- Meissel's "First" Expansion, specified to higher-order by Chishtie et al. 2005. 
   -- good for  x<<n 
******************************************************************************************/

double BesselJ_Meissel_First( double n, double x ) 
{
  double U, Z, eps, ninv, fn, z, Vsum1, Vsum2, exp_val, factor, exp2, invZp1; 
  double t1, t2, t3, retval, loggamma_exp;

  fn =  n ; 
  z = x/fn;
  eps = (fn-x)/fn;
  Z = sqrt(eps*(1+z));
  ninv = 1./fn;
  U = 1./(n*Z*Z*Z);
  t1 = z * z;
  t2 =  ninv * ninv;

   
    // Part 1 of V_n sum 
  Vsum1  = (U * (860160 + 1290240 * t1 + ((-2580480 - 645120 * t1) * t1 + (-28672 + (2709504 + (6547968 + 672000 * t1) * t1) * t1 + ((-2580480 + (-23224320 + (-18708480 - 1048320 * t1) * t1) * t1) * t1 + (-8192 + (-2519040 + (-60518400 + (-151828480 + (-61254720 - 2163168 * t1) * t1) * t1) * t1) * t1 + ((2580480 + (138700800 + (800163840 + (940423680 + (228049920 + 5537280 * t1) * t1) * t1) * t1) * t1) * t1 + (6144 + (-2644992 + (-299351808 + (-3405435264 + (-8653594320 + (-5897669400 + (-954875250 - 16907985 * t1) * t1) * t1) * t1) * t1) * t1) * t1 + (2580480 + (625766400 + (12841758720 + (60631119360 + (86387857920 + (38435160960 + (4450158720 + 59968440 * t1) * t1) * t1) * t1) * t1) * t1) * t1) * t1 * U) * U) * U) * U) * U) * U) * U)) / 0.10321920e8;


  // Part 2 of V_n  sum 
  Vsum2 = -(ninv * (420 + (-14 + (-4 + 3 * t2) * t2) * t2)) / 0.5040e4;

  // I substitute  Gamma(n+1) with (n+1)*Gamma(n)  in the denominator:
  factor  = 1./( (fn+1.) * sqrt(Z) );
  

  // 1/(1+Z) to eigth order
  if( Z < 1.e-3 ) { 
    invZp1 = 1 + (-1 + (1 + (-1 + (1 + (-1 + (1 - Z) * Z) * Z) * Z) * Z) * Z) * Z;
  } 
  else { 
    invZp1 = 1./(1.+Z);
  }

  if( (eps < 1.e-4) && (fn > 1.e3)  ) { 
    t3 = t2*t2;
    loggamma_exp =  (ninv * (-420 + 14 * t2 - 4 * t3 + 3 * t3 * t2)) / 0.5040e4;  // O(1/n^k) part of -log(gamma(n))
    exp2 = -fn * sqrt(2.*eps) * eps * (0.984023040e9 + (0.442810368e9 + (0.303114240e9 + (0.233192960e9 + (0.190139040e9 + (0.160692840e9 + 0.139204065e9 * eps) * eps) * eps) * eps) * eps) * eps) / 0.1476034560e10;
    exp_val = 0.5*log(0.5*fn/C_pi) + loggamma_exp + exp2 - Vsum1 - Vsum2;
  }
  else { 
    exp_val = n*(  log(x*invZp1) - (1 - Z) ) - Vsum1 - Vsum2 - lgamma(n) ;
  }

  retval = exp_factor( factor, exp_val ) ;

  return( retval ) ; 
}
  
/******************************************************************************************/
/******************************************************************************************
  BesselJ_Meissel_Second():
  -----------------------
   -- Meissel's "Second" Expansion, specified to higher-order by Chishtie et al. 2005. 
   -- good for  x>>n 
******************************************************************************************/

double BesselJ_Meissel_Second( double n, double x ) 
{
  double U, Z, eps, fn, z, Qt, Qsum, exp_val, factor; 
  double t1, t2, retval;

  fn = n;
  z = x/fn;
  eps = (x-fn)/fn;
  Z = sqrt(eps*(1+z));
  U = 1./(n*Z*Z*Z);
  t1 = z * z;
  t2 = U * U;
   
    // P_n sum  ( exponent = -Psum  , so this is  (-Psum) )
  exp_val = (t1 * t2 * (-3072 - 768 * t1 + (3072 + (27648 + (22272 + 1248 * t1) * t1) * t1 + (-3072 + (-165120 + (-952576 + (-1119552 + (-271488 - 6592 * t1) * t1) * t1) * t1) * t1 + (3072 + (744960 + (15287808 + (72179904 + (102842688 + (45756144 + (5297808 + 71391 * t1) * t1) * t1) * t1) * t1) * t1) * t1) * t2) * t2) * t2)) / 0.12288e5;


  // Part 1 of Q_n  sum 
  Qt = n * ( Z - acosl(fn/x) );

  // Part 2 of Q_n  sum 
  Qsum = -(U * (860160 + 1290240 * t1 + (28672 + (-2709504 + (-6547968 - 672000 * t1) * t1) * t1 + (8192 + (2519040 + (60518400 + (151828480 + (61254720 + 2163168 * t1) * t1) * t1) * t1) * t1 + (-6144 + (2644992 + (299351808 + (3405435264 + (8653594320 + (5897669400 + (954875250 + 16907985 * t1) * t1) * t1) * t1) * t1) * t1) * t1) * t2) * t2) * t2)) / 0.10321920e8;


  factor = sqrt( 2 / ( C_pi*n*Z ) ) * cosl( Qsum + Qt - 0.25*C_pi )  ;

  retval = exp_factor( factor, exp_val ) ;

  return( retval ) ; 
}
  
/******************************************************************************************/
/******************************************************************************************
  BesselJ_Asympt1():
  ------------------
  -- valid for x/n << 1 and  n>>1   
  --  leading order terms are  O(1/n) and O(x/n);
 ******************************************************************************************/
double BesselJ_Asympt1( double n, double x)
{
  double z, factor, exp_val;

  z = x/n;
  factor = 1./sqrt(2*C_pi*(n+1));
  exp_val = n*( 1 + log(0.5*z) );
  return( exp_factor( factor, exp_val ) );
}  
  
/******************************************************************************************/
/******************************************************************************************
  BesselJ_bigx():
  ---------------
  -- valid for x/n >> 1 
  --  leading order terms are  O(1/x) 
 ******************************************************************************************/
double BesselJ_bigx( double n, double x)
{
  return(  sqrt(2./(C_pi*x)) * cosl( x - 0.5*C_pi*( n + 0.5 ) ) ) ;
}  
  
/******************************************************************************************/
/******************************************************************************************
    exp_factor():
    -----------------
       -- given f_exp, and f_factor, returns with   f_factor * exp( f_exp ) ;

       -- performs this calculation more accurately than the straightforward way;

*******************************************************************************************/
double exp_factor(  double f_factor, double f_exp )  
{
  double log_f, sign_f, fabs_exp, x;

  if( f_factor == 0. ) { return( 0. ); } 

  fabs_exp = fabs(f_exp);

  /* If it is small, then use a 9th-order Taylor expansion of the exponent : */
  if ( fabs_exp < 1.e-3 ) { 
    x = f_exp; 
    return(  f_factor 
             * ( 1. + ((40320.+(20160.+(6720.+(1680.+(336.+(56.+(8.+x)*x)*x)*x)*x)*x)*x)*x / 40320.) ) ) ;  
  }
  /* We want to prevent overflows, calculate exponential wisely if we are close: */
  else if( fabs_exp > 690. ) { 
    sign_f = ( f_factor < 0. ) ? -1. : 1. ;
    log_f = log( fabs(f_factor)  );
    if( log_f * f_exp < 0. ) { 
      return( sign_f * exp( log_f + f_exp ) ); 
    } 
    else { 
      return( f_factor * exp( f_exp ) );
    }
  }
  /* Standard expression: */
  else { 
    return( f_factor * exp( f_exp ) );
  }
}

#undef C_pi  
#undef BESSEL_EPSILON_ORDER
#undef SLOPE1 
#undef SLOPE2 
#undef SLOPE3 
#undef INTERCEPT1 
#undef INTERCEPT2 
#undef INTERCEPT3 
