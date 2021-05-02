#include "koral_coords.h"
#include "simcoords.h"
#include "coordinates.h"

#include <math.h>
#include <assert.h>

#define TNY 96
#define RCYL 30.
#define NCYL 1.
#define MINY (0.0025)
#define MAXY (1.-0.0025)
#define RMIN (0.8*(1.+sqrt(1-a*a)))

static double x2cyl, rmidcyl, thetaCYL, sinthetaCYL, thetaAX, sinthetaAX;

double KORAL_theta_diskjet(double r, double x2, void *params);

int set_cyl_params()
{
   jetcoords_params tpar;
   x2cyl = MINY + 0.5*NCYL/((double)TNY);
   rmidcyl = 0.5 * (RCYL + RMIN);
    
   tpar.r0 = mp_koral_jetcoords.mksr0;
   tpar.rbrk = mp_koral_jetcoords.rbrk;
   tpar.fdisk = mp_koral_jetcoords.fdisk;
   tpar.fjet = mp_koral_jetcoords.fjet;
   tpar.runi = mp_koral_jetcoords.runi;
   tpar.rcoll_jet = mp_koral_jetcoords.rcoll_jet;
   tpar.rcoll_disk = mp_koral_jetcoords.rcoll_disk;
   tpar.rdecoll_jet = mp_koral_jetcoords.rdecoll_jet;
   tpar.rdecoll_disk = mp_koral_jetcoords.rdecoll_disk;
   tpar.alpha_1 = mp_koral_jetcoords.alpha_1;
   tpar.alpha_2 = mp_koral_jetcoords.alpha_2;

   thetaCYL = KORAL_theta_diskjet(RCYL, x2cyl, &tpar);
   sinthetaCYL = sin(thetaCYL);
   thetaAX = KORAL_theta_diskjet(RCYL, MAXY, &tpar);
   sinthetaAX = sin(thetaAX);
   
   return 0;
}

double sinth1(double r, double x2, void* params);
double sinth2(double r, double x2, void* params);

double sinth0(double r, double x2, void* params)
{
  double sinth0 = RCYL * sin(thetaCYL)/r;
  return sinth0;
}

double minn(double a, double b, double df)
{
  double delta = (b-a)/df;
  return b - psi_smooth(delta)*df;
}

double maxx(double a, double b, double df)
{
  return -1 * minn(-a, -b, df);
}

double psi_smooth(double x)
{
  if (x < -1) {
    return 0.;
  } else if(x>=1) {
    return x;
  } else {
    double xout = (-35.*cos(0.5*M_PI*x) - (5./6.)*cos(1.5*M_PI*x) + 0.1*cos(2.5*M_PI*x))/(32.*M_PI);
    return xout + 0.5*(x+1.);
  }
}

double theta_smooth(double x)
{
  if (x < -1) {
    return 0.;
  } else if (x >= 1) {
    return 1.;
  } else {
    return 0.5 + (70.*sin(0.5*M_PI*x) + 5*sin(1.5*M_PI*x) - sin(2.5*M_PI*x))/128.;
  }
}

double theta_disk_or_jet(double r, double x2, double rdecoll, double rcoll, double runi, double a1, double a2)
{
  double r1 = minn(r, rdecoll, 0.5*rdecoll)/runi;
  double r2 = minn(r/(r1*runi), rcoll/rdecoll, 0.5*rcoll/rdecoll);
  double y = pow(r2, a2)*tan(0.5*x2*M_PI);
  double x = pow(r1, a1); //opposite sign convention for alpha1 from ressler 2017!
  return 0.5*M_PI + atan2(y,x);
}

double wjet(double x2, double fdisk, double fjet)
{
  //NOTE! fjet and fdisk must both be positive and sum to < 1. 
  //NOTE! fjet is NOT defined as in Ressler 2017: their fjet = 1 - (our fjet)
  double delta = 2*(fabs(x2) - fdisk)/(1 - fjet - fdisk) - 1.;
  return theta_smooth(delta);
}

double KORAL_theta_diskjet(double r, double x2, void *params)
{
  jetcoords_params *par = (jetcoords_params *)params;
  double theta_disk, theta_jet, wfrac, theta;

  theta_disk = theta_disk_or_jet(r, x2, par->rdecoll_disk, par->rcoll_disk, par->runi,
			         par->alpha_1, par->alpha_2);
  theta_jet = theta_disk_or_jet(r, x2, par->rdecoll_jet, par->rcoll_jet, par->runi,
			        par->alpha_1, par->alpha_2);
  wfrac = wjet(x2, par->fdisk, par->fjet);
  theta = wfrac*theta_jet + (1-wfrac)*theta_disk;

  return theta;
}

double to1stquad(double x2)
{
  double ntimes = floor(0.25*(x2+2));
  double x2out = x2 - 4*ntimes;

  if (x2out > 0) {
    x2out = -x2out;
  }
  
  if (x2out < -1) {
    x2out = - 2 - x2out;
  }

  return x2out;
}

double f2func(double r, double x2, void* params)
{
  jetcoords_params *par = (jetcoords_params *)params;

  double s1in = sinth1(r, x2, par);
  double s2in = sinth2(r, x2, par);
  
  double s1ax = sinth1(r, MAXY, par); 
  double s2ax = sinth2(r, MAXY, par);
  double df = fabs(s2ax - s1ax) + 1.e-16; 

  if (r >= RCYL) {
    return maxx(s1in, s2in, df);
  } else {
    return minn(s1in, s2in, df);
  }
}

double sinth1(double r, double x2, void* params)
{
  jetcoords_params *par = (jetcoords_params *) params;
  double theta1 = KORAL_theta_diskjet(RCYL, x2, par);
  double sinth1 =  RCYL*sin(theta1)/r;
  return sinth1;
}

double sinth2(double r, double x2, void* params)
{
  jetcoords_params *par = (jetcoords_params *) params;

  double theta = KORAL_theta_diskjet(r, x2, par);
  double theta2 = KORAL_theta_diskjet(r, x2cyl, par);

  double thetamid = 0.5*M_PI;

  double thetaA = asin(sinth0(r,x2,par));
  double thetaB = (theta-theta2)*(thetamid-thetaA)/(thetamid-theta2);
  double sinth2 = sin(thetaA + thetaB);

  return sinth2;
}

double KORAL_cylindrify(double r, double x2, void *params)
{
  jetcoords_params *par = (jetcoords_params *) params;

  double thin = KORAL_theta_diskjet(r, x2, par);

  double x2mir = to1stquad(x2);
  double thmir = KORAL_theta_diskjet(r, x2mir, par);

  double f1 = sin(thmir);
  double f2 = f2func(r, x2mir, par);

  double thmid = KORAL_theta_diskjet(rmidcyl, x2mir, par);
  double f1mid = sin(thmid);
  double f2mid = f2func(rmidcyl, x2mir, par);

  double df = f2mid - f1mid;

  double thout = asin(maxx(r*f1, r*f2, r*fabs(df) + 1.e-16)/r);
  if (x2 != x2mir) {
    thout = thin + thmir - thout;
  }

  return thout;
}




void fwd_KORAL_JETCOORDS(double *xJET, double *xKS)
{
  // get coordinate system parameters
  double r0 = mp_koral_jetcoords.mksr0;
  double rbrk = mp_koral_jetcoords.rbrk;
  double alpha_1 = mp_koral_jetcoords.alpha_1;
  double alpha_2 = mp_koral_jetcoords.alpha_2;
  double fdisk = mp_koral_jetcoords.fdisk;
  double fjet = mp_koral_jetcoords.fjet;
  double runi = mp_koral_jetcoords.runi;
  double rcoll_jet = mp_koral_jetcoords.rcoll_jet;
  double rcoll_disk = mp_koral_jetcoords.rcoll_disk;
  double rdecoll_jet = mp_koral_jetcoords.rdecoll_jet;
  double rdecoll_disk = mp_koral_jetcoords.rdecoll_disk;

  /*
  double x1in = mp_koral_jetcoords.hypx1in;     // TODO
  double x1out = mp_koral_jetcoords.hypx1out;   // TODO
  double x1brk = mp_koral_jetcoords.hypx1brk;   // TODO
   */
  double x1in = 0.;
  double x1out = 0.;
  double x1brk = 0.;

  // so many parameters, so little time
  jetcoords_params tpar;
  tpar.r0 = r0;
  tpar.rbrk = rbrk;
  tpar.runi = runi;
  tpar.rcoll_jet = rcoll_jet;
  tpar.rcoll_disk = rcoll_disk;
  tpar.rdecoll_jet = rdecoll_jet;
  tpar.rdecoll_disk = rdecoll_disk;
  tpar.alpha_1 = alpha_1;
  tpar.alpha_2 = alpha_2;
  tpar.fdisk = fdisk;
  tpar.fjet = fjet;

  // pass time through
  xKS[0] = xJET[0];

  // get radial coordinate
  double x1sc = x1in + xJET[1] * (x1out-x1in);  // scale out of 0-1 range
  if (x1sc < x1brk) {
    xKS[1] = exp(x1sc) + r0;
  } else {
    // from hyperexp_func
    double x1brk = log(rbrk - r0);
    xKS[1] = r0 + exp(xJET[1] + 4. * pow(xJET[1] - x1brk, 4));
  }

  // get elevation coordinate
  if (mp_koral_jetcoords.cylindrify) {
    xKS[2] = KORAL_cylindrify(xKS[1], xJET[2], &tpar);
  } else {
    // return theta_diskjet(r, x2, par);
    assert(1==0);  // WARNING, THIS HAS NOT BEEN VALIDATED BY ITSELF!
    xKS[2] = KORAL_theta_diskjet(xKS[1], xJET[2], &tpar);
  }

  // pass phi through
  xKS[3] = xJET[3];
}


