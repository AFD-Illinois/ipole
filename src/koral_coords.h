#ifndef KORAL_COORDS
#define KORAL_COORDS

double minn(double a, double b, double df);
double maxx(double a, double b, double df);
double psi_smooth(double x);
double theta_smooth(double x);
double theta_disk_or_jet(double r, double x2, double rdecoll, double rcoll, double runi, double a1, double a2);
double wjet(double x2, double fdisk, double fjet);
double KORAL_theta_diskjet(double r, double x2, void *params);
double KORAL_cylindrify(double r, double x2, void *params);

typedef struct jetcoords_params_struct {
  double r_test, theta_test;
  double r0, rbrk, runi;
  double rdecoll_jet, rcoll_jet, rdecoll_disk, rcoll_disk;
  double alpha_1, alpha_2;
  double fdisk, fjet;
} jetcoords_params;

#endif // KORAL_COORDS
