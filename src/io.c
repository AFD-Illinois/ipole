/*
 * io.c
 *
 *  Created on: Sep 10, 2019
 *      Author: bprather
 */

#include "decs.h"
#include "coordinates.h"
#include "geometry.h"
#include "radiation.h"
#include "par.h"
#include "hdf5_utils.h"

#include "model.h"
#include "model_radiation.h"

#include <unistd.h>
#include <string.h>

void write_header(double scale, double cam[NDIM],
    double fovx, double fovy, Params *params);

void read_restart(const char *fname, double *tA, double *tB, double *last_img_target,
                  int *nopenimgs, int *nimg, int nconcurrentimgs, int s2,
                  double *target_times, int *valid_imgs, struct of_image *dimages) {

  hdf5_open(fname);

  hdf5_set_directory("/");
  hdf5_read_single_val(tA, "tA", H5T_IEEE_F64LE);
  hdf5_read_single_val(tB, "tB", H5T_IEEE_F64LE);
  hdf5_read_single_val(last_img_target, "last_img_target", H5T_IEEE_F64LE);
  hdf5_read_single_val(nimg, "nimg", H5T_STD_I32LE);
  hdf5_read_single_val(nopenimgs, "nopenimgs", H5T_STD_I32LE);

  int *tint = calloc(s2, sizeof(*tint));
  double *tdbl = calloc(s2, sizeof(*tdbl));

  for (int i=0; i<nconcurrentimgs; ++i) {
    tdbl[i] = target_times[i];
    tint[i] = valid_imgs[i];
  }

  hsize_t dims[1] = {nconcurrentimgs};
  hdf5_read_full_array(target_times, "target_times", 1, dims, H5T_IEEE_F64LE);
  hdf5_read_full_array(valid_imgs, "valid_images", 1, dims, H5T_STD_I32LE);

  dims[0] = s2;
  hdf5_set_directory("/dimg/");
  hdf5_read_full_array(tint, "nstep", 1, dims, H5T_STD_I32LE);
  for (int i=0; i<s2; ++i) dimages[i].nstep = tint[i];

  hdf5_read_full_array(tdbl, "intensity", 1, dims, H5T_IEEE_F64LE);
  for (int i=0; i<s2; ++i) dimages[i].intensity = tdbl[i];

  hdf5_read_full_array(tdbl, "tau", 1, dims, H5T_IEEE_F64LE);
  for (int i=0; i<s2; ++i) dimages[i].tau = tdbl[i];

  hdf5_read_full_array(tdbl, "tauF", 1, dims, H5T_IEEE_F64LE);
  for (int i=0; i<s2; ++i) dimages[i].tauF = tdbl[i];

  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      char tgt[20];
      snprintf(tgt, 18, "Nr%d%d", mu, nu);
      hdf5_read_full_array(tdbl, tgt, 1, dims, H5T_IEEE_F64LE);
      for (int i=0; i<s2; ++i) dimages[i].N_coord[mu][nu] = tdbl[i];
      snprintf(tgt, 18, "Ni%d%d", mu, nu);
      hdf5_read_full_array(tdbl, tgt, 1, dims, H5T_IEEE_F64LE);
      for (int i=0; i<s2; ++i) dimages[i].N_coord[mu][nu] += tdbl[i] * _Complex_I;
    }
  }

  free(tint);
  free(tdbl);

  hdf5_close();

}

void write_restart(const char *fname, double tA, double tB, double last_img_target,
        int nopenimgs, int nimg, int nconcurrentimgs, int s2,
        double *target_times, int *valid_imgs, struct of_image *dimages) {

  hdf5_create(fname);

  hid_t dtype_version = hdf5_make_str_type(strlen(xstr(VERSION)));
  hdf5_add_attr(xstr(VERSION), "githash", "/", dtype_version);
  hdf5_write_single_val(&tA, "/tA", H5T_IEEE_F64LE);
  hdf5_write_single_val(&tB, "/tB", H5T_IEEE_F64LE);
  hdf5_write_single_val(&last_img_target, "/last_img_target", H5T_STD_I32LE);
  hdf5_write_single_val(&nimg, "/nimg", H5T_STD_I32LE);
  hdf5_write_single_val(&nopenimgs, "/nopenimgs", H5T_STD_I32LE);

  hsize_t dims[1] = {nconcurrentimgs};
  hdf5_write_full_array(target_times, "/target_times", 1, dims, H5T_IEEE_F64LE);
  hdf5_write_full_array(valid_imgs, "/valid_images", 1, dims, H5T_STD_I32LE);

  dims[0] = s2;
  hdf5_make_directory("dimg");
  hdf5_set_directory("/dimg/");
  // save dimg struct
  int *tint = calloc(s2, sizeof(*tint));
  double *tdbl = calloc(s2, sizeof(*tdbl));

  for (int i=0; i<s2; ++i) tint[i] = dimages[i].nstep;
  hdf5_write_full_array(tint, "nstep", 1, dims, H5T_STD_I32LE);

  for (int i=0; i<s2; ++i) tdbl[i] = dimages[i].intensity;
  hdf5_write_full_array(tdbl, "intensity", 1, dims, H5T_IEEE_F64LE);

  for (int i=0; i<s2; ++i) tdbl[i] = dimages[i].tau;
  hdf5_write_full_array(tdbl, "tau", 1, dims, H5T_IEEE_F64LE);

  for (int i=0; i<s2; ++i) tdbl[i] = dimages[i].tauF;
  hdf5_write_full_array(tdbl, "tauF", 1, dims, H5T_IEEE_F64LE);

  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      char tgt[20];
      snprintf(tgt, 18, "Nr%d%d", mu, nu);
      for (int i=0; i<s2; ++i) tdbl[i] = creal(dimages[i].N_coord[mu][nu]);
      hdf5_write_full_array(tdbl, tgt, 1, dims, H5T_IEEE_F64LE);
      snprintf(tgt, 18, "Ni%d%d", mu, nu);
      for (int i=0; i<s2; ++i) tdbl[i] = cimag(dimages[i].N_coord[mu][nu]);
      hdf5_write_full_array(tdbl, tgt, 1, dims, H5T_IEEE_F64LE);
    }
  }

  free(tint);
  free(tdbl);
  hdf5_close();
}

void write_header(double scale, double cam[NDIM],
    double fovx, double fovy, Params *params)
{
  hid_t dtype_version = hdf5_make_str_type(20);
  hdf5_add_attr(xstr(VERSION), "githash", "/", dtype_version);

  hdf5_make_directory("header");
  hdf5_set_directory("/header/");
  hdf5_write_single_val(VERSION_STRING, "version", dtype_version);
  hdf5_write_single_val(xstr(VERSION), "githash", dtype_version);
  // Make locals for things we'll use later
  double freqcgs = params->freqcgs;
  double dsource = params->dsource;
  hdf5_write_single_val(&freqcgs, "freqcgs", H5T_IEEE_F64LE);
  hdf5_write_single_val(&scale, "scale", H5T_IEEE_F64LE);
  hdf5_write_single_val(&dsource, "dsource", H5T_IEEE_F64LE);

  char conv[2];
  if (params->qu_conv == 0) {
    strncpy(conv,"N",2);
  } else {
    strncpy(conv,"W",2);
  }
  hid_t dtype_evpa = hdf5_make_str_type(2);
  hdf5_write_single_val(conv, "evpa_0", dtype_evpa);

  hdf5_make_directory("camera");
  hdf5_set_directory("/header/camera/");
  int nx = params->nx;
  int ny = params->ny;
  hdf5_write_single_val(&nx, "nx", H5T_STD_I32LE);
  hdf5_write_single_val(&ny, "ny", H5T_STD_I32LE);
  hdf5_write_single_val(&(params->dx), "dx", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(params->dy), "dy", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(params->fovx_dsource), "fovx_dsource", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(params->fovy_dsource), "fovy_dsource", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(params->rcam), "rcam", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(params->thetacam), "thetacam", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(params->phicam), "phicam", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(params->rotcam), "rotcam", H5T_IEEE_F64LE);
  hdf5_write_single_val(&fovx, "fovx", H5T_IEEE_F64LE);
  hdf5_write_single_val(&fovy, "fovy", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(params->xoff), "xoff", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(params->yoff), "yoff", H5T_IEEE_F64LE);
  hsize_t vec_dim[1] = {NDIM};
  hdf5_write_full_array(cam, "x", 1, vec_dim, H5T_IEEE_F64LE);

  //fprintf(stderr, "Wrote header\n");

  // allow model to output its parameters
  output_hdf5();
}

void dump_check(double image[], double imageS[], double taus[], Params *params, int skip)
{
  int unpol = params->only_unpolarized;

  size_t nx = params->nx;
  size_t ny = params->ny;

  if (access(params->outf, F_OK) != -1) {
    hdf5_append(params->outf);
  } else {
    fprintf(stderr, "oh no :(\nno file to check.\n\n");
    return;
  }

  hsize_t unpol_mdim[2] = {nx, ny};
  hsize_t pol_mdim[3] = {nx, ny, NIMG};
  hsize_t unpol_fdim[2] = { (int)ceil(nx/skip), (int)ceil(ny/skip) };
  hsize_t pol_fdim[3] = { (int)ceil(nx/skip), (int)ceil(ny/skip), NIMG };
  hsize_t unpol_start[2] = {0, 0};
  hsize_t pol_start[3] = {0, 0, 0};

  hdf5_write_array(image, "unpol_check", 2, unpol_fdim, unpol_start, unpol_fdim, unpol_mdim, unpol_start, H5T_IEEE_F64LE);
  hdf5_write_array(taus, "tau_check", 2, unpol_fdim, unpol_start, unpol_fdim, unpol_mdim, unpol_start, H5T_IEEE_F64LE);
  if (!unpol) hdf5_write_array(imageS, "pol_check", 3, pol_fdim, pol_start, pol_fdim, pol_mdim, pol_start, H5T_IEEE_F64LE);

  hsize_t unpol_mstart[2] = { (int)ceil(nx/skip), 0 };
  hsize_t pol_mstart[3] = { unpol_mstart[0], 0, 0 };

  hdf5_write_array(image, "unpol_recheck", 2, unpol_fdim, unpol_start, unpol_fdim, unpol_mdim, unpol_mstart, H5T_IEEE_F64LE);
  hdf5_write_array(taus, "tau_recheck", 2, unpol_fdim, unpol_start, unpol_fdim, unpol_mdim, unpol_mstart, H5T_IEEE_F64LE);
  if (!unpol) hdf5_write_array(imageS, "pol_recheck", 3, pol_fdim, pol_start, pol_fdim, pol_mdim, pol_mstart, H5T_IEEE_F64LE);

  // housekeeping
  hdf5_close();
}

void dump(double image[], double imageS[], double taus[],
    const char *fname, double scale, double cam[NDIM],
    double fovx, double fovy, size_t nx, size_t ny, Params *params)
{
  int unpol = params->only_unpolarized;

  hdf5_create(fname);

  write_header(scale, cam, fovx, fovy, params);

  double freqcgs = params->freqcgs;
  double dsource = params->dsource;

  // processing
  double Ftot_unpol=0., Ftot=0.;
  for (int i=0; i < params->nx; ++i) {
    for (int j=0; j < params->ny; ++j) {
      Ftot_unpol += image[i*ny+j] * scale;
      if (!unpol) Ftot += imageS[(i*ny+j)*NIMG+0] * scale;
    }
  }

  hdf5_set_directory("/");
  // output stuff
  hdf5_write_single_val(&Ftot_unpol, "Ftot_unpol", H5T_IEEE_F64LE);
  if (!unpol) hdf5_write_single_val(&Ftot, "Ftot", H5T_IEEE_F64LE);
  double nuLnu = 4. * M_PI * Ftot * dsource * dsource * JY * freqcgs;
  hdf5_write_single_val(&nuLnu, "nuLnu", H5T_IEEE_F64LE);
  nuLnu = 4. * M_PI * Ftot_unpol * dsource * dsource * JY * freqcgs;
  hdf5_write_single_val(&nuLnu, "nuLnu_unpol", H5T_IEEE_F64LE);

  hsize_t unpol_mdim[2] = {nx, ny};
  hsize_t pol_mdim[3] = {nx, ny, NIMG};
  hsize_t unpol_fdim[2] = {params->nx, params->ny};
  hsize_t pol_fdim[3] = {params->nx, params->ny, NIMG};
  hsize_t unpol_start[2] = {0, 0};
  hsize_t pol_start[3] = {0, 0, 0};

  hdf5_write_array(image, "unpol", 2, unpol_fdim, unpol_start, unpol_fdim, unpol_mdim, unpol_start, H5T_IEEE_F64LE);
  hdf5_write_array(taus, "tau", 2, unpol_fdim, unpol_start, unpol_fdim, unpol_mdim, unpol_start, H5T_IEEE_F64LE);
  if (!unpol) hdf5_write_array(imageS, "pol", 3, pol_fdim, pol_start, pol_fdim, pol_mdim, pol_start, H5T_IEEE_F64LE);

  // housekeeping
  hdf5_close();
}

/*
 * Given a path, dump a variable computed along that path into a file.
 * Note this is most definitely *not* thread-safe, so it gets called from an 'omp critical'
 * 
 * TODO output of I,Q,U,V, as parallel-transported to camera, along traces
 */
void dump_var_along(int i, int j, int nstep, struct of_traj *traj, int nx, int ny,
                    double scale, double cam[NDIM], double fovx, double fovy, Params *params)
{
  if (access(params->trace_outf, F_OK) != -1) {
    hdf5_append(params->trace_outf);
  } else {
    hdf5_create(params->trace_outf);
    write_header(scale, cam, fovx, fovy, params);
  }

  int nprims = 8;
  double *prims = calloc(nprims*nstep, sizeof(double));

  double *b = calloc(nstep, sizeof(double));
  double *ne = calloc(nstep, sizeof(double));
  double *thetae = calloc(nstep, sizeof(double));
  double *nu = calloc(nstep, sizeof(double));
  double *mu = calloc(nstep, sizeof(double));

  double *dl = calloc(nstep, sizeof(double));
  double *X = calloc(NDIM*nstep, sizeof(double));
  // Vectors aren't really needed.  Put back in behind flag if it comes up
  double *Kcon = calloc(NDIM*nstep, sizeof(double));
//  double *Kcov = calloc(NDIM*nstep, sizeof(double));

  double *r = calloc(nstep, sizeof(double));
  double *th = calloc(nstep, sizeof(double));
  double *phi = calloc(nstep, sizeof(double));

//  double *Ucon = calloc(NDIM*nstep, sizeof(double));
//  double *Ucov = calloc(NDIM*nstep, sizeof(double));
//  double *Bcon = calloc(NDIM*nstep, sizeof(double));
//  double *Bcov = calloc(NDIM*nstep, sizeof(double));

  // TODO NDIM and NSTOKES are not the same thing
  double *j_inv = calloc(NDIM*nstep, sizeof(double));
  double *alpha_inv = calloc(NDIM*nstep, sizeof(double));
  double *rho_inv = calloc(NDIM*nstep, sizeof(double));
  double *unpol_inv = calloc(2*nstep, sizeof(double));

  for (int i=0; i<nstep; i++) {
    get_model_primitives(traj[i].X, &(prims[i*nprims]));
    double Ucont[NDIM], Ucovt[NDIM], Bcont[NDIM], Bcovt[NDIM];
    get_model_fourv(traj[i].X, traj[i].Kcon, Ucont, Ucovt, Bcont, Bcovt);
    jar_calc(traj[i].X, traj[i].Kcon, &(j_inv[i*NDIM]), &(j_inv[i*NDIM+1]), &(j_inv[i*NDIM+2]), &(j_inv[i*NDIM+3]),
             &(alpha_inv[i*NDIM]), &(alpha_inv[i*NDIM+1]), &(alpha_inv[i*NDIM+2]), &(alpha_inv[i*NDIM+3]),
             &(rho_inv[i*NDIM+1]), &(rho_inv[i*NDIM+2]), &(rho_inv[i*NDIM+3]), params);
    get_jkinv(traj[i].X, traj[i].Kcon, &(unpol_inv[i*2+0]), &(unpol_inv[i*2+1]), params);

    b[i] = get_model_b(traj[i].X);
    ne[i] = get_model_ne(traj[i].X);
    thetae[i] = get_model_thetae(traj[i].X);
    nu[i] = get_fluid_nu(traj[i].Kcon, Ucovt);
    mu[i] = get_bk_angle(traj[i].X, traj[i].Kcon, Ucovt, Bcont, Bcovt);

    dl[i] = traj[i].dl;
    MULOOP {
      X[i*NDIM+mu] = traj[i].X[mu];
      Kcon[i*NDIM+mu] = traj[i].Kcon[mu];
    }
//    double Gcov[NDIM][NDIM];
//    gcov_func(traj[i].X, Gcov);
//    lower(traj[i].Kcon, Gcov, &(Kcov[i*NDIM]));

    bl_coord(traj[i].X, &(r[i]), &(th[i]));
    phi[i] = traj[i].X[3];
  }

  hdf5_set_directory("/");

  // PICTURES: Anything with one value for every pixel
  hsize_t fdims_p[] = { nx, ny };
  hsize_t fstart_p[] = { i, j };
  hsize_t fcount_p[] = { 1, 1 };
  hsize_t mdims_p[] =  { 1, 1 };
  hsize_t mstart_p[] = { 0, 0 };

  // Could just flat record the final emission values here...
  hdf5_write_array(&nstep, "nstep", 2, fdims_p, fstart_p, fcount_p, mdims_p, mstart_p, H5T_STD_I32LE);

  // SCALARS: Anything with one value for every geodesic step
  hsize_t fdims_s[] = { nx, ny, params->maxnstep };
  hsize_t chunk_s[] =  { 1, 1, 200 };
  hsize_t fstart_s[] = { i, j, 0 };
  hsize_t fcount_s[] = { 1, 1, nstep };
  hsize_t mdims_s[] =  { 1, 1, nstep };
  hsize_t mstart_s[] = { 0, 0, 0 };

  hdf5_write_chunked_array(b, "b", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(ne, "ne", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(thetae, "thetae", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(nu, "nu", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(mu, "b_k_angle", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);

  hdf5_write_chunked_array(dl, "dl", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(r, "r", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(th, "th", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(phi, "phi", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);

  // VECTORS: Anything with N values per geodesic step
  hsize_t fdims_v[] = { nx, ny, params->maxnstep, 4 };
  hsize_t chunk_v[] =  { 1, 1, 200, 4 };
  hsize_t fstart_v[] = { i, j, 0, 0 };
  hsize_t fcount_v[] = { 1, 1, nstep, 4 };
  hsize_t mdims_v[] =  { 1, 1, nstep, 4 };
  hsize_t mstart_v[] = { 0, 0, 0, 0 };

  hdf5_write_chunked_array(X, "X", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(Kcon, "Kcon", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
//  hdf5_write_chunked_array(Kcov, "Kcov", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);

//  hdf5_write_chunked_array(Ucon, "Ucon", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
//  hdf5_write_chunked_array(Ucov, "Ucov", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
//  hdf5_write_chunked_array(Bcon, "Bcon", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
//  hdf5_write_chunked_array(Bcov, "Bcov", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);

  hdf5_write_chunked_array(j_inv, "j_inv", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(alpha_inv, "alpha_inv", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(rho_inv, "rho_inv", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);

  fdims_v[3] = chunk_v[3] = fcount_v[3] = mdims_v[3] = 8;
  hdf5_write_chunked_array(prims, "prims", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);

  // Unpolarized coefficients j and k
  fdims_v[3] = chunk_v[3] = fcount_v[3] = mdims_v[3] = 2;
  hdf5_write_chunked_array(unpol_inv, "jk", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);

  free(b); free(ne); free(thetae);
  free(nu); free(mu);

  free(X);
  //free(Kcon); free(Kcov);
  free(r); free(th); free(phi);

  //free(Ucon); free(Ucov); free(Bcon); free(Bcov);
  free(j_inv); free(alpha_inv); free(rho_inv);
  free(unpol_inv);

  free(prims);

  hdf5_close();
}
