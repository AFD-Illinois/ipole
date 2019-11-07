/*
 * restart.c
 *
 *  Created on: Sep 10, 2019
 *      Author: bprather
 */

#include "decs.h"
#include "hdf5_utils.h"
#include "model.h"
#include "model_radiation.h"
#include "par.h"

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

  int *tint = malloc(s2 * sizeof(*tint));
  double *tdbl = malloc(s2 * sizeof(*tdbl));

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
  int *tint = malloc(s2 * sizeof(*tint));
  double *tdbl = malloc(s2 * sizeof(*tdbl));

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
  hid_t dtype_version = hdf5_make_str_type(strlen(xstr(VERSION)));
  hdf5_add_attr(xstr(VERSION), "githash", "/", dtype_version);

  hdf5_make_directory("header");
  hdf5_set_directory("/header/");
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

  hdf5_set_directory("/header/");
  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&M_unit, "M_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Te_unit, "Thetae_unit", H5T_IEEE_F64LE);

  // allow model to output its parameters
  output_hdf5();
}

void dump(double image[], double imageS[], double taus[],
    const char *fname, double scale, double cam[NDIM],
    double fovx, double fovy, Params *params)
{
  hdf5_create(fname);

  write_header(scale, cam, fovx, fovy, params);

  int nx = params->nx;
  int ny = params->ny;
  double freqcgs = params->freqcgs;
  double dsource = params->dsource;

  // processing
  double Ftot_unpol=0., Ftot=0.;
  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      Ftot_unpol += image[i*ny+j] * scale;
      Ftot += imageS[(i*ny+j)*NIMG+0] * scale;
    }
  }

  hdf5_set_directory("/");
  // output stuff
  hdf5_write_single_val(&Ftot, "Ftot", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Ftot_unpol, "Ftot_unpol", H5T_IEEE_F64LE);
  double nuLnu = 4. * M_PI * Ftot * dsource * dsource * JY * freqcgs;
  hdf5_write_single_val(&nuLnu, "nuLnu", H5T_IEEE_F64LE);
  nuLnu = 4. * M_PI * Ftot_unpol * dsource * dsource * JY * freqcgs;
  hdf5_write_single_val(&nuLnu, "nuLnu_unpol", H5T_IEEE_F64LE);

  hsize_t unpol_dim[2] = {nx, ny};
  hsize_t pol_dim[3] = {nx, ny, NIMG};
  hdf5_write_full_array(image, "unpol", 2, unpol_dim, H5T_IEEE_F64LE);
  hdf5_write_full_array(taus, "tau", 2, unpol_dim, H5T_IEEE_F64LE);
  hdf5_write_full_array(imageS, "pol", 3, pol_dim, H5T_IEEE_F64LE);

  // housekeeping
  hdf5_close();
}

/*
 * Given a path, dump a variable computed along that path into a file.
 * Note this is most definitely *not* thread-safe
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

  double *X = calloc(NDIM*nstep, sizeof(double));
  double *K = calloc(NDIM*nstep, sizeof(double));
  double *Ucon = calloc(NDIM*nstep, sizeof(double));
  double *Ucov = calloc(NDIM*nstep, sizeof(double));
  double *Bcon = calloc(NDIM*nstep, sizeof(double));
  double *Bcov = calloc(NDIM*nstep, sizeof(double));

  double *j_inv = calloc(NDIM*nstep, sizeof(double));
  double *alpha_inv = calloc(NDIM*nstep, sizeof(double));
  double *rho_inv = calloc(NDIM*nstep, sizeof(double));

  for (int i=0; i<nstep; i++) {
    get_model_primitives(traj[i].X, &(prims[i*nprims]));
    get_model_fourv(traj[i].X, &(Ucon[i*NDIM]), &(Ucov[i*NDIM]), &(Bcon[i*NDIM]), &(Bcov[i*NDIM]));
    jar_calc(traj[i].X, traj[i].Kcon, &(j_inv[i*NDIM]), &(j_inv[i*NDIM+1]), &(j_inv[i*NDIM+2]), &(j_inv[i*NDIM+3]),
             &(alpha_inv[i*NDIM]), &(alpha_inv[i*NDIM+1]), &(alpha_inv[i*NDIM+2]), &(alpha_inv[i*NDIM+3]),
             &(rho_inv[i*NDIM+1]), &(rho_inv[i*NDIM+2]), &(rho_inv[i*NDIM+3]));

    b[i] = get_model_b(traj[i].X);
    ne[i] = get_model_ne(traj[i].X);
    thetae[i] = get_model_thetae(traj[i].X);

    MULOOP {
      X[i*NDIM+mu] = traj[i].X[mu];
      K[i*NDIM+mu] = traj[i].Kcon[mu];
    }
  }

  // PICTURES: Anything with one value for every pixel
  hsize_t fdims_p[] = { nx, ny };
  hsize_t fstart_p[] = { i, j };
  hsize_t fcount_p[] = { 1, 1 };
  hsize_t mdims_p[] =  { 1, 1 };
  hsize_t mstart_p[] = { 0, 0 };

  // Could just flat record the final emission values here...
  hdf5_write_array(&nstep, "nstep", 2, fdims_p, fstart_p, fcount_p, mdims_p, mstart_p, H5T_STD_I32LE);

  // SCALARS: Anything with one value for every geodesic step
  hsize_t fdims_s[] = { nx, ny, MAXNSTEP };
  hsize_t chunk_s[] =  { 1, 1, 200 };
  hsize_t fstart_s[] = { i, j, 0 };
  hsize_t fcount_s[] = { 1, 1, nstep };
  hsize_t mdims_s[] =  { 1, 1, nstep };
  hsize_t mstart_s[] = { 0, 0, 0 };

  hdf5_write_chunked_array(b, "b", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(ne, "ne", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(thetae, "thetae", 3, fdims_s, fstart_s, fcount_s, mdims_s, mstart_s, chunk_s, H5T_IEEE_F64LE);

  // VECTORS: Anything with N values per geodesic step
  hsize_t fdims_v[] = { nx, ny, MAXNSTEP, 4 };
  hsize_t chunk_v[] =  { 1, 1, 200, 4 };
  hsize_t fstart_v[] = { i, j, 0, 0 };
  hsize_t fcount_v[] = { 1, 1, nstep, 4 };
  hsize_t mdims_v[] =  { 1, 1, nstep, 4 };
  hsize_t mstart_v[] = { 0, 0, 0, 0 };

  hdf5_write_chunked_array(X, "X", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(K, "K", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);

  hdf5_write_chunked_array(Ucon, "Ucon", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(Ucov, "Ucov", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(Bcon, "Bcon", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(Bcov, "Bcov", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);

  hdf5_write_chunked_array(j_inv, "j_inv", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(alpha_inv, "alpha_inv", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);
  hdf5_write_chunked_array(rho_inv, "rho_inv", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);

  fdims_v[3] = chunk_v[3] = fcount_v[3] = mdims_v[3] = 8;
  hdf5_write_chunked_array(prims, "prims", 4, fdims_v, fstart_v, fcount_v, mdims_v, mstart_v, chunk_v, H5T_IEEE_F64LE);

  free(b); free(ne); free(thetae);

  free(X); free(K);
  free(Ucon); free(Ucov); free(Bcon); free(Bcov);
  free(j_inv); free(alpha_inv); free(rho_inv);

  free(prims);

  hdf5_close();
}
