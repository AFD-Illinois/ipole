/*
 * restart.c
 *
 *  Created on: Sep 10, 2019
 *      Author: bprather
 */

#include "decs.h"
#include "hdf5_utils.h"
#include "model.h"

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


void dump(double image[], double imageS[], double taus[],
    const char *fname, double scale, double Dsource, double cam[NDIM], double DX,
    double DY, double fovx, double fovy, double rcam, double thetacam, double phicam,
    double rotcam, double xoff, double yoff)
{
  hdf5_create(fname);
  hid_t dtype_version = hdf5_make_str_type(strlen(xstr(VERSION)));
  hdf5_add_attr(xstr(VERSION), "githash", "/", dtype_version);

  hdf5_make_directory("header");
  hdf5_set_directory("/header/");
  hdf5_write_single_val(&freqcgs, "freqcgs", H5T_IEEE_F64LE);
  hdf5_write_single_val(&scale, "scale", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Dsource, "dsource", H5T_IEEE_F64LE);

  char *conv;
  if (QU_CONVENTION == 0) {
    conv = "N";
  } else {
    conv = "W";
  }
  hdf5_write_single_val(conv, "evpa_0", H5T_C_S1);

  hdf5_make_directory("camera");
  hdf5_set_directory("/header/camera/");
  int nx = NX;
  int ny = NY;
  hdf5_write_single_val(&nx, "nx", H5T_STD_I32LE);
  hdf5_write_single_val(&ny, "ny", H5T_STD_I32LE);
  hdf5_write_single_val(&DX, "dx", H5T_IEEE_F64LE);
  hdf5_write_single_val(&DY, "dy", H5T_IEEE_F64LE);
  hdf5_write_single_val(&fovx, "fovx", H5T_IEEE_F64LE);
  hdf5_write_single_val(&fovy, "fovy", H5T_IEEE_F64LE);
  hdf5_write_single_val(&rcam, "rcam", H5T_IEEE_F64LE);
  hdf5_write_single_val(&thetacam, "thetacam", H5T_IEEE_F64LE);
  hdf5_write_single_val(&phicam, "phicam", H5T_IEEE_F64LE);
  hdf5_write_single_val(&rotcam, "rotcam", H5T_IEEE_F64LE);
  hdf5_write_single_val(&xoff, "xoff", H5T_IEEE_F64LE);
  hdf5_write_single_val(&yoff, "yoff", H5T_IEEE_F64LE);
  hsize_t vec_dim[1] = {NDIM};
  hdf5_write_full_array(cam, "x", 1, vec_dim, H5T_IEEE_F64LE);

  hdf5_set_directory("/header/");
  hdf5_make_directory("units");
  hdf5_set_directory("/header/units/");
  hdf5_write_single_val(&L_unit, "L_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&M_unit, "M_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&T_unit, "T_unit", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Te_unit, "Thetae_unit", H5T_IEEE_F64LE);

  // processing
  double Ftot_unpol=0., Ftot=0.;
  for (int i=0; i<NX; ++i) {
    for (int j=0; j<NY; ++j) {
      Ftot_unpol += image[i*NY+j] * scale;
      Ftot += imageS[(i*NY+j)*NIMG+0] * scale;
    }
  }

  hdf5_set_directory("/");
  // output stuff
  hdf5_write_single_val(&Ftot, "Ftot", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Ftot_unpol, "Ftot_unpol", H5T_IEEE_F64LE);
  double nuLnu = 4. * M_PI * Ftot * Dsource * Dsource * JY * freqcgs;
  hdf5_write_single_val(&nuLnu, "nuLnu", H5T_IEEE_F64LE);
  nuLnu = 4. * M_PI * Ftot_unpol * Dsource * Dsource * JY * freqcgs;
  hdf5_write_single_val(&nuLnu, "nuLnu_unpol", H5T_IEEE_F64LE);

  hsize_t unpol_dim[2] = {NX, NY};
  hsize_t pol_dim[3] = {NX, NY, NIMG};
  hdf5_write_full_array(image, "unpol", 2, unpol_dim, H5T_IEEE_F64LE);
  hdf5_write_full_array(taus, "tau", 2, unpol_dim, H5T_IEEE_F64LE);
  hdf5_write_full_array(imageS, "pol", 3, pol_dim, H5T_IEEE_F64LE);

  // allow model to output
  output_hdf5();

  // housekeeping
  hdf5_close();
}
