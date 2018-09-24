
#include "h5io.h"

void h5io_add_group(hid_t fid, const char *path) 
{
  hid_t group_id = H5Gcreate2(fid, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);
}

void h5io_add_attribute_int(hid_t fid, const char *path, const char *name, int attribute)
{
  hid_t dataspace_id = H5Screate(H5S_SCALAR);
  hid_t attribute_id = H5Acreate_by_name(fid, path, name, H5IO_FMT_INT, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_INT, &attribute);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);
}

void h5io_add_attribute_dbl(hid_t fid, const char *path, const char *name, double attribute)
{
  hid_t dataspace_id = H5Screate(H5S_SCALAR);
  hid_t attribute_id = H5Acreate_by_name(fid, path, name, H5IO_FMT_DBL, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &attribute);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);
}

void h5io_add_attribute_str(hid_t fid, const char *path, const char *name, const char *attribute)
{ 
  hid_t dtype_id = H5Tcopy(H5IO_FMT_STR);
  H5Tset_size(dtype_id, strlen(attribute)+1);
  H5Tset_strpad(dtype_id, H5T_STR_NULLTERM);
  hid_t dataspace_id = H5Screate(H5S_SCALAR);
  hid_t attribute_id = H5Acreate_by_name(fid, path, name, dtype_id, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, dtype_id, attribute);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);
  H5Tclose(dtype_id);
}

void h5io_add_data_int(hid_t fid, const char *path, int data) 
{
  hid_t dataspace_id = H5Screate(H5S_SCALAR);
  hid_t dataset_id = H5Dcreate2(fid, path, H5IO_FMT_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}


void h5io_add_data_dbl(hid_t fid, const char *path, double data) 
{
  hid_t dataspace_id = H5Screate(H5S_SCALAR);
  hid_t dataset_id = H5Dcreate2(fid, path, H5IO_FMT_DBL, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}

void h5io_add_data_str(hid_t fid, const char *path, const char *data)
{
  hid_t dtype_id = H5Tcopy(H5IO_FMT_STR);
  H5Tset_size(dtype_id, strlen(data)+1);
  H5Tset_strpad(dtype_id, H5T_STR_NULLTERM);
  hid_t dataspace_id = H5Screate(H5S_SCALAR);
  hid_t dataset_id = H5Dcreate2(fid, path, dtype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Tclose(dtype_id);
}

void h5io_add_data_dbl_1d(hid_t fid, const char *path, hsize_t n1, double data[n1])
{
  hsize_t dims[1] = { n1 };
  hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
  hid_t dataset_id = H5Dcreate2(fid, path, H5IO_FMT_DBL, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}

void h5io_add_data_dbl_2d(hid_t fid, const char *path, hsize_t n1, hsize_t n2, double data[n1][n2])
{
  hsize_t dims[2] = { n1, n2 };
  hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
  hid_t dataset_id = H5Dcreate2(fid, path, H5IO_FMT_DBL, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id); 
}

void h5io_add_data_dbl_3d(hid_t fid, const char *path, hsize_t n1, hsize_t n2, hsize_t n3, double data[n1][n2][n3])
{
  hsize_t dims[3] = { n1, n2, n3 };
  hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
  hid_t dataset_id = H5Dcreate2(fid, path, H5IO_FMT_DBL, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id); 
}

