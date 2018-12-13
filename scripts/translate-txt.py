
import sys
import h5py


def write_vtk(fname, hfp):
  """ """

  print("making {0:s}".format(fname))

  imgu = hfp['unpol'][:,:]
  imgp = hfp['pol'][:,:,:]
  nx = imgu.shape[0]
  ny = imgu.shape[1]

  fp = open(fname,'w')

  for j in range(ny):
    for i in range(nx):
      fp.write("{0:d} {1:d} {2:g} {3:g} {4:g} {5:g} {6:g}\n".format(i,j, imgu[i,j], imgp[i,j,0], imgp[i,j,1], imgp[i,j,2], imgp[i,j,3]))

  fp.close()


if __name__ == "__main__":

  fnames = [ x for x in sys.argv[1:] if x[-3:] == ".h5" ]

  for fname in fnames:

    ofname = fname.replace(".h5",".dat")

    hfp = h5py.File(fname,'r')
    write_vtk(ofname, hfp)
    hfp.close()

