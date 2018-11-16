import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.ndimage.interpolation import rotate

SMALL = 1.e-30

if len(sys.argv) < 3 :
	print "usage: ipole.py ipole.dat phi(deg)"
	quit()

fil = sys.argv[1]
phi = float(sys.argv[2])

# read in header
with open(fil, 'r') as f:
  line = f.readline().split(' ')
  NX = int(line[0])
  NY = int(line[1])
  DX = float(line[2])
  DY = float(line[3])
  scale = float(line[4])
  L_unit = float(line[5])
  M_unit = float(line[6])

print 'NX = %i' % NX
print 'NY = %i' % NY

# read in data 
i0, j0, x, y, Ia, Is, Qs, Us, Vs, tauF = np.loadtxt(fil, unpack=True, skiprows=1)

# set image size
#ImRes = int(round(np.sqrt(len(i0))))
#print "Image resolution: ", ImRes

# set render size
FOV = 40.

# size of single pixel in rad: M/pixel . muas/pix . rad/muas
FOV = 40.
print "rendering FOV = ", FOV, "GM/c^2"
#
#da = FOV/ImRes * 5. / (1.e6 * 206265.)
# solid angle subtended by pixel
#dO = da*da
#print dO
#print scale
#print sum(Is)*scale
#sodifj
#Jy = 1.e-23  # cgs
#flux = dO*sum(Is)/Jy
flux = sum(Is)*scale
print "Flux [Jy]: ", flux

# recast indices into offset in units of M
i = (np.reshape(i0, (NX,NY))+1)*DX/NX - DY/2
j = (np.reshape(j0, (NX,NY))+1)*DY/NY - DY/2

# LP plot
ax = plt.subplot(2,2,2)
lpfrac = 100.*np.sqrt(Qs*Qs + Us*Us)/(Is + SMALL)
z = rotate(np.reshape(lpfrac, (NX,NY)), phi, reshape=False)
plt.pcolormesh(i,j,z,cmap='hot', vmin = 0., vmax = 100.)
plt.title('LP [%]')
plt.axis([-20,20,-20,20])
plt.colorbar()
ax.set_aspect('equal')

# EVPA plot
ax = plt.subplot(2,2,3)
evpa = (180./3.14159)*0.5*np.arctan2(Us,Qs)
z = rotate(np.reshape(evpa, (NX,NY)), phi, reshape=False)
plt.pcolormesh(i,j,z,cmap='jet')
plt.title('EVPA [deg]')
plt.axis([-20,20,-20,20])
plt.colorbar()
ax.set_aspect('equal')

# CP plot
ax = plt.subplot(2,2,4)
cpfrac = 100.*Vs/(Is + SMALL)
z = rotate(np.reshape(cpfrac, (NX,NY)), phi, reshape=False)
plt.pcolormesh(i,j,z,cmap='jet', vmin = -5, vmax = 5.)
plt.title('CP [%]')
plt.axis([-20,20,-20,20])
plt.colorbar()
ax.set_aspect('equal')

# total intensity 
ax = plt.subplot(2,2,1)
z = rotate(np.reshape(Is, (NX,NY)), phi, reshape=False)
plt.pcolormesh(i,j,z,cmap='afmhot', vmin=0, vmax=z.max())
#plt.pcolormesh(i,j,np.log10(z/z.max()),cmap='afmhot', vmin=-4, vmax=0.)
plt.colorbar()
plt.title('Stokes I [cgs]')
plt.axis([-20,20,-20,20])
ax.set_aspect('equal')

# superpose EV on total intensity using quiver
Qb = sum(Qs)
Ub = sum(Us)
LP = np.sqrt(Qb*Qb + Ub*Ub)/sum(Is)
print "LP [%]: ", LP*100.
CP = sum(Vs)/sum(Is)
print "CP [%]: ", CP*100.
amp = np.sqrt(Qs*Qs + Us*Us)
scal = max(amp)
#print(scal)
vxp = np.sqrt(Qs*Qs + Us*Us)*np.cos(evpa*3.14159/180.)/scal
vyp = np.sqrt(Qs*Qs + Us*Us)*np.sin(evpa*3.14159/180.)/scal
vx = rotate(np.reshape(vxp, (NX,NY)), phi, reshape=False)
vy = rotate(np.reshape(vyp, (NX,NY)), phi, reshape=False)
skip = 16*NX/256
plt.quiver(i[::skip, ::skip],j[::skip, ::skip],vx[::skip, ::skip],vy[::skip, ::skip], 
	headwidth=1, headlength=1, 
	width=0.005,
	color='green', 
	units='width', 
	scale=16)

plt.subplots_adjust(wspace=0.3,hspace=0.3)

# show, or save
plt.show()
#plt.savefig('tst')

