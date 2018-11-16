import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) is not 3:
	print "usage: ipole.py ipole.dat filenum"
	quit()

fil = sys.argv[1]
n = int(sys.argv[2])

# read in data 
i0, j0, x, y, Ia, Is, Qs, Us, Vs, tauF = np.loadtxt(fil, unpack=True, skiprows=1)

# set image size
ImRes = int(round(np.sqrt(len(i0))))
print "Image resolution: ", ImRes

# size of single pixel in rad: M/pixel . muas/pix . rad/muas
FOV = 40.
print "assuming FOV = ", FOV, "GM/c^2"
#
da = FOV/ImRes * 5. / (1.e6 * 206265.)
# solid angle subtended by pixel
dO = da*da
Jy = 1.e-23  # cgs
flux = dO*sum(Is)/Jy
print "Flux [Jy]: ", flux

# recast indices into offset in units of M
i = (np.reshape(i0, (ImRes,ImRes))+1)*40./ImRes - 20.
j = (np.reshape(j0, (ImRes,ImRes))+1)*40./ImRes - 20.

fig = plt.figure(figsize=(9,7))

# LP plot
ax = plt.subplot(2,2,2)
lpfrac = 100.*np.sqrt(Qs*Qs + Us*Us)/Is
z = np.reshape(lpfrac, (ImRes,ImRes))
plt.pcolormesh(i,j,z,cmap='jet', vmin = 0., vmax = 100.)
plt.title('LP [%]')
plt.axis([-20,20,-20,20])
plt.colorbar()
ax.set_aspect('equal')

# EVPA plot
ax = plt.subplot(2,2,3)
evpa = (180./3.14159)*0.5*np.arctan2(Us,Qs)
z = np.reshape(evpa, (ImRes,ImRes))
plt.pcolormesh(i,j,z,cmap='jet', vmin=-90,vmax=90)
plt.title('EVPA [deg]')
plt.axis([-20,20,-20,20])
plt.colorbar()
ax.set_aspect('equal')

# tauF plot
ax = plt.subplot(2,2,4)
z = np.reshape(tauF, (ImRes,ImRes))
plt.pcolormesh(i,j,np.log10(z),cmap='RdBu', vmin = 0, vmax = 4.)
plt.title('log10 tauF')
plt.axis([-20,20,-20,20])
plt.colorbar()
ax.set_aspect('equal')

# total intensity 
ax = plt.subplot(2,2,1)
z = np.reshape(Is, (ImRes,ImRes))
plt.pcolormesh(i,j,z,cmap='afmhot', vmin=0., vmax=0.0007)
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
vx = np.reshape(vxp, (ImRes,ImRes))
vy = np.reshape(vyp, (ImRes,ImRes))
skip = 16
plt.quiver(i[::skip, ::skip],j[::skip, ::skip],vx[::skip, ::skip],vy[::skip, ::skip], 
	headwidth=1, headlength=1, 
	width=0.005,
	color='green', 
	units='width', 
	scale=16)

plt.subplots_adjust(wspace=0.3,hspace=0.3)

print 'about to save'
plt.suptitle('t = %04d' % (5*n))

# show, or save
#plt.show()
plt.savefig('frame_%08d.png' % n, bbox_inches='tight')

