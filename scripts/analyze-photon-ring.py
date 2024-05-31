import numpy as np 
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt 
import pandas as pd 
import h5py
import csv
import sys
from numpy import unravel_index
from scipy import ndimage, datasets
from scipy.ndimage import gaussian_filter

#GO TO LINE 166 IN ORDER TO CHANGE THE THEORY LABEL BEFORE RUNNING!

if __name__ == "__main__":
  allRuns=[]
  for fname in sys.argv[1:]:
    
    if fname[-3:] != ".h5": continue
    print("plotting {0:s}".format(fname))

    do_fullscreen_unpol = False
    if '--fullscreen' in sys.argv:
      do_fullscreen_unpol = True

    # load
    hfp = h5py.File(fname,'r')    
    dx = hfp['header']['camera']['dx'][()]
    dsource = hfp['header']['dsource'][()]
    lunit = hfp['header']['units']['L_unit'][()]
    fov_muas = dx / dsource * lunit * 2.06265e11
    scale = hfp['header']['scale'][()]
    evpa_0 = 'W'
    if 'evpa_0' in hfp['header']:
      evpa_0 = hfp['header']['evpa_0'][()]
    unpol = np.copy(hfp['unpol']).transpose((1,0))
    if 'pol' in hfp:
      no_pol = False
      do_fullscreen_unpol = False
      imagep = np.copy(hfp['pol']).transpose((1,0,2))
      I = imagep[:,:,0]
      Q = imagep[:,:,1]
      U = imagep[:,:,2]
      V = imagep[:,:,3]
    else:
      no_pol = True
      I = unpol
      Q = np.zeros((1, 1))
      U = np.zeros((1, 1))
      V = np.zeros((1, 1))
    hfp.close()

    npix = I.shape[0]
    xs = np.linspace(-fov_muas/2,fov_muas/2,npix)
    Xs,Ys = np.meshgrid(xs,xs)

    #Get diameter
    diameter=0
    #if care about diameter positions uncomment this
    #diamPositions = []
    
    #rotate the image by some amount and get the diameter at each point
    #the function rotates via degrees
    #we go from 0 to 89 (the 90 isn't included bc that would be a repeat of 0) 
    for phi in np.arange(0,90,1):
        rotatedRing = np.copy(I)
        rotatedRing = ndimage.rotate(rotatedRing,phi,reshape=False)
        #in order to minimize computing time I rotate only from 0 to 90 but get both the vertical and horizontal diameter to ensure that every point is still considered
        vertical = rotatedRing[:,int(npix/2)]
        horizontal = rotatedRing[int(npix/2),:]

        vertspikeOne = np.argmax(vertical)
        horspikeOne = np.argmax(horizontal)

        vertical[vertspikeOne-15:vertspikeOne+15]=0
        horizontal[horspikeOne-15:horspikeOne+15]=0

        vertspikeTwo = np.argmax(vertical)
        horspikeTwo = np.argmax(horizontal)
        
        
        diamy = np.abs(xs[vertspikeTwo]-xs[vertspikeOne])
        diamx = np.abs(xs[horspikeTwo]-xs[horspikeOne])
        if(diamx<20 or diamy<20):
          print("Tiny diameter!")
          print(diamx)
          print(diamy)
        elif(diamx>55 or diamy>55):
          print("Huge diameter!")
          print(diamx)
          print(diamy)
        else:
          diameter = diameter + (diamx+diamy)/2
          #if want to include diameter positions uncomment below
          #diamPositions.append([xs[horspikeOne],xs[horspikeTwo])
          #diamPositions will be a (2,181) array

          #if only want to plot the diamPositions for 0 rotation then uncomment lines below
          #if(phi==0):
          #   diamPositions = [xs[horspikeOne],xs[horspikeTwo]]
    diameter = diameter/90

    #ALTERNATIVE VERSION OF GETTING THE DIAMETER

    #Isumx=I.sum(axis=0)
    #spikeOnex = np.argmax(Isumx)
    #Isum[spikeOnex-20:spikeOnex+20]=0
    #spikeTwox = np.argmax(Isumx)
    #xdiameter = xs[spikeTwox]-xs[spikeOnex]

    #Isumy=I.sum(axis=1)
    #spikeOney = np.argmax(Isumy)
    #Isumy[spikeOney-20:spikeOney+20]=0
    #spikeTwoy = np.argmax(Isumy)
    #ydiameter = xs[spikeTwoy]-xs[spikeOney]
    #print(diameter)

    #Get width
    #the arctan2 just handles negative phi's so that's why i use it over regular arctan 
    phi = np.arctan2(Ys,Xs)
    #rounding phi to first decimal place bc the floats make it so no exact phi will match, but this will give us a small sliver of the ring all with the same phi
    phi=np.round(phi,1)
    unique = np.unique(phi)
    widthDistr =[]
    
    for ph in unique:
      Islice = np.copy(I)
      Islice[phi!=ph]=0
        
      spikeOne = np.unravel_index(np.argmax(Islice),Islice.shape)
        
      difference = np.abs(Islice-np.max(Islice)/2)

      leftHalf = np.copy(difference)
        #for r bigger than that of Imax put diff as = 10 so it's effectively not considered & we just look at inside portion
      leftHalf[np.sqrt(Xs**2+Ys**2)>=np.sqrt(Xs[spikeOne]**2+Ys[spikeOne]**2)] =10
       
      rightHalf = np.copy(difference)
        #for r less than that of Imax put diff as = 10 so it's effectively not considered & we just look at outside portion
      rightHalf[np.sqrt(Xs**2+Ys**2)<=np.sqrt(Xs[spikeOne]**2+Ys[spikeOne]**2)]=10
        

      half_max_firstPoint = np.unravel_index(np.argmin(leftHalf),leftHalf.shape)
      half_max_secondPoint = np.unravel_index(np.argmin(rightHalf),rightHalf.shape)

      wtemp=np.sqrt(Xs[half_max_secondPoint]**2 + Ys[half_max_secondPoint]**2)-np.sqrt(Xs[half_max_firstPoint]**2 + Ys[half_max_firstPoint]**2)

        
      edges =[0,npix-1]
        #if the half max point is on the edge of the array/image at all we say it's wrong bc the pic is big enough that the black hole should never be on the edges
      if (np.any(np.isin(half_max_firstPoint,edges)) or np.any(np.isin(half_max_secondPoint,edges))):
        print("Bad half_max point @ phi= ")
        print(ph)
      else:
        widthDistr.append([ph,wtemp])
    width=float(np.sum(widthDistr,axis=0)[1]/unique.shape)

    #Get total intensity
    totI = I.sum()

    #Append diameter, spike positions, width, widthDistribution, & total intensity
    #allRuns = np.concatenate((allRuns,np.asarray([[diameter,xs[spikeOne], xs[spikeTwo], width, widthDistr, totI]])),axis=0)
    #print(np.asarray([[diameter,xs[spikeOne], xs[spikeTwo], width, widthDistr, totI]]))

    #CHANGE THE THEORY NAME BASED ON WHICH ONE YOU ARE LOOKING AT
    mydict = [{'theory':'GR','diameter':diameter,'width':width}]
    fields=['theory','diameter','width']

    #NOTE: this appends the photon ring analysis to whichever csv file you have written
    with open("./allRuns.csv", 'a') as csvfile:
      # creating a csv writer object
      #csvwriter = csv.writer(csvfile)
      writer = csv.DictWriter(csvfile, fieldnames=fields)
 
      # writing headers (field names)
      #writer.writeheader()
 
      # writing the data rows
      writer.writerows(mydict)
      #csvwriter.writerows([diameter,width])
      csvfile.close()
    
