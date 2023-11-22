# selects and stores orientations for each galaxy that span the range of axis ratios

import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--z") # redshift
parser.add_argument("--galaxy") # name of galaxy
args = parser.parse_args()

z = args.z
galaxy = args.galaxy

codePath = '/home/ntf229/nihao2/'
resultPath = '/scratch/ntf229/nihao2/' # store results here

# Directory structure stores important parameters
sampledPath = resultPath+'sampledAxisRatios/z'+z+'/'+galaxy+'/'
selectedPath = resultPath+'selectedOrientations/z'+z+'/'+galaxy+'/'

sampledIncAzAR = np.load(sampledPath+'sampledIncAzAR.npy')

selectedIncAzAr = np.zeros((10, 3))
os.system('mkdir -p '+selectedPath) 
incs = sampledIncAzAR[:,0]
azs = sampledIncAzAR[:,1]
axisRatios = sampledIncAzAR[:,2]
axisRatioBins = np.linspace(np.amin(axisRatios), 1, num=10)
for j in range(len(axisRatioBins)):
    # find axis ratios closest to bins
    index = np.abs(axisRatios - axisRatioBins[j]).argmin()
    selectedIncAzAr[j,0] = incs[index]
    selectedIncAzAr[j,1] = azs[index]
    selectedIncAzAr[j,2] = axisRatios[index]
    # remove selections to avoid duplicates
    incs = np.delete(incs, index)
    azs = np.delete(azs, index) 
    axisRatios = np.delete(axisRatios, index) 
    # copy ellipse images
    incStr = str(selectedIncAzAr[j,0])
    azStr = str(selectedIncAzAr[j,1])
    if incStr == '0.0':
        incStr = '0'
    if azStr == '0.0':
        azStr = '0'
    axisRatioStr = format(selectedIncAzAr[j,2], '.4f')
    os.system('cp '+sampledPath+'axisRatio'+axisRatioStr+'_inc'+incStr+'_az'+
                    azStr+'_ellipse.png '+selectedPath)
# save selections as numpy array
np.save(selectedPath+'selectedIncAzAR.npy', selectedIncAzAr)


