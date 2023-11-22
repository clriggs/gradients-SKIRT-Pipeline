import argparse
import os
from os.path import expanduser
from timeit import default_timer as timer
import numpy as np
import sys
import os
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from matplotlib.patches import Ellipse
import datetime

def getData(galaxy, fileName):
	inc = fileName.split('inc')[1].split('_az')[0]
	az = fileName.split('az')[1].split('_total')[0]
	images_file = fits.open(SKIRTPath+fileName)
	images = np.asarray(images_file[0].data) # cube of broadbands in MJy/sr
	r = images[0, :, :] # r band
	return inc, az, r

def plotEllipse(galaxy, inc, az, r, center, phi, ba):
	os.system('mkdir -p '+savePath)
	fig = plt.figure()
	sizes = np.shape(r)
	fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward = False)
	ax = plt.Axes(fig, [0., 0., 1., 1.])
	#ax.set_axis_off()
	fig.add_axes(ax)
	ax.imshow(r,interpolation='none', origin='lower')
	ax.add_patch(Ellipse((center[1], center[0]), width=100, height=100*ba, angle=90 - 1*phi,
		edgecolor='black',
		facecolor='none',
		linewidth=0.2))
	axisRatioStr = format(ba, '.4f')
	plt.savefig(savePath+'axisRatio'+axisRatioStr+'_inc'+inc+'_az'+az+'_ellipse.png', dpi=sizes[0])
	plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("--num") # number of orientations 
parser.add_argument("--z") # redshift 
parser.add_argument("--galaxy") # name of galaxy
args = parser.parse_args()

num = int(args.num)-1 # first one is inc=0 az=0, already included in original ski file
z = args.z
galaxy = args.galaxy

origDir = os.getcwd()
codePath=expanduser('~')+'/nihao2/'
resultPath = '/scratch/ntf229/nihao2/' # store results here

# Directory structure stores important parameters
SKIRTPath = resultPath+'sampleOrientations_SKIRT/z'+z+'/'+galaxy+'/'
savePath = resultPath+'sampledAxisRatios/z'+z+'/'+galaxy+'/'

start = timer()

numCenters = 50 # get average position of numCenters brightest pixels
numOrientations = int(args.num)

stretch = 0.01 # decreasing stretch makes dim pixels brighter
stretch_image = 0.001

minCut = 0.01

incAzAR = np.zeros((numOrientations, 3)) # inc, az, axis ratio

length = 250 # number of image pixels
positions = np.linspace(0, length-1, num=length)
xGrid = np.zeros((length, length), dtype=int)
yGrid = np.zeros((length, length), dtype=int)
for j in range(length):
	for k in range(length):
		xGrid[j,k] = int(positions[j])
		yGrid[j,k] = int(positions[k])

fileNames = np.asarray(os.listdir(SKIRTPath))
fitsMask = np.char.find(fileNames, '.fits') != -1
fileNames = fileNames[fitsMask] # only includes fits files
inc = []
az = []
for i in range(numOrientations):
	inc, az, r = getData(galaxy, fileNames[i])
	temp_r = r.copy()
	centers = np.zeros((numCenters, 2), dtype=int)
	for j in range(numCenters):
		centers[j] = np.unravel_index(np.argmax(temp_r), temp_r.shape)
		maxCut = temp_r[centers[j,0], centers[j,1]] # last loop will set maximum allowed value
		temp_r[centers[j,0], centers[j,1]] = 0 # so we don't get the same index over and over
	center = [int(np.mean(centers[:,0])), int(np.mean(centers[:,1]))]
	norm_r = r/np.amax(r)
	maxCut /= np.amax(r)
	norm_r[norm_r > maxCut] = maxCut
	norm_r = norm_r/np.amax(norm_r) # re-normalize
	norm_r[norm_r<minCut] = 0
	scaled_r = np.arcsinh(norm_r/stretch)
	image_r = np.arcsinh((r/np.amax(r))/stretch_image) # no cuts
	radius = np.zeros((length, length))
	# shift grid to be centered at center
	xGridShift = xGrid - center[0]
	yGridShift = yGrid - center[1]
	radius = np.sqrt(xGridShift**2 + yGridShift**2)
	centerMask = radius > 7
	Mxx = np.sum(scaled_r[centerMask] * xGridShift[centerMask]**2 / radius[centerMask]**2) / np.sum(scaled_r[centerMask])
	Mxy = np.sum(scaled_r[centerMask] * xGridShift[centerMask] * yGridShift[centerMask] / radius[centerMask]**2) / np.sum(scaled_r[centerMask])
	Myy = np.sum(scaled_r[centerMask] * yGridShift[centerMask]**2 / radius[centerMask]**2) / np.sum(scaled_r[centerMask])
	q = 2 * Mxx - 1
	u = 2 * Mxy
	phi = 0.5 * np.arctan2(u, q) * 180. / np.pi
	yy = q**2 + u**2
	ba = (1. + yy - 2. * np.sqrt(yy)) / (1 - yy)
	plotEllipse(galaxy, inc, az, image_r, center, phi, ba)
	incAzAR[i, 0] = float(inc)
	incAzAR[i, 1] = float(az)
	incAzAR[i, 2] = float(ba)

np.save(savePath+'sampledIncAzAR.npy', incAzAR)

end = timer()
time_SKIRT = end - start
time_SKIRT = str(datetime.timedelta(seconds=time_SKIRT))
print('Time to finish:', time_SKIRT)  

