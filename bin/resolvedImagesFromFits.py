import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import datetime
import matplotlib.patches as mpatches
import fitsio

def directoryStructure(tauClear, dustFraction):
	particlePath = resultPath+'Particles/'
	SKIRTPath = resultPath+'bestParamsRT/'
	noDustSKIRTPath = resultPath+'bestParamsRT/'
	plotPath = resultPath+'bestParamsPlots/'
	if eval(args.ageSmooth):
	    particlePath += 'ageSmooth/'
	    SKIRTPath += 'ageSmooth/'
	    noDustSKIRTPath += 'ageSmooth/'
	    plotPath += 'ageSmooth/'
	else:
	    particlePath += 'noAgeSmooth/'
	    SKIRTPath += 'noAgeSmooth/'
	    noDustSKIRTPath += 'noAgeSmooth/'
	    plotPath += 'noAgeSmooth/'
	if eval(args.SF):
	    particlePath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	    SKIRTPath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	    noDustSKIRTPath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	    plotPath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	else:
	    particlePath += 'noSF/'
	    SKIRTPath += 'noSF/'
	    noDustSKIRTPath += 'noSF/'
	    plotPath += 'noSF/'
	noDustPlotPath = plotPath
	SKIRTPath += 'dust/dustFraction'+str(dustFraction)+'/maxTemp'+args.maxTemp+'/'
	plotPath += 'dust/dustFraction'+str(dustFraction)+'/maxTemp'+args.maxTemp+'/'
	noDustSKIRTPath += 'noDust/'
	noDustPlotPath += 'noDust/'
	SKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	plotPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	noDustPlotPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	return SKIRTPath, plotPath, noDustSKIRTPath, noDustPlotPath, particlePath

def makeImages(galaxy, z):
    os.system('mkdir -p '+plotPath+'Images/composite/z'+z+'/'+galaxy+'/')
    for i in range(len(flux[:,0,0,0])): # loop over orientations
        #for j in range(numOrientations):
        # RGB = z, r, g:
        r_grid = flux[i,6,:,:] # sdss_z band (see band_names for indicies)
        g_grid = flux[i,4,:,:] # sdss_r band
        b_grid = flux[i,3,:,:] # sdss_g band
        fuv_grid = flux[i,0,:,:]
        w4_grid = flux[i,13,:,:]
        # set brightest pixels to value of 50th brightest pixel
        numVisable = len(r_grid[r_grid > np.amax(r_grid)*1e-5])
        numBrightest = int(numVisable*0.001)
        print('numBrightest:',numBrightest)
        max_r = r_grid.flatten()[np.argsort(r_grid.flatten())][-numBrightest]
        max_g = g_grid.flatten()[np.argsort(g_grid.flatten())][-numBrightest]
        max_b = b_grid.flatten()[np.argsort(b_grid.flatten())][-numBrightest]
        max_fuv = fuv_grid.flatten()[np.argsort(fuv_grid.flatten())][-numBrightest]
        max_w4 = w4_grid.flatten()[np.argsort(w4_grid.flatten())][-numBrightest]
        r_grid[r_grid > max_r] = max_r
        g_grid[g_grid > max_g] = max_g
        b_grid[b_grid > max_b] = max_b
        fuv_grid[fuv_grid > max_fuv] = max_fuv
        w4_grid[w4_grid > max_w4] = max_w4
        tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)])
        stretch_image = 0.005 # increase to make dim pixels dimmer
        stretch_image2 = 0.001
        image_r = np.arcsinh((r_grid/np.amax(r_grid))/stretch_image) 
        image_g = np.arcsinh((g_grid/np.amax(g_grid))/stretch_image) 
        image_b = np.arcsinh((b_grid/np.amax(b_grid))/stretch_image) 
        image_fuv = np.arcsinh((fuv_grid/np.amax(fuv_grid))/stretch_image2) 
        image_w4 = np.arcsinh((w4_grid/np.amax(w4_grid))/stretch_image2) 
        fig = plt.figure()
        sizes = np.shape(r_grid)
        fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward = False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        red = (image_r*0.8+image_w4*0.2)/np.amax(image_r*0.8+image_w4*0.2)
        green = image_g/np.amax(image_g)*0.8
        blue = ((image_b*0.90+image_fuv*0.1)/np.amax(image_b*0.90+image_fuv*0.1))*0.8
        image = np.transpose(np.asarray([red,green,blue]))
        ax.imshow(image, interpolation='none')
        plt.savefig(plotPath+'Images/composite/z'+z+'/'+galaxy+'/axisRatio'+str(axisRatios[i])+'.png', dpi=sizes[0])
        plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
#parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles parameter)
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
#parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
parser.add_argument("--SSP") # simple stellar population model including IMF after underscore (SKIRT parameter)
args = parser.parse_args()

codePath = '/home/ntf229/nihao2/'
resultPath = '/scratch/ntf229/nihao2/' # store results here
#fitsPath = '/scratch/ntf229/nihao2/bestParamsFits/ageSmooth/SF/'
fitsPath = '/scratch/ntf229/nihao2/bestParamsFits/'

if eval(args.ageSmooth):
    fitsPath += 'ageSmooth/'
else:
    fitsPath += 'noAgeSmooth/'
if eval(args.SF):
    fitsPath += 'SF/'
else:
    fitsPath += 'noSF/'

# Best parameters
tauClear = 2.5
dustFraction = 0.1

SKIRTPath, plotPath, noDustSKIRTPath, noDustPlotPath, particlePath = directoryStructure(tauClear, dustFraction)

# SKIRT broadband photometry names 
band_names = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', '2MASS_H', '2MASS_KS', 'W1', 
                'W2', 'W3', 'W4', 'PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500']

galaxies = fitsio.read(fitsPath+'nihao-integrated-seds.fits', ext='GALAXIES')

names = galaxies['name']
redshifts = list(galaxies['redshift'])

print('reshifts:', redshifts)

strRedshifts = []
for i in range(len(redshifts)):
    print('redshifts[i]:', redshifts[i])
    if redshifts[i] == 2.:
        strRedshifts.append('2.0')
    else:
        strRedshifts.append('3.6')

# flux shape: (10, 20, 500, 500)

for i in range(len(names)):
    summary = fitsio.read(fitsPath+'resolved/z'+strRedshifts[i]+'/'+names[i]+'_nihao-resolved-photometry.fits', ext='SUMMARY')
    bands = summary['bands'][0]
    axisRatios = summary['axis_ratio']
    flux = summary['flux']
    makeImages(names[i], strRedshifts[i])

print('done')
