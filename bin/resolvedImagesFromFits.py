import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import datetime
import matplotlib.patches as mpatches
import fitsio
import pickle

def directoryStructure(dustFraction):
	particlePath = resultPath+'Particles/'
	SKIRTPath = resultPath+'bestParamsRT/'
	noDustSKIRTPath = resultPath+'bestParamsRT/'
	plotPath = resultPath+'bestParamsPlots/'
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
	SKIRTPath += 'numPhotons'+args.numPhotons+'/'#+args.SSP+'/'
	plotPath += 'numPhotons'+args.numPhotons+'/'#+args.SSP+'/'
	noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/'#+args.SSP+'/'
	noDustPlotPath += 'numPhotons'+args.numPhotons+'/'#+args.SSP+'/'
	return SKIRTPath, plotPath, noDustSKIRTPath, noDustPlotPath, particlePath

def makeImages(galaxy, z):
    os.system('mkdir -p '+plotPath+'Images/'+args.instrument+'/z'+str(z)+'/'+galaxy+'/')
    for i in range(len(flux[:,0,0,0])): # loop over orientations
        #for j in range(numOrientations):
        # RGB = z, r, g:
        if args.instrument == 'composite':
            r_grid = flux[i,6,:,:] # sdss_z band (see band_names for indicies)
            g_grid = flux[i,4,:,:] # sdss_r band
            b_grid = flux[i,3,:,:] # sdss_g band
            fuv_grid = flux[i,0,:,:]
            w4_grid = flux[i,13,:,:]
            numVisable = len(r_grid[r_grid > np.amax(r_grid)*1e-5])
        elif args.instrument == 'spitzer':
            i1_grid = flux[i,20,:,:]
            i2_grid = flux[i,21,:,:]
            i3_grid = flux[i,22,:,:]
            i4_grid = flux[i,23,:,:]
            numVisable = len(i1_grid[i1_grid > np.amax(i1_grid)*1e-5])
        # set brightest pixels to value of 50th brightest pixel
        
        # numVisable = len(i1_grid[i1_grid > np.amax(i1_grid)*1e-5])
        numBrightest = int(numVisable*0.001)
        print('numBrightest:',numBrightest)
        if args.instrument == 'composite':
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
        elif args.instrument == 'spitzer':
            max_i1 = i1_grid.flatten()[np.argsort(i1_grid.flatten())][-numBrightest]
            max_i2 = i2_grid.flatten()[np.argsort(i2_grid.flatten())][-numBrightest]
            max_i3 = i3_grid.flatten()[np.argsort(i3_grid.flatten())][-numBrightest]
            max_i4 = i4_grid.flatten()[np.argsort(i4_grid.flatten())][-numBrightest]
            i1_grid[i1_grid > max_i1] = max_i1
            i2_grid[i2_grid > max_i2] = max_i2
            i3_grid[i3_grid > max_i3] = max_i3
            i4_grid[i4_grid > max_i4] = max_i4
            tot_max = np.amax([np.amax(i1_grid),np.amax(i2_grid),np.amax(i3_grid)])
            stretch_image = 0.005 # increase to make dim pixels dimmer
            stretch_image2 = 0.001
            image_i1 = np.arcsinh((i1_grid/np.amax(i1_grid))/stretch_image) 
            image_i2 = np.arcsinh((i2_grid/np.amax(i2_grid))/stretch_image)
            image_i3 = np.arcsinh((i3_grid/np.amax(i3_grid))/stretch_image)
            image_i4 = np.arcsinh((i4_grid/np.amax(i4_grid))/stretch_image)
        fig = plt.figure()
        if args.instrument == 'composite':
            sizes=np.shape(r_grid)
        elif args.instrument == 'spitzer':
            sizes = np.shape(i1_grid)
        fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward = False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        if args.instrument == 'composite':
            red = (image_r*0.8+image_w4*0.2)/np.amax(image_r*0.8+image_w4*0.2)
            green = image_g/np.amax(image_g)*0.8
            blue = ((image_b*0.90+image_fuv*0.1)/np.amax(image_b*0.90+image_fuv*0.1))*0.8
            image = np.transpose(np.asarray([red,green,blue]))
        elif args.instrument == 'spitzer':
            image = image_i1
            image = np.transpose(image) #rotates the image!
        ax.imshow(image, interpolation='none')
        plt.savefig(plotPath+'Images/'+args.instrument+'/z'+str(z)+'/'+galaxy+'/axisRatio'+str(axisRatios[i])+'.png', dpi=sizes[0])
        plt.close()

parser = argparse.ArgumentParser()
# parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles parameter)
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
parser.add_argument("--sim") # sim to use
parser.add_argument("--sim_dict_path") # path to the dictionary with sim and halo info
parser.add_argument("--halo") # halo to use
parser.add_argument("--instrument")
args = parser.parse_args()

sim_dict = pickle.load(open(args.sim_dict_path, 'rb'))

codePath = '/data/riggs/gradients-SKIRT-Pipeline/'
resultPath = '/data/riggs/SKIRT/'+sim_dict[args.sim]['class']+'/'+args.sim+'/'+str(args.halo)+'/' # store results here
fitsPath = resultPath+'bestParamsFits/'

# if eval(args.ageSmooth):
#     fitsPath += 'ageSmooth/'
# else:
#     fitsPath += 'noAgeSmooth/'
if eval(args.SF):
    fitsPath += 'SF/'
else:
    fitsPath += 'noSF/'



# # Best parameters
tauClear = int(args.tauClear)
dustFraction = args.dustFraction

SKIRTPath, plotPath, noDustSKIRTPath, noDustPlotPath, particlePath = directoryStructure(dustFraction)

# SKIRT broadband photometry names 
band_names = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', '2MASS_H', '2MASS_KS', 'W1', 
                'W2', 'W3', 'W4', 'PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500',
                'SPITZER_IRAC_I1', 'SPITZER_IRAC_I2', 'SPITZER_IRAC_I3', 'SPITZER_IRAC_I4']
# band_names = np.asarray(['SPITZER_IRAC_I1', 'SPITZER_IRAC_I2', 'SPITZER_IRAC_I3', 'SPITZER_IRAC_I4'])

galaxies = fitsio.read(fitsPath+'integrated-seds.fits', ext='GALAXIES')

names = galaxies['name']
redshifts = list(galaxies['redshift'])

print('reshifts:', redshifts)

strRedshifts = []
for i in range(len(redshifts)):
    print('redshifts[i]:', redshifts[i])
    if redshifts[i] == 2.:
        strRedshifts.append('2.0')
    else:
        strRedshifts.append(redshifts[i])

# flux shape: (10, 20, 500, 500)

for i in range(len(names)):
    summary = fitsio.read(fitsPath+'resolved/_resolved-photometry.fits', ext='SUMMARY')
    bands = summary['bands'][0]
    axisRatios = summary['axis_ratio']
    flux = summary['flux']
    makeImages(names[i], strRedshifts[i])

print('done')
