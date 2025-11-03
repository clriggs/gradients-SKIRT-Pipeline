import numpy as np
import os
import argparse
from astropy.io import fits
from os.path import expanduser
import fitsio
import pickle

def getStellarMass(galaxy):
    if os.path.isfile(massPath+'/stellarMass.npy'):
        SKIRT_stellarMass = float(np.load(massPath+'/stellarMass.npy'))
    else:
        #stars = np.load(particlePath+'/stars.npy')
        #youngStars = np.load(particlePath+'/youngStars.npy')
        #SKIRT_stellarMass = np.sum(stars[:,7]) + (np.sum(youngStars[:,7] * 1.e7)) # in Msun
        stars = np.load(particlePathNoSF+'/stars.npy')
        SKIRT_stellarMass = np.sum(stars[:,7]) # in Msun
        os.system('mkdir -p '+massPath+'/')
        np.save(massPath+'/stellarMass.npy', SKIRT_stellarMass)
    return SKIRT_stellarMass

def getSFR(galaxy):
    if os.path.isfile(SFRPath+'/SFR.npy'):
        SKIRT_SFR = float(np.load(SFRPath+'/SFR.npy'))
    else:
        #stars = np.load(particlePath+'/stars.npy')
        #youngStars = np.load(particlePath+'/youngStars.npy')
        #youngMask = stars[:,9] < 1.e8 # younger than 100 Myrs
        #SKIRT_SFR = (np.sum(stars[youngMask,7]) / 1.e8) + np.sum(youngStars[:,7]) # in Msun per year
        stars = np.load(particlePathNoSF+'/stars.npy')
        youngMask = stars[:,9] < 1.e8 # younger than 100 Myrs
        SKIRT_SFR = (np.sum(stars[youngMask,7]) / 1.e8) # in Msun per year
        os.system('mkdir -p '+SFRPath+'/')
        np.save(SFRPath+'/SFR.npy', SKIRT_SFR)
    return SKIRT_SFR

def getSFH(galaxy):
    # star formation history (1D histogram)
    # ages on log scale in years
    stars = np.load(particlePathNoSF+'/stars.npy')
    ages = stars[:,9] # in years 
    masses = stars[:,7] # in M_sun
    binGrid = np.logspace(np.log10(np.amin(ages)), np.log10(np.amax(ages)), num=numBins+1)
    counts, bins = np.histogram(ages, bins=binGrid, weights=masses, density=False)
    ages = bins[:-1]
    SFH = counts
    return SFH, ages

def getCEH(galaxy):
    # chemical evolution history (2D histogram)
    # ages and metallicities on log scale
    stars = np.load(particlePathNoSF+'/stars.npy')
    ages = stars[:,9] # in years 
    masses = stars[:,7] # in M_sun
    metals = np.float64(stars[:,8])
    metals[metals==0] = np.amin(metals[metals>0]) # set to smallest nonzero value so we can log scale
    xbinGrid = np.logspace(np.log10(np.amin(ages)), np.log10(np.amax(ages)), num=numBins+1)
    ybinGrid = np.logspace(np.log10(np.amin(metals)), np.log10(np.amax(metals)), num=numBins+1)
    counts, xbins, ybins = np.histogram2d(ages, metals, bins=[xbinGrid,ybinGrid], weights=masses, density=False)
    #CEH = np.zeros((3, numBins, numBins))
    #xMesh, yMesh = np.meshgrid(xbins[:-1], ybins[:-1])
    #CEH[0] = xMesh
    #CEH[1] = yMesh
    #CEH[2] = counts
    CEH = counts
    metals = ybins[:-1]
    return CEH, metals

def getDustMass(galaxy):
    if os.path.isfile(massPath+'/metalMass.npy'):
        metalMass = float(np.load(massPath+'/metalMass.npy'))
    else:
        gas = np.load(particlePath+'/gas.npy')
        tempMask = gas[:,6] < float(16000.)
        ghostMask = np.asarray(gas[tempMask][:,4] > 0, dtype=bool) # mask out negative mass ghost particles
        metalMass = np.sum(gas[tempMask, 4][ghostMask] * gas[tempMask, 5][ghostMask]) # in Msun
        os.system('mkdir -p '+massPath+'/')
        np.save(massPath+'/metalMass.npy', metalMass)
    SKIRT_dustMass = metalMass * dustFraction
    return SKIRT_dustMass

def getAv(galaxy, instName):
    att_mask = (wave >= 912) & (wave <= 2e4)
    dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
    noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
    attenuation = dustMags[att_mask] - noDustMags[att_mask]
    Av_index = np.abs(wave[att_mask] - 5500).argmin() # find wave index closest to 5500 angstroms (V)
    Av = attenuation[Av_index]
    return Av

def reduceImageSize(image): # (x,y,band)
    newImage = np.zeros((len(image[:,0,0]), int(len(image[0,:,0])/4), int(len(image[0,0,:])/4)))
    for i in range(len(newImage[:,0,0])): # loop over bands
        for j in range(4):
            for k in range(4): 
                newImage[i,:,:] += image[i,j::4,k::4] # stride 4
        newImage[i,:,:] /= 16
    return newImage

def attenuationCurves():
    att_mask = (wave >= 912) & (wave <= 2e4)
    attenuationWave = wave[att_mask] # angstroms
    dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
    noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
    attenuationMags = dustMags[att_mask] - noDustMags[att_mask]
    return attenuationWave, attenuationMags

def energyBalance():
    c = 2.998e18 # speed of light in Anstroms per second
    freq = c / wave # in Hz
    att_mask = wave <= 2e4
    emit_mask = wave > 2e4
    attenuation = noDustSpec[att_mask] - spec[att_mask] # in Janskys
    emission = spec[emit_mask] - noDustSpec[emit_mask]
    attEnergy = -1*np.trapezoid(attenuation, freq[att_mask]) # 10^(−23) erg * s^(−1) * cm^(−2)⋅
    emitEnergy = -1*np.trapezoid(emission, freq[emit_mask]) # 10^(−23) erg * s^(−1) * cm^(−2)⋅
    return attEnergy, emitEnergy

def getSize(galaxy):
    # calculate size of galaxy image from text files
    stars = np.load(particlePath+'/stars.npy')
    gas = np.load(particlePath+'/gas.npy')
    xLengthStars = (np.amax(stars[:,0]) - np.amin(stars[:,0]))
    yLengthStars = (np.amax(stars[:,1]) - np.amin(stars[:,1]))
    zLengthStars = (np.amax(stars[:,2]) - np.amin(stars[:,2]))
    if len(gas) == 0:
        xLengthGas = 0.
        yLengthGas = 0.
        zLengthGas = 0.
    else:
        xLengthGas = (np.amax(gas[:,0]) - np.amin(gas[:,0]))
        yLengthGas = (np.amax(gas[:,1]) - np.amin(gas[:,1]))
        zLengthGas = (np.amax(gas[:,2]) - np.amin(gas[:,2]))
    maxLength = np.amax([xLengthStars, yLengthStars, zLengthStars, xLengthGas, yLengthGas, zLengthGas])
    return maxLength

parser = argparse.ArgumentParser()
# parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles $
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
parser.add_argument("--sim") # sim name (ex: cptmarvel, h148, r431, etc.) 
parser.add_argument("--sim_dict_path") # path to pickle file where simulation info is stored. EX: /resources/marvel_dcjl_sim_dict.pickle
parser.add_argument("--halo") # halo number
args = parser.parse_args()
sim_dict=pickle.load(open(args.sim_dict_path, 'rb'))

storeImages = True

numPixels = int(args.pixels)
reducedPixels = numPixels#int(numPixels/4)

dustFraction = float(0.1)

origDir = os.getcwd()
codePath='/data/riggs/gradients-SKIRT-Pipeline/'
resultPath = '/data/riggs/SKIRT/'+sim_dict[args.sim]['class']+'/'+args.sim+'/'+str(args.halo)+'/'  # store results here

selectedPath = resultPath+'selectedOrientations/'
massPath = resultPath+'GlobalProps/stellarMasses/'
SFRPath = resultPath+'GlobalProps/SFR/'
SFHPath = resultPath+'GlobalProps/SFH/'
CEHPath = resultPath+'GlobalProps/CEH/'
savePath = resultPath+'bestParamsFits/'

# Directory structure stores important parameters
SKIRTPath = resultPath+'bestParamsRT/'
particlePath = resultPath+'Particles/'
particlePathNoSF = resultPath+'Particles/' # for SFHs

if eval(args.SF):
    SKIRTPath += 'SF/tauClear'+args.tauClear+'/'
    particlePath += 'SF/tauClear'+args.tauClear+'/'
    savePath += 'SF/'
else:
    SKIRTPath += 'noSF/'
    particlePath += 'noSF/'
    savePath += 'noSF/'
particlePathNoSF += 'noSF/'

noDustSKIRTPath = SKIRTPath+'noDust/'
SKIRTPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
# savePath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'

SKIRTPath += 'numPhotons'+args.numPhotons+'/'
noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/'
# savePath += 'numPhotons'+args.numPhotons+'/'

os.system('mkdir -p '+savePath)
if storeImages:
    os.system('mkdir -p '+savePath+'resolved/')

band_names = np.asarray(['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', 
              '2MASS_H', '2MASS_KS', 'W1', 'W2', 'W3', 'W4', 'PACS70', 
              'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500',
              'SPITZER_IRAC_I1', 'SPITZER_IRAC_I2', 'SPITZER_IRAC_I3',
              'SPITZER_IRAC_I4'])
# band_names = np.asarray(['SPITZER_IRAC_I1', 'SPITZER_IRAC_I2', 'SPITZER_IRAC_I3',
              # 'SPITZER_IRAC_I4'])

# make names 2d: [[name1, z1], [name2, z2], ... ]
# names[:,0] gives all names, names[:,1] gives all redshifts
names = np.asarray([[args.sim+str(args.halo), 'z0.0']])

numOrientations = 10
numBins = 200 # number of bins for stored SFH / CEH

galaxies_dtype = [('name', str, 20),
                  ('redshift', np.float32),
                  ('stellar_mass', np.float32),
                  ('sfr', np.float32),
                  ('sfh', np.float32, numBins),
                  ('ceh', np.float32, (numBins, numBins)),
                  ('ages', np.float32, numBins),
                  ('metals', np.float32, numBins),
                  ('dust_mass', np.float32),
                  ('axis_ratios', np.float32, numOrientations),
                  ('size', np.float32)]

summary_dtype = [('name', str, 20),
                 ('redshift', np.float32),
                 ('stellar_mass', np.float32),
                 ('sfr', np.float32),
                 ('dust_mass', np.float32),
                 ('axis_ratio', np.float32),
                 ('Av', np.float32),
                 ('attenuated_energy', np.float32),
                 ('emitted_energy', np.float32),
                 ('bands', band_names.dtype, len(band_names)),
                 ('flux', np.float32, len(band_names)),
                 ('flux_nodust', np.float32, len(band_names)),
                 ('size', np.float32)]

image_summary_dtype = [('name', str, 20),
                 ('redshift', np.float32),
                 ('stellar_mass', np.float32),
                 ('sfr', np.float32),
                 ('dust_mass', np.float32),
                 ('axis_ratio', np.float32),
                 ('Av', np.float32),
                 ('bands', band_names.dtype, len(band_names)),
                 ('flux', np.float32, (len(band_names), reducedPixels, reducedPixels)),
                 ('flux_nodust', np.float32, (len(band_names), reducedPixels, reducedPixels)),
                 ('size', np.float32)]

galaxies = np.zeros(len(names), dtype=galaxies_dtype)
summary = np.zeros(len(names) * numOrientations, dtype=summary_dtype)

wave = None
attenuation_mags = None

indx = 0

for i in range(len(names[:,0])):
    if storeImages:
        image_galaxies = np.zeros(1, dtype=galaxies_dtype)
        image_summary = np.zeros(numOrientations, dtype=image_summary_dtype)
    stellarMass = getStellarMass(names[i,:])
    SFR = getSFR(names[i,:])
    SFH, ages = getSFH(names[i,:])
    CEH, metals = getCEH(names[i,:])
    dustMass = getDustMass(names[i,:])
    selections = np.load(selectedPath+'/selectedIncAzAR.npy') # [inc, az, axisRatio]
    axisRatios = selections[:,2]
    size = getSize(names[i,:])
    galaxies['name'][i] = names[i,0]
    galaxies['redshift'][i] = float(names[i,1].split('z')[1])
    galaxies['stellar_mass'][i] = stellarMass
    galaxies['sfr'][i] = SFR
    galaxies['sfh'][i] = SFH
    galaxies['ceh'][i] = CEH
    galaxies['ages'][i] = ages
    galaxies['metals'][i] = metals
    galaxies['dust_mass'][i] = dustMass
    galaxies['axis_ratios'][i] = axisRatios
    galaxies['size'][i] = size
    if storeImages:
        image_galaxies['name'] = names[i,0]
        image_galaxies['redshift'] = float(names[i,1].split('z')[1])
        image_galaxies['stellar_mass'] = stellarMass
        image_galaxies['sfr'] = SFR
        image_galaxies['sfh'] = SFH
        image_galaxies['ceh'] = CEH
        image_galaxies['ages'] = ages
        image_galaxies['metals'] = metals
        image_galaxies['dust_mass'] = dustMass
        image_galaxies['axis_ratios'] = axisRatios
        image_galaxies['size'] = size
    for j in range(numOrientations-1): # loop through orientations
        instName = 'axisRatio' + str(np.round(axisRatios[j+1], decimals=4)) # orientation naming system
        # Spatially integrated
        # print('SKIRTPath=',SKIRTPath)
        # print('noDustSKIRTPath=',noDustSKIRTPath)
        BB = np.loadtxt(SKIRTPath+'sph_broadband_'+instName+'_sed.dat', unpack = True)[1] # spatially integrated broadband fluxes in Janskys
        
        noDustBB = np.loadtxt(noDustSKIRTPath+'sph_broadband_'+instName+'_sed.dat', unpack = True)[1]  # spatially integrated broadband fluxes in Janskys
        sed = np.loadtxt(SKIRTPath+'sph_SED_'+instName+'_sed.dat', unpack = True)
        if (wave is None): # only need to do once
            wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
            spectrum = np.zeros((len(names) * numOrientations, len(wave)), dtype=np.float32)
            spectrum_nodust = np.zeros((len(names) * numOrientations, len(wave)), dtype=np.float32)
        spec = sed[1] # spatially integrated SED fluxes in Janskys
        noDustSed = np.loadtxt(noDustSKIRTPath+'sph_SED_'+instName+'_sed.dat', unpack = True)
        noDustSpec = noDustSed[1]
        Av = getAv(names[i], instName)
        attenuationWave, attenuationMags = attenuationCurves()
        if (attenuation_mags is None):
            attenuation_mags = np.zeros((len(names) * numOrientations, len(attenuationWave)), dtype=np.float32)
        attEnergy, emitEnergy = energyBalance()
        # Store integrated data
        summary['name'][indx] = galaxies['name'][i]
        summary['redshift'][indx] = galaxies['redshift'][i]
        summary['stellar_mass'][indx] = galaxies['stellar_mass'][i]
        summary['sfr'][indx] = galaxies['sfr'][i]
        summary['dust_mass'][indx] = galaxies['dust_mass'][i]
        summary['axis_ratio'][indx] = galaxies['axis_ratios'][i, j]
        summary['Av'][indx] = Av
        summary['bands'][indx] = band_names
        summary['flux'][indx, :] = BB
        summary['flux_nodust'][indx, :] = noDustBB
        summary['attenuated_energy'][indx] = attEnergy
        summary['emitted_energy'][indx] = emitEnergy
        summary['size'][indx] = galaxies['size'][i]
        spectrum[indx, :] = spec
        spectrum_nodust[indx, :] = noDustSpec
        attenuation_mags[indx, :] = attenuationMags
        # Spatially resolved
        if storeImages:
            imageFile = fits.open(SKIRTPath+'sph_broadband_'+instName+'_total.fits')
            imageBB = np.asarray(imageFile[0].data) # (20, 2000, 2000) bands in MJy/st
            # imageBB = reduceImageSize(imageBB) # (20, 500, 500)
            noDustImageFile = fits.open(noDustSKIRTPath+'sph_broadband_'+instName+'_total.fits')
            noDustImageBB = np.asarray(noDustImageFile[0].data) # (20, 2000, 2000) bands in MJy/st
            # noDustImageBB = reduceImageSize(noDustImageBB) # (20, 500, 500)
            # Store resolved data
            image_summary['name'][j] = galaxies['name'][i]
            image_summary['redshift'][j] = galaxies['redshift'][i]
            image_summary['stellar_mass'][j] = galaxies['stellar_mass'][i]
            image_summary['sfr'][j] = galaxies['sfr'][i]
            image_summary['dust_mass'][j] = galaxies['dust_mass'][i]
            image_summary['axis_ratio'][j] = galaxies['axis_ratios'][i, j]
            image_summary['Av'][j] = Av
            image_summary['bands'][j] = band_names
            image_summary['flux'][j, :] = imageBB
            image_summary['flux_nodust'][j, :] = noDustImageBB
            image_summary['size'][j] = galaxies['size'][i]
        # Increase indx
        indx += 1
        if storeImages:
            os.system('mkdir -p '+savePath+'resolved/')
            fitsio.write(savePath+'resolved/'+'_resolved-photometry.fits', image_galaxies, extname='GALAXIES', clobber=True)
            fitsio.write(savePath+'resolved/'+'_resolved-photometry.fits', image_summary, extname='SUMMARY', clobber=False)

fitsio.write(savePath+'integrated-seds.fits', galaxies, extname='GALAXIES', clobber=True)
fitsio.write(savePath+'integrated-seds.fits', summary, extname='SUMMARY', clobber=False)
fitsio.write(savePath+'integrated-seds.fits', wave, extname='WAVE', clobber=False)
fitsio.write(savePath+'integrated-seds.fits', spectrum, extname='SPEC', clobber=False)
fitsio.write(savePath+'integrated-seds.fits', spectrum_nodust, extname='SPECNODUST', clobber=False)
fitsio.write(savePath+'integrated-seds.fits', attenuationWave, extname='ATTENUATION_WAVE', clobber=False)
fitsio.write(savePath+'integrated-seds.fits', attenuation_mags, extname='ATTENUATION_MAGS', clobber=False)

print('done')
