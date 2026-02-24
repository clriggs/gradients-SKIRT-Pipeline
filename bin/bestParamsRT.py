# runs full RT from selected orientations which span each galaxie's range of axis ratios
# need to first run sampleOrientations.py and axisRatioDist.py to create map from 
# inc and az angles to axis ratios for each galaxy 
# this script runs SKIRT TWICE!!!!! one runs a version with dust and the other runs a version without dust!!!! 
import argparse
import os
import shutil
from os.path import expanduser
from timeit import default_timer as timer
import subprocess
import datetime
import numpy as np
import sys
import xml.etree.ElementTree as ET
from copy import deepcopy
import pickle

parser = argparse.ArgumentParser()
# parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles parameter)

parser.add_argument("--sim") # sim name (ex: cptmarvel, h148, r431, etc.) 
parser.add_argument("--sim_dict_path") # path to pickle file where simulation info is stored. EX: /resources/marvel_dcjl_sim_dict.pickle
parser.add_argument("--halo") # halo number
parser.add_argument("--distance") #distance to system in Mpc 
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp")
parser.add_argument("--FoV") #field of view of instrument in pc (SKIRT parameter)

args = parser.parse_args()

sim_dict=pickle.load(open(args.sim_dict_path, 'rb'))

origDir = os.getcwd()
codePath='/data/riggs/gradients-SKIRT-Pipeline/'
resultPath = '/data/riggs/SKIRT/'+sim_dict[args.sim]['class']+'/'+args.sim+'/'+str(args.halo)+'/' # store results here
selectedPath = resultPath+'selectedOrientations/'

# Directory structure stores important parameters
particlePath = resultPath+'Particles/'
SKIRTPath = resultPath+'bestParamsRT/'
noDustSKIRTPath = resultPath+'bestParamsRT/'

if eval(args.SF):
    particlePath += 'SF/tauClear'+args.tauClear+'/'
    SKIRTPath += 'SF/tauClear'+args.tauClear+'/'
    noDustSKIRTPath += 'SF/tauClear'+args.tauClear+'/'
    selectedPath += 'SF/tauClear'+args.tauClear+'/'
else:
    particlePath += 'noSF/'
    SKIRTPath += 'noSF/'
    noDustSKIRTPath += 'noSF/'
    selectedPath += 'noSF/'
    
SKIRTPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
noDustSKIRTPath += 'noDust/'

SKIRTPath += 'numPhotons'+args.numPhotons+'/'
noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/'

start = timer()
# calculate size of galaxy image from particles

stars = np.load(particlePath+'stars.npy') 
gas = np.load(particlePath+'gas.npy')
xLengthStars = (np.amax(stars[:,0]) - np.amin(stars[:,0]))
yLengthStars = (np.amax(stars[:,1]) - np.amin(stars[:,1]))
zLengthStars = (np.amax(stars[:,2]) - np.amin(stars[:,2]))
if len(gas) != 0:
    xLengthGas = (np.amax(gas[:,0]) - np.amin(gas[:,0]))
    yLengthGas = (np.amax(gas[:,1]) - np.amin(gas[:,1]))
    zLengthGas = (np.amax(gas[:,2]) - np.amin(gas[:,2]))
    
else:
    xLengthGas = 0
    yLengthGas = 0
    zLengthGas = 0

maxLength = np.amax([xLengthStars, yLengthStars, zLengthStars, xLengthGas, yLengthGas, zLengthGas])

print("xlength stars: ",str(xLengthStars))
print("ylength stars: ",str(yLengthStars))
print("zlength stars: ",str(zLengthStars))

print("xlength gas: ",str(xLengthGas))
print("ylength gas: ",str(yLengthGas))
print("zlength gas: ",str(zLengthGas))
selections = np.load(selectedPath+'selectedIncAzAR.npy')

selectedInc = selections[:,0]
selectedAz = selections[:,1]
selectedAxisRatio = selections[:,2]

instName = 'axisRatio'+str(np.round(selectedAxisRatio[0], decimals = 4))
# instDist = sim_dict[args.sim]['dist'][args.halo]

# including dust
if os.path.isfile(SKIRTPath+'sph_SED_'+instName+'_sed.dat'):
    print('skipping dust run')
else:
    os.system('mkdir -p '+SKIRTPath) #-p allows nested directories to be created if they don't already exist. 
    # save stars and gas text files in SKIRT directory
    np.savetxt(SKIRTPath+'stars.txt', stars)
    np.savetxt(SKIRTPath+'gas.txt', gas)
    
    if eval(args.SF):
        np.savetxt(SKIRTPath+'youngStars.txt', np.load(particlePath+'youngStars.npy'))
    else:
        os.system('touch '+SKIRTPath+'youngStars.txt') # create empty text file
    # move ski file to SKIRT directory
    os.system('cp '+codePath+'resources/customWLG.ski '+SKIRTPath+'sph.ski')
    os.system('cp '+codePath+'resources/combinedWavelengths.txt '+SKIRTPath+'combinedWavelengths.txt')
    print('copied .ski file to the SKIRTPath location')
    # change ski file values including first instrument inc and az values
    os.system('python '+codePath+'python/modify_ski.py --filePath='+SKIRTPath+
            'sph.ski --inc='+str(selectedInc[0])+' --az='+str(selectedAz[0])+
            ' --BBinstrument=broadband_'+instName+' --SEDinstrument=SED_'+instName+' --numPhotons='+args.numPhotons+' --pixels='+args.pixels+' --size='+str(maxLength)+' --dustFraction='+args.dustFraction+' --maxTemp='+args.maxTemp+' --distance='+args.distance+' --FoVX='+str(maxLength)+' --FoVY='+str(maxLength))
    # create new instruments with remaining inc and az values
    tree = ET.parse(SKIRTPath+'sph.ski')
    root = tree.getroot()
    for child in root.iter('FullInstrument'): 
        fullBB = deepcopy(child)
    for child in root.iter('SEDInstrument'): 
        fullSED = deepcopy(child)
    for child in root.iter('instruments'):
        for i in range(len(selectedInc)-1):
            instName = 'axisRatio'+str(np.round(selectedAxisRatio[i+1], decimals = 4)) 
            fullBB.set('inclination', str(selectedInc[i+1])+' deg')
            fullBB.set('azimuth', str(selectedAz[i+1])+' deg')
            fullBB.set('instrumentName', 'broadband_'+instName)
            child.insert(0,deepcopy(fullBB))
            fullSED.set('inclination', str(selectedInc[i+1])+' deg')
            fullSED.set('azimuth', str(selectedAz[i+1])+' deg')
            fullSED.set('instrumentName', 'SED_'+instName)
            child.insert(0,deepcopy(fullSED))
    tree.write(SKIRTPath+'sph.ski', encoding='UTF-8', xml_declaration=True)
    # go to SKIRT directory and run, then cd back
    os.chdir(SKIRTPath)
    os.system('skirt sph.ski') #runs SKIRT!!!
    # delete text files
    os.system('rm stars.txt')
    os.system('rm gas.txt')
    os.system('rm youngStars.txt')
    os.system('rm combinedWavelengths.txt')

instName = 'axisRatio'+str(np.round(selectedAxisRatio[0], decimals = 4))

# no dust 
if os.path.isfile(noDustSKIRTPath+'sph_SED_'+instName+'_sed.dat'):
    print('skipping no dust run')
else:
    os.system('mkdir -p '+noDustSKIRTPath)
    os.system('touch '+noDustSKIRTPath+'gas.txt') # create empty text file -- THIS IS WHAT MAKES THE RUN HAVE NO DUST!!!!!!!
    # copy stars text files to SKIRT directory
    np.savetxt(noDustSKIRTPath+'stars.txt', stars)
    if eval(args.SF):
        np.savetxt(noDustSKIRTPath+'youngStars.txt', np.load(particlePath+'youngStars.npy'))
    else:
        os.system('touch '+noDustSKIRTPath+'youngStars.txt') # create empty text file
    # move ski file to SKIRT directory
    os.system('cp '+codePath+'resources/customWLG.ski '+noDustSKIRTPath+'sph.ski')
    os.system('cp '+codePath+'resources/combinedWavelengths.txt '+noDustSKIRTPath+'combinedWavelengths.txt')
    # change ski file values including first instrument inc and az values
    os.system('python '+codePath+'python/modify_ski.py --filePath='+noDustSKIRTPath+
            'sph.ski --inc='+str(selectedInc[0])+' --az='+str(selectedAz[0])+
            ' --BBinstrument=broadband_'+instName+' --SEDinstrument=SED_'+instName+' --numPhotons='+args.numPhotons+' --pixels='+args.pixels+' --size='+str(args.FoV)+' --dustFraction='+args.dustFraction+' --maxTemp='+args.maxTemp+' --distance='+args.distance+' --FoVX='+str(args.FoV)+' --FoVY='+str(args.FoV))
    # create new instruments with remaining inc and az values
    tree = ET.parse(noDustSKIRTPath+'sph.ski')
    root = tree.getroot()
    for child in root.iter('FullInstrument'): 
        fullBB = deepcopy(child)
    for child in root.iter('SEDInstrument'): 
        fullSED = deepcopy(child)
    for child in root.iter('instruments'):
        for i in range(len(selectedInc)-1):
            instName = 'axisRatio'+str(np.round(selectedAxisRatio[i+1], decimals = 4))
            fullBB.set('inclination', str(selectedInc[i+1])+' deg')
            fullBB.set('azimuth', str(selectedAz[i+1])+' deg')
            fullBB.set('instrumentName', 'broadband_'+instName)
            child.insert(0,deepcopy(fullBB))
            fullSED.set('inclination', str(selectedInc[i+1])+' deg')
            fullSED.set('azimuth', str(selectedAz[i+1])+' deg')
            fullSED.set('instrumentName', 'SED_'+instName)
            child.insert(0,deepcopy(fullSED))
    tree.write(noDustSKIRTPath+'sph.ski', encoding='UTF-8', xml_declaration=True)
    # go to SKIRT directory and run, then cd back
    os.chdir(noDustSKIRTPath) #this cds back
    os.system('skirt sph.ski')
    # delete text files
    os.system('rm stars.txt')
    os.system('rm gas.txt')
    os.system('rm youngStars.txt')
    os.system('rm combinedWavelengths.txt')

end = timer()
time_SKIRT = end - start
time_SKIRT = str(datetime.timedelta(seconds=time_SKIRT))
print('Time to finish:', time_SKIRT)  

