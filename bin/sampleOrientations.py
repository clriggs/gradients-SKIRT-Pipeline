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
parser.add_argument("--num") # number of orientations 
parser.add_argument("--sim") # sim name (ex: cptmarvel, h148, r431, etc.)
parser.add_argument("--sim_dict_path") # path to pickle file where sim info is stored
parser.add_argument("--halo")
args = parser.parse_args()

num = int(args.num)-1 # first one is inc=0 az=0, already included in original ski file

#load in .pickle file storing simulation information
sim_dict=pickle.load(open(args.sim_dict_path, 'rb'))

origDir = os.getcwd()
codePath='/data/riggs/NIHAO-SKIRT-Pipeline-main/'
resultPath = '/data/riggs/SKIRT/' # store results here

# Directory structure stores important parameters
particlePath = resultPath+sim_dict[args.sim]['class']+'/'+args.sim+'/'+str(args.halo)+'/Particles/noAgeSmooth/noSF/' #the particle path is where the particle info created by "makeParticles.py" is stored
SKIRTPath = resultPath+sim_dict[args.sim]['class']+'/'+args.sim+'/'+str(args.halo)+'/sampleOrientations_SKIRT/' #where the results from this script will go

particlePath += 'z'+args.z+'/'
SKIRTPath += 'z'+args.z+'/'
particlePath += args.galaxy+'/'
SKIRTPath += args.galaxy+'/'

start = timer()

# sample orientations uniformly on the sphere 
a = np.random.uniform(low=-1.0, high=1.0, size=num)
inc = np.arccos(a) * 180 / np.pi
az = np.random.uniform(low=0.0, high=360.0, size=num)

# round inc and az to 2 decimals
inc = np.round_(inc, decimals = 2)
az = np.round_(az, decimals = 2)

# calculate size of galaxy image from text files
stars = np.load(particlePath+'stars.npy') 
xLengthStars = (np.amax(stars[:,0]) - np.amin(stars[:,0]))
yLengthStars = (np.amax(stars[:,1]) - np.amin(stars[:,1]))
zLengthStars = (np.amax(stars[:,2]) - np.amin(stars[:,2]))
maxLength = np.amax([xLengthStars, yLengthStars, zLengthStars])

os.system('mkdir -p '+SKIRTPath)
# save stars and gas text files in SKIRT directory
np.savetxt(SKIRTPath+'stars.txt', stars)
# move ski file to SKIRT directory
os.system('cp '+codePath+'resources/sampleOrientations_sph_template.ski '+SKIRTPath+'sph.ski')
# create instruments with sampled inc and az 
tree = ET.parse(SKIRTPath+'sph.ski')
root = tree.getroot()
for child in root.iter('FullInstrument'): 
    fullBB = deepcopy(child)
for child in root.iter('instruments'):
    for i in range(num):
        fullBB.set('inclination', str(inc[i])+' deg')
        fullBB.set('azimuth', str(az[i])+' deg')
        fullBB.set('instrumentName', 'inc'+str(inc[i])+'_az'+str(az[i]))
        child.insert(0,deepcopy(fullBB))
tree.write(SKIRTPath+'sph.ski', encoding='UTF-8', xml_declaration=True)
# change parameter  values in newly created .ski file 
os.system('python '+codePath+'python/sampleOrientations_modify_ski.py --filePath='+
          SKIRTPath+'sph.ski --size='+str(maxLength))
# go to SKIRT directory and run, then cd back
os.chdir(SKIRTPath)
os.system('skirt sph.ski')        
# delete radiation text files
os.system('rm stars.txt')

end = timer()
time_SKIRT = end - start
time_SKIRT = str(datetime.timedelta(seconds=time_SKIRT))
print('Time to finish:', time_SKIRT)  

