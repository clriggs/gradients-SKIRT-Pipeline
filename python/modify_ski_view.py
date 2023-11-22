# changes .ski file
# begin line with 'SKIRT/' for skirt parameters
# use a unique parent header after the '/' to point to each parameter
# example: SKIRT/GeometricSource/scaleLength "2000 pc"

import xml.etree.ElementTree as ET
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filePath") # path to .ski file
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
parser.add_argument("--size") # length scale of region (used for field of views and spatial grid)
parser.add_argument("--dustFraction") # dust-to-metal ratio 
parser.add_argument("--maxTemp") # maximum temperature at which dust can form
parser.add_argument("--SSP") # simple stellar population model including IMF after underscore (SKIRT parameter)
args = parser.parse_args()

# Update SSP model and IMF
ssp = args.SSP.split('_')[0]
imf = args.SSP.split('_')[1]
# Read in the file
with open(args.filePath, 'r') as file:
    filedata = file.read()
# Replace the target string
filedata = filedata.replace('FSPSSEDFamily', ssp+'SEDFamily')
if ssp == 'Bpass':
    filedata = filedata.replace('imf="Chabrier"','') # BpassSEDFamily has no imf option, always uses Chabrier
else:
    filedata = filedata.replace('Chabrier', imf)
# Write the file out again
with open(args.filePath, 'w') as file:
    file.write(filedata)

tree = ET.parse(args.filePath)
root = tree.getroot()

# MonteCarloSimulation/numPackets 1e7
# FullInstrument/inclination 0_deg
# radiationFieldWLG/LogWavelengthGrid/numWavelengths
# dustEmissionWLG/LogWavelengthGrid/numWavelengths
# defaultWavelengthGrid/LogWavelengthGrid/numWavelengths

minXYZ = -float(args.size) / 2
maxXYZ = float(args.size) / 2

FoV = float(args.size) * 2 / 3 # smaller field of view to account for rotations (factor 1/sqrt(2) smaller) 

d = {
    'FullInstrument/numPixelsX' : str(args.pixels),
    'FullInstrument/numPixelsY' : str(args.pixels),
	'MonteCarloSimulation/numPackets' : str(args.numPhotons),
    'PolicyTreeSpatialGrid/minX' : str(minXYZ)+'_pc',
    'PolicyTreeSpatialGrid/maxX' : str(maxXYZ)+'_pc',
    'PolicyTreeSpatialGrid/minY' : str(minXYZ)+'_pc',
    'PolicyTreeSpatialGrid/maxY' : str(maxXYZ)+'_pc',
    'PolicyTreeSpatialGrid/minZ' : str(minXYZ)+'_pc',
    'PolicyTreeSpatialGrid/maxZ' : str(maxXYZ)+'_pc',
    'FullInstrument/fieldOfViewX' : str(FoV)+'_pc',
    'FullInstrument/fieldOfViewY' :	str(FoV)+'_pc',
    'ParticleMedium/massFraction' : str(args.dustFraction),
    'ParticleMedium/maxTemperature' : str(args.maxTemp)+'_K'
}


# store config.txt inputs as dictionary 
#d = {}
#with open(args.projectPath+"/config.txt") as f:
#	for line in f:
#		(key, val) = line.split()
#		d[key] = val
#print(d)


for name, value in d.items():
	s = name.split('/')
	print('split:', s)
	print('length:', len(s))
	print(s[-1], value)

	for item in root.iter(s[0]):
		if len(s) == 2:
			item.set(s[-1], value.replace("_", " "))
		if len(s) == 3:
			for sub_item in item:
				sub_item.set(s[-1], value.replace("_", " "))

tree.write(args.filePath, encoding='UTF-8', xml_declaration=True)

