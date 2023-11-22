# changes .ski file

import xml.etree.ElementTree as ET
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filePath") # path to .ski file
parser.add_argument("--inc") # inclination angle (SKIRT parameter)
parser.add_argument("--az") # azimuth angle (SKIRT parameter)
parser.add_argument("--size") # length scale of region (used for field of views and spatial grid)
args = parser.parse_args()

tree = ET.parse(args.filePath)
root = tree.getroot()

FoV = float(args.size) * 2 / 3 # smaller field of view to account for rotations (factor 1/sqrt(2) smaller) 

d = {
	'FullInstrument/inclination' : str(args.inc)+'_deg',
	'FullInstrument/azimuth' : str(args.az)+'_deg',
    'FullInstrument/fieldOfViewX' : str(FoV)+'_pc',
    'FullInstrument/fieldOfViewY' :	str(FoV)+'_pc',
}

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

