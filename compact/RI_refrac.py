import xml.etree.ElementTree  as ET
import argparse
import os
import numpy as np

# This file is checking if last RefractiveIndecis from SingleCrystal_cosmic_ray 
# are different from DualTestBeam xml files
# to run: 
#     python RI_refrac.py -f1 DRSingleCrystalCosmicRay.xml -f2 DRDualTestBeam.xml


argParser = argparse.ArgumentParser()
argParser.add_argument("-f1", "--xml_filename1", help="SingleCrystal_cosmic_ray file")
argParser.add_argument("-f2", "--xml_filename2", help="Dualtestbeam file")
args = argParser.parse_args()

print(args.xml_filename1, args.xml_filename2)

names = []
values = {}

def diff_RIfiles(xml_file):
  tree = ET.parse(xml_file)
  root = tree.getroot()
  for child in root:
    if child.tag!='properties': continue
    for node in child:
      if node.attrib['name'].startswith("RI"):
        name = node.attrib['name']
        names.append(name)
        v= node.attrib['values'].split(' ')
        v = [s for s in v if s]
        v = [float(s) for s in v if 'eV' not in s]
        values[name] = np.array(v)
  return names, values

print(os.getcwd() + '/' + args.xml_filename2)
name_file1, values_file1 = diff_RIfiles(os.getcwd() + '/' + args.xml_filename1)
name_file2, values_file2 = diff_RIfiles(os.getcwd() + '/' + args.xml_filename2)

for name in name_file1:
    diff_val = np.subtract(values_file1[name], values_file2[name])
    if diff_val.any() !=0: print(name)
