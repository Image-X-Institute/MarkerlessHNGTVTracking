# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 12:02:38 2021

Create a rotation file from one patient info file and one rotation template

@author: mgar5380
"""

import json
import sys
import os

def CreateRotationInfo(InputFile,OutputFile):
    
    FileParts = OutputFile.split('/')
    #sep = os.path.sep
    sep = '/'       
    
    BaseDir = sep.join(FileParts[:-2])
    inputDirectory = BaseDir + "/dicom"
    
    #intermediateDirectory = BaseDir + "/DeepLearningNew"
    #OutputDirectory = BaseDir + "/DeepLearningNew"    
    intermediateDirectory = BaseDir + sep + FileParts[-2]
    OutputDirectory = BaseDir + sep + FileParts[-2]
    
    with open(InputFile) as json_file:
        VolTemplateData = json.load(json_file) 
    
    #Read info from CT file
    InfoFile = BaseDir + "/VolumeDeformations/CTScanInfo.json"
    
    with open(InfoFile) as json_file:
        CTScanInfo = json.load(json_file)
     
    for i in range(0,len(VolTemplateData['axes'])):
        AxisArr = VolTemplateData['axes'][i]
    
        if abs(AxisArr[0]) == 1:
            RotationPts = CTScanInfo['C1-C2']            
        elif abs(AxisArr[1]) == 1:
            RotationPts = CTScanInfo['C2-C3']
        elif abs(AxisArr[2]) == 1:
            RotationPts = CTScanInfo['Oc-C1']
    
        VolTemplateData['point_of_rotation'][i] = RotationPts
    
    VolTemplateData['coordinates_cutoff'] = CTScanInfo['coordinates_cutoff']
    
    #data['coordinates_cutoff_lateral'] = CTScanInfo['coordinates_cutoff_lateral']
    
    VolTemplateData['name'] = CTScanInfo['name']
    VolTemplateData['input_directory'] = inputDirectory
    VolTemplateData['intermediate_directory'] = intermediateDirectory
    VolTemplateData['output_directory'] = OutputDirectory

    with open(OutputFile,'w') as json_file:
        json.dump(VolTemplateData,json_file)  

if __name__ == '__main__':
    #Map command line arguments to function arguments.
    CreateRotationInfo(*sys.argv[1:])

#InputFile = "C://Users/mgar5380/Documents/Writing/MarkerlessHNTracking/RotationInfo.json"
#OutputFile = "Z://2RESEARCH/2_ProjectData/RemoveTheMask/CTData/HNSCC/HNSCC-01-0018/03-01-2009-NA-RT SIMULATION-16942/DeepLearningMultiVol/RotationInfo.json"

#CreateRotationInfo(InputFile, OutputFile)