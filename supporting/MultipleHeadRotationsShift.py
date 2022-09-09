# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 16:34:05 2021

like MultipleHeadRotations but does the GTV shift after the rigid/non-rigid registration

@author: mgar5380
"""

# Basics
from pathlib import Path
import numpy as np
import SimpleITK as sitk
import re
import os
import sys
import json
import ntpath

sys.path.append('C://Users/mgar5380/Documents/GitHub/platipy')

# Tools for conversion from 
#from platipy.dicom.dicom_directory_crawler.conversion_utils import process_dicom_directory
from platipy.dicom.io.crawl import process_dicom_directory

# Visualisation tool
#from platipy.imaging.visualisation.tools import ImageVisualiser

# Center of mass
#from platipy.imaging.utils.tools import get_com

# Deformation field operations
#from platipy.imaging.registration.registration import apply_field
from platipy.imaging.registration.utils import apply_transform
#from platipy.imaging.deformation_fields.deformation_field_operations import (get_external_mask, get_bone_mask, generate_field_radial_bend)
from platipy.imaging.generation.mask import (
    get_bone_mask,
    get_external_mask)
from platipy.imaging.generation.dvf import generate_field_shift

#Perform field rotations
from headMaskFn import get_head_mask
from tryRotate import generate_field_rotation

# Tools for conversion back to DICOM
#from platipy.dicom.nifti_to_rtstruct.convert import convert_nifti
from platipy.dicom.io.nifti_to_rtstruct import convert_nifti
#from platipy.dicom.nifti_to_series.convert import convert_nifti_to_dicom_series
from platipy.dicom.io.nifti_to_series import convert_nifti_to_dicom_series

def CombineDVFs(DVFUnSmoothedImg,DVFSmoothedImg,coordinates_cutoff,RotationPt):
    
    DVFUnSmoothed = sitk.GetArrayFromImage(DVFUnSmoothedImg)
    DVFSmoothed = sitk.GetArrayFromImage(DVFSmoothedImg)
    #ImageCTArray = sitk.GetArrayFromImage(image_ct)
    ImageMask = np.zeros(DVFSmoothed.shape)
    
    m1 = (RotationPt[-1][0] - coordinates_cutoff[0][0])/(RotationPt[-1][1] - coordinates_cutoff[0][1])
    m2 = (RotationPt[-1][0] - coordinates_cutoff[1][0])/(RotationPt[-1][1] - coordinates_cutoff[1][1])
    m3 = (coordinates_cutoff[0][0] - coordinates_cutoff[1][0])/(coordinates_cutoff[0][1] - coordinates_cutoff[1][1])
    
    b1 = coordinates_cutoff[0][0] - m1*coordinates_cutoff[0][1]
    b2 = coordinates_cutoff[1][0] - m2*coordinates_cutoff[1][1]
    b3 = coordinates_cutoff[1][0] - m3*coordinates_cutoff[1][1]    
    
    for i in range(0,ImageMask.shape[0]):
        for j in range(0,ImageMask.shape[1]):
            if i < m1*j + b1 and i < m2*j + b2: #and i > m3*j + b3:
                ImageMask[i,j,:,:] = DVFSmoothed[i,j,:,:]
            else:
                ImageMask[i,j,:,:] = DVFUnSmoothed[i,j,:,:]   

    DVFNew = sitk.GetImageFromArray(ImageMask)    
    DVFNew.CopyInformation(DVFSmoothedImg)
    return DVFNew

def RotateAroundPoint(RTform,RotationPt,RotationOrigin):
    RotationVector = np.array(RotationPt) - np.array(RotationOrigin)
    NewPts = RTform.apply(RotationVector) + RotationOrigin
    return NewPts.tolist()
    pass

def main(input_file="Z://2RESEARCH/2_ProjectData/RemoveTheMask/CTData/HNSCC/HNSCC-01-0040/12-12-1998-RT SIMULATION-35323/DeepLearningMultiVol2/RotationInfo.json",noise=True):

    #default values
    output_dcm_dir = './output_dcm'
    intermediate_dir = None
    #gauss_smooth = 10
    
    # if len(sys.argv) < 2:
    #     print("Not enough arguments. Need json file")
    #     sys.exit(0)
    # input_file = str(sys.argv[1])
    
    #input_file = "Z://2RESEARCH/2_ProjectData/RemoveTheMask/CTData/HNSCC/VolumeDeformations/Rotation1.json"
    
    # read json file
    with open(input_file) as json_file:
        data = json.load(json_file)
    
    
    """
    Sample json file:  
    {
    	"name": "HNSCC-01-0033",
    	"input_directory": "C:/Users/cche2532/Documents/HNSCC/HNSCC-01-0033/dicom", 
    	"axes":"[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]",
    	"angles":"[1,2,3,4,5]",
    	"point_of_rotation":"[66,247,259]",
    	"intermediate_directory": "C:/Users/cche2532/Documents/PlatipyProject/platipy/platipy/examples/generate_synthetic_deformations",
    	"output_directory": "C:/Users/cche2532/Documents/PlatipyProject/platipy/platipy/examples/generate_synthetic_deformations",
    	"gaussian_smooth": "0", 
    	"coordinates_cutoff": "[66,247,259],[44,198,259]"
    }
    """
    #check existence of essential json values
    try:
        patient_id = data['name']
    except Exception as e:
        print('Import patient name failed. Error was:')
        print(e)
        sys.exit(0)
    try:
        input_dcm_dir = data['input_directory']
    except Exception as e:
        print('Import input dicom directory failed. Error was:')
        print(e)
        sys.exit(0)
    
    # patient name
    patient_id = data['name']
    #input dicom directory
    
    
    #check input dcm directory exists
    if not os.path.isdir(input_dcm_dir):
        print('Error with input dicom directory (May not exist). Exiting')
        sys.exit(0)
        
    #load axes
    #axes = []  
    try:   
        axes = data['axes']
        #temp_str = data['axes'].strip().split('[')
        # for a in temp_str[1:len(temp_str)]:
        #     aint = list(map(int, a.rstrip('],').split(','))) 
        #     axes.append(aint)
        
    except Exception as e:
        print('Reading axes failed due to:')
        print(e)
        sys.exit(0)
    
    #load angles
    #angles = []  
    try:   
        angles = data['angles']
        #temp_str = data['angles'].strip().split('[')
        #for a in temp_str[1:len(temp_str)]:
        #    angles = list(map(float, a.rstrip('],').split(','))) 
    except Exception as e:
        print('Reading angles failed due to:')
        print(e)
        sys.exit(0)
        
    #load points of rotation
    #point_of_rotation = []
    try: 
        point_of_rotation = data['point_of_rotation']
        #temp_str = data['point_of_rotation'].strip().split('[')
        #for a in temp_str[1:len(temp_str)]:
            #point_of_rotation = list(map(int, a.rstrip('],').split(','))) 
        #for a in temp_str[1:len(temp_str)]:
        #    aint = list(map(int, a.rstrip('],').split(','))) 
        #    point_of_rotation.append(aint)            
    except Exception as e:
        print('Reading point of rotation failed due to:')
        print(e)
        sys.exit(0)
        
    coordinates_cutoff = []
    #coordinates_cutoff.append(point_of_rotation)
        
    #non essential json values
    try:
        output_dcm_dir = data['output_directory']
    except Exception as e:
        print('No output dicom directory provided. Using default ./output_dcm')
    try:
        intermediate_dir = data['intermediate_directory']
    except Exception as e:
        print('No intermediate dicom directory provided. Using default .')
    try:
        gauss_smooth = float(data["gaussian_smooth"])
    except Exception as e:
        print('No gaussian smoothing parameter provided. Using default of 0.')
    temp_str = None
    try: 
        coordinates_cutoff = data['coordinates_cutoff']
    except Exception as e:
        print('No coordinates provided for head mask cutoff. Using default (with point of rotation).')
    # if temp_str is not None:
    #     try: 
            
    #         #coordinates_cutoff = []
    #         #temp_str = temp_str.strip().split('[')
    #         #for a in temp_str[1:len(temp_str)]:
    #         #    aint = list(map(int, a.rstrip('],').split(','))) 
    #         #    coordinates_cutoff.append(aint)
    #     except Exception as e: 
    #         print('Reading axes failed due to:')
    #         print(e)
    #         sys.exit(0)
            
            
    if intermediate_dir is None:
        nifti_directory = Path("./nifti")
    else:
        nifti_directory = Path(intermediate_dir + "/nifti")
    
    patientunderscore = patient_id.replace('-','_')
    
    #ct_path = next(iter(nifti_directory.glob(f"{patientunderscore}/IMAGES/*CT*.nii.gz")))

    # convert dcm to nifti
    if not Path(nifti_directory).is_dir():        
        process_dicom_directory(input_dcm_dir,output_directory=nifti_directory)
    ct_path = next(iter(nifti_directory.glob(f"{patientunderscore}/IMAGES/*CT*.nii.gz")))
    
    # Read in image    
    image_ct = sitk.ReadImage(str(ct_path))
    
    # Read in structures
    structure_paths = nifti_directory.glob(f"{patientunderscore}/STRUCTURES/*.nii.gz")
    structures = {re.findall(".*_RTSTRUCT_(.*).nii.gz", p.name)[0]: sitk.ReadImage(str(p)) for p in structure_paths}
    
    # Automatically generate external mask
    struct_external = get_external_mask(image_ct)
    
    #Change to read multiple points of rotation
    #tup = tuple(point_of_rotation)
    
    #generate head mask using list of coordinates as cutoff points
    head_mask = get_head_mask(struct_external,coordinates_cutoff, image_ct)
    
    #if len(coordinates_cutoff_Lateral) > 0:
    #    head_mask = get_head_mask(struct_external,coordinates_cutoff_Lateral, image_ct)
    #    pass
    
    RotationTFormTotal = []
    
    #FinalTransform = sitk.Transform(3,sitk.sitkComposite)
    
    image_ct_Loop = []
    
    RotationPt = [];
    
    Mask_deformed = [];
    dvf_transform = [];
    dvf_field = sitk.Image(image_ct.GetSize(), sitk.sitkVectorFloat64);
    sitk.Image.CopyInformation(dvf_field,image_ct)
    
    if np.count_nonzero(angles) > 0:
    
        try: 
            len(axes[0])
            LoopNum = len(axes)
        except:
            LoopNum = 1
        
        try:
            angles[0]
        except:
            angles = [angles]
        
        for i in range(0,LoopNum):
            #Transform new axex to match frame of reference of old axis
            #if i == 0:
            try:
                len(axes[i])
                axisLoop = axes[i]
            except:
                axisLoop = axes
            
            try:
                len(point_of_rotation[i])
                RotationPt.append(point_of_rotation[i])
                
            except:
                RotationPt.append(point_of_rotation)
            # if len(point_of_rotation[i]) == 1:
            #     RotationPt.append(point_of_rotation)
            # else:
            #     RotationPt.append(point_of_rotation[i])
            #RotationPt.append(point_of_rotation[i])
            Mask_deformed_loop, dvf_transform_loop, dvf_field_loop, RotationTForm = generate_field_rotation(
                image_ct,
                head_mask,
                tuple(RotationPt[i]),
                axis_of_rotation=axisLoop,
                angle=angles[i]*np.pi/180,
                gaussian_smooth=0
                #mask_bend_from_reference_point = coordinates_cutoff
            )
            Mask_deformed.append(Mask_deformed_loop)
            dvf_transform.append(dvf_transform_loop)
            #dvf_field.append(dvf_field_loop)
            dvf_field = dvf_field + dvf_field_loop
          
            #FinalTransform = sitk.CompositeTransform([FinalTransform,dvf_transform_loop])         
    
    try:
        if np.count_nonzero(data['Structure_Shift']) > 0:
            structShiftFlag = 1
        else:
            structShiftFlag = 0   
    except:
        structShiftFlag = 0      
    
    if structShiftFlag:
        if np.count_nonzero(data['Structure_Shift']) > 0:
            if np.size(data['Structure_Shift'][0]) > 1:
                for i in range(0,len(data['Structure_Shift'])):
                    
                    Mask_shift, dvf_shift, dvf_field_shift = generate_field_shift(
                        mask_image = structures[data['Structure_Names'][i]],
                        vector_shift = tuple(data['Structure_Shift'][i]),
                        gaussian_smooth = 2
                    )
                    #dvf_transform.append(dvf_shift)
                    #dvf_New = dvf_New + sitk.Cast(dvf_field_shift, sitk.sitkVectorFloat64)
                    dvf_field = dvf_field + dvf_field_shift
                    dvf_transform.append(dvf_shift)
                    Mask_deformed.append(Mask_shift)
            else:
                Mask_shift, dvf_shift, dvf_field_shift = generate_field_shift(
                    mask_image = structures[data['Structure_Names']],
                    vector_shift = tuple(data['Structure_Shift']),
                    gaussian_smooth = 2
                )                
                dvf_transform.append(dvf_shift)
                Mask_deformed.append(Mask_shift)
                
                dvfLoop = sitk.Cast(dvf_field_shift, sitk.sitkVectorFloat32)  
                dvfLoop.CopyInformation(dvf_field)     
                
                try:
                    dvf_field = dvf_field + dvf_field_shift
                except:
                    dvf_field = dvf_field + sitk.Cast(dvfLoop,sitk.sitkVectorFloat64)   
                    # dvf_transform.append(dvf_shift)
                    # dvfLoop = sitk.Cast(dvf_field_shift, sitk.sitkVectorFloat32)  
                    # dvfLoop.CopyInformation(dvf_New)
                    # try:
                    #     dvf_New = dvf_New + dvfLoop
                    # except:
                    #     dvf_New = dvf_New + sitk.Cast(dvfLoop,sitk.sitkVectorFloat64)      
    
    #gauss_smooth = 5
    
    #dvf_field_smoothed = sitk.SmoothingRecursiveGaussian(dvf_field, gauss_smooth)
      
    #FinalTransform = sitk.DisplacementFieldTransform(
    #    sitk.Cast(dvf_field_smoothed, sitk.sitkVectorFloat64)
    #)
    
    #image_ct_deformed = apply_field(image_ct, FinalTransform, structure=False, interp=sitk.sitkLinear)

    #dvf_Final = CombineDVFs(dvf_field,dvf_field_smoothed,coordinates_cutoff,RotationPt)
    #dvf_Final = dvf_field    
    
    if np.count_nonzero(angles) > 0 or structShiftFlag:
    
        InitialTransform = sitk.DisplacementFieldTransform(
            sitk.Cast(dvf_field, sitk.sitkVectorFloat64)
        )
        
        #image_ct_deformed = apply_field(image_ct, InitialTransform, structure=False, interp=sitk.sitkLinear)
        image_ct_deformed = apply_transform(image_ct, transform=InitialTransform, interpolator=sitk.sitkLinear)
        
        #Prepare data for rigid/non-rigid registration
        #Normalise and write initial and deformed CT vols
        image_ct_norm = sitk.Normalize(image_ct)
        sitk.WriteImage(image_ct_norm,str(nifti_directory) + '/CTNorm.nii.gz')
    
        image_ct_deformed_norm = sitk.Normalize(image_ct_deformed)
        sitk.WriteImage(image_ct_deformed_norm,str(nifti_directory) + '/CTWarpedNorm.nii.gz')
    
        #Create Bone Masks
        BoneMaskInit = get_bone_mask(image_ct)
        BoneMaskWarped = get_bone_mask(image_ct_deformed)
        BoneMaskInitNorm = sitk.Normalize(BoneMaskInit)
        BoneMaskWarpedNorm = sitk.Normalize(BoneMaskWarped)
        sitk.WriteImage(BoneMaskInitNorm,str(nifti_directory) + '/BoneMaskInit.nii.gz')
        sitk.WriteImage(BoneMaskWarpedNorm,str(nifti_directory) + '/BoneMaskWarped.nii.gz')
    
        #Do registration
        FixedImage = str(nifti_directory) +  '/CTWarpedNorm.nii.gz'
        MovingImage = str(nifti_directory) +  '/CTNorm.nii.gz'
        parameterFilePath = intermediate_dir + '/Elastix_BSpline_OpenCL_RigidPenalty.txt'
        outputPath = str(nifti_directory)
    
        os.system('elastix -f "{}" -m "{}" -p "{}" -out "{}"'.format(FixedImage,MovingImage,parameterFilePath,outputPath))
    
        ParamFile = outputPath + '/TransformParameters.0.txt'
    
        #os.system('transformix -in "{}" -tp "{}" -out "{}"'.format(ct_path,ParamFile,ct_deformed_File))
    
        #Create DVF from registration results
        os.system('transformix -def all -tp "{}" -out "{}"'.format(ParamFile, outputPath))
    
        if Path(outputPath + '/deformationField.mha'):
            dvf_New = sitk.ReadImage(outputPath + '/deformationField.mha')
        elif Path(outputPath + '/deformationField.nii.gz'):
            dvf_New = sitk.ReadImage(outputPath + '/deformationField.nii.gz')
        else:
            assert False, 'Deformation File failed to be read'
    else:
        dvf_New = dvf_field
    #dvf_New = dvf_field

    #Do GTV Shift
    
    # try:
    #     data['Structure_Shift']
    #     structShiftFlag = 1
    # except:
    #     structShiftFlag = 0      
    
    # if structShiftFlag:
    #     if np.count_nonzero(data['Structure_Shift']) > 0:
    #         if np.size(data['Structure_Shift'][0]) > 1:
    #             for i in range(0,len(data['Structure_Shift'])):
                    
    #                 Mask_shift, dvf_shift, dvf_field_shift = generate_field_shift(
    #                     mask_image = structures[data['Structure_Names'][i]],
    #                     vector_shift = tuple(data['Structure_Shift'][i]),
    #                     gaussian_smooth = 2
    #                 )
    #                 #dvf_transform.append(dvf_shift)
    #                 dvf_New = dvf_New + sitk.Cast(dvf_field_shift, sitk.sitkVectorFloat64)
    #         else:
    #             Mask_shift, dvf_shift, dvf_field_shift = generate_field_shift(
    #                 mask_image = structures[data['Structure_Names']],
    #                 vector_shift = tuple(data['Structure_Shift']),
    #                 gaussian_smooth = 2
    #             )
    #             dvf_transform.append(dvf_shift)
    #             dvfLoop = sitk.Cast(dvf_field_shift, sitk.sitkVectorFloat32)  
    #             dvfLoop.CopyInformation(dvf_New)
    #             try:
    #                 dvf_New = dvf_New + dvfLoop
    #             except:
    #                 dvf_New = dvf_New + sitk.Cast(dvfLoop,sitk.sitkVectorFloat64)    

    FinalTransform = sitk.DisplacementFieldTransform(
        sitk.Cast(dvf_New, sitk.sitkVectorFloat64)
    )
    
    #Use DVF to deform CT
    #image_ct_deformed2 = apply_field(image_ct, FinalTransform, structure=False, interp=sitk.sitkLinear)
    image_ct_deformed2 = apply_transform(image_ct, transform=FinalTransform, interpolator=sitk.sitkLinear)

    image_ct_deformed_norm2 = sitk.Normalize(image_ct_deformed2)
    sitk.WriteImage(image_ct_deformed_norm2,str(nifti_directory) + '/CTWarpedNormFinal.nii.gz')

    strname = ntpath.basename(input_file)
    str2 = strname.split('.')[0]
    
    
    #apply dvf field to deformed_structures
    deformed_structures = {}
    for struct in structures:
        #deformed_structures[struct] = apply_field(structures[struct], dvf_transform, structure=True, default_value=0,interp=sitk.sitkLinear)
        #deformed_structures[struct] = apply_field(structures[struct], FinalTransform, structure=True, default_value=0,interp=sitk.sitkLinear)
        deformed_structures[struct] = apply_transform(structures[struct], transform=FinalTransform, default_value=0,interpolator=sitk.sitkLinear)
    # struct_external_def = apply_field(struct_external, dvf_transform, structure=True, default_value=0, interp=sitk.sitkLinear)

    #nifti structure files output dir
    if intermediate_dir is None:
        output_dir_STRUCT = Path("./output_nifti/"+patient_id + "/" + str2 + "/rtstruct")
    else:
        output_dir_STRUCT = Path(intermediate_dir + "/output_nifti/"+patient_id + "/" + str2 + "/rtstruct")
    #output_dir_STRUCT = Path(intermediate_dir + '/TestDir/rtstruct')
    output_dir_STRUCT.mkdir(exist_ok=True, parents=True)
    
    output_dir_CT = Path(intermediate_dir + "/output_nifti/"+patient_id + "/" + str2 + "/ct")
    output_dir_CT.mkdir(exist_ok=True, parents=True)
    
    #write CT Volume to nifti format
    sitk.WriteImage(image_ct_deformed2,str(output_dir_CT) + '/ct.nii.gz')
    
    #write structure to nifti format 
    for struct in deformed_structures:
        sitk.WriteImage(deformed_structures[struct],str(output_dir_STRUCT / str(struct + ".nii.gz")))
    
    #ct dcm files output dir
    output_dir_dcm_CT = Path(output_dcm_dir+"/"+patient_id + "/" + str2 + "/ct")
    #output_dir_dcm_CT = Path(intermediate_dir + '/TestDir/ct')
    output_dir_dcm_CT.mkdir(exist_ok=True, parents=True)
    
    # if noise:
    #     print("Adding Noise to deformed image")
    #     image_ct_deformed2 = sitk.AdditiveGaussianNoise(image_ct_deformed,
    #                                                    standardDeviation=200.0,
    #                                                    mean=0.0)
    # else:
    #     image_ct_deformed2 = image_ct_deformed
    #     #pass
    
    #write deformed ct to dcm
    convert_nifti_to_dicom_series(
        image=image_ct_deformed2,
        reference_dcm=input_dcm_dir+"/CT/",
        output_directory=str(output_dir_dcm_CT)
    )
    
    
    # rtstruct output dir 
    output_dir_dcm_STRUCT = Path(output_dcm_dir+"/"+patient_id + "/" + str2 +  "/rtstruct")
    output_dir_dcm_STRUCT.mkdir(exist_ok=True, parents=True)
    
    # dictionary containing path of nifti files
    masks = {}
    for m in os.listdir(output_dir_STRUCT):
        name = m.split('.')[0]
        mask_path = str(output_dir_STRUCT / m)
        masks[name] = mask_path

    #convert structure format from nifti to dcm 
    # convert_nifti(
    #     dcm_file = str(output_dir_dcm_CT / "0.dcm"),
    #     mask = masks, 
    #     out_rt_filename = str(output_dir_dcm_STRUCT / "struct.dcm")
    # )   
    
    #dcm_path, mask_input, output_file
    
    convert_nifti(
        dcm_path = str(output_dir_dcm_CT / "0.dcm"),
        mask_input = masks, 
        output_file = str(output_dir_dcm_STRUCT / "struct.dcm")
    )    
        
    if np.count_nonzero(angles) > 0:
        #delete temp files used in transformation
        os.remove(Path(outputPath + '/deformationField.mha'))
        os.remove(Path(outputPath + '/result.0.mha'))
        os.remove(Path(ParamFile))
    
    # for axis in axes:
    #     for ang in angles:
    #         # generate deformation field
    #         image_ct_deformed, dvf_transform, dvf_field = generate_field_rotation(
    #             image_ct,
    #             head_mask,
    #             tup,
    #             axis_of_rotation=axis,
    #             angle=ang*np.pi/180,
    #             gaussian_smooth=gauss_smooth
    #         )
            
    #         # resample image and structures with linear interpolation
    
    #         image_ct_deformed = apply_field(image_ct, dvf_transform, structure=False, interp=sitk.sitkLinear)
    #         str2 = 'axis_' + '_'.join(str(x) for x in axis) + '_angle_' + str(ang)
            
    #         #apply dvf field to deformed_structures
    #         deformed_structures = {}
    #         for struct in structures:
    #             deformed_structures[struct] = apply_field(structures[struct], dvf_transform, structure=True, default_value=0,interp=sitk.sitkLinear)
    #         # struct_external_def = apply_field(struct_external, dvf_transform, structure=True, default_value=0, interp=sitk.sitkLinear)
    
    #         #nifti structure files output dir
    #         if intermediate_dir is None:
    #             output_dir_STRUCT = Path("./output_nifti/"+patient_id + "/" + str2 + "/rtstruct")
    #         else:
    #             output_dir_STRUCT = Path(intermediate_dir + "/output_nifti/"+patient_id + "/" + str2 + "/rtstruct")
            
    #         output_dir_STRUCT.mkdir(exist_ok=True, parents=True)
            
    #         #write structure to nifti format 
    #         for struct in deformed_structures:
    #             sitk.WriteImage(deformed_structures[struct],str(output_dir_STRUCT / str(struct + ".nii.gz")))
            
    #         #ct dcm files output dir
    #         output_dir_dcm_CT = Path(output_dcm_dir+"/"+patient_id + "/" + str2 + "/ct")
    #         output_dir_dcm_CT.mkdir(exist_ok=True, parents=True)
            
    #         #write deformed ct to dcm
    #         convert_nifti_to_dicom_series(
    #             image=image_ct_deformed,
    #             reference_dcm=input_dcm_dir+"/CT/",
    #             output_directory=str(output_dir_dcm_CT)
    #         )
            
    #         # rtstruct output dir 
    #         output_dir_dcm_STRUCT = Path(output_dcm_dir+"/"+patient_id + "/" + str2 +  "/rtstruct")
    #         output_dir_dcm_STRUCT.mkdir(exist_ok=True, parents=True)
            
    #         # dictionary containing path of nifti files
    #         masks = {}
    #         for m in os.listdir(output_dir_STRUCT):
    #             name = m.split('.')[0]
    #             mask_path = str(output_dir_STRUCT / m)
    #             masks[name] = mask_path
    
    #         #convert structure format from nifti to dcm 
    #         convert_nifti(
    #             dcm_file = str(output_dir_dcm_CT / "0.dcm"),
    #             mask = masks, 
    #             out_rt_filename = str(output_dir_dcm_STRUCT / "struct.dcm")
    #         )

if __name__ == '__main__':
    # Map command line arguments to function arguments.
    main(*sys.argv[1:])