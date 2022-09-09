function GenerateMultipleVolumes(BaseFile,InputFile,DeformationParamFile)
%% GenerateMultipleVolumes(BaseFile)
% ------------------------------------------
% FILE   : GenerateMultipleVolumes.m
% AUTHOR : Mark Gardner, The University of Sydney
% DATE   : 2022-09-09  Created.
% ------------------------------------------
% PURPOSE
%   Deforms the base planning CT volumes.
% ------------------------------------------
% INPUT
%   BaseFile:           The path to the directory where the HNSCC files are
%                       stored
%   InputjsonFile:      Path to template Json file that contains the info for creating
%                       the rotations.
%   elastixParamFile:   Path to elastix parameter file

    if nargin < 1
        BaseFile = uigetdir(matlabroot,'Select directory where the HNSCC files are stored');    
    end
    if nargin < 2
        InputFile = fullfile(BaseFile,'RotationInfo.json');
    end
    if nargin < 3
        DeformationParamFile = fullfile(BaseFile,'Elastix_BSpline_OpenCL_RigidPenalty.txt');
    end

    close('')
    
    %Vols = No motion(1), Tumour shift up down in AP and up down in SI
    %(4?). Head rotations (6?)
    
    VolumeNum = 15;     %Number of volumes created, 11 training plus two test volumes
    
    % Alter to adapt to data directory layout
    Files = {'HNSCC-01-0001/12-05-1998-RT SIMULATION-43582',...
    'HNSCC-01-0003/11-20-2001-RT SIMULATION-30560',...
    'HNSCC-01-0004/08-24-1996-RT SIMULATION-72882',...
    'HNSCC-01-0005/01-19-1998-NA-RT SIMULATION-39358',...
    'HNSCC-01-0007/04-29-1997-RT SIMULATION-32176',...
    'HNSCC-01-0040/12-12-1998-RT SIMULATION-35323',...
    'HNSCC-01-0039/12-14-1998-RT SIMULATION-26503',...
    'HNSCC-01-0070/05-10-1999-RT SIMULATION-21211',...
    'HNSCC-01-0101/09-25-1999-NA-RT SIMULATION-93304',...
    'HNSCC-01-0094/09-08-1999-NA-RT SIMULATION-68058',...    
    'HNSCC-01-0002/11-20-2001-RT SIMULATION-57976',...
    'HNSCC-01-0018/03-01-2009-NA-RT SIMULATION-16942',...    
    'HNSCC-01-0122/04-18-2000-RT SIMULATION-20757',...
    'HNSCC-01-0123/03-24-2000-RT SIMULATION-47852',...    
    'HNSCC-01-0130/05-29-2000-NA-RT SIMULATION-72920'   
    };

    %Rotation components for head/neck to be written to InputFile. In
    %degrees
    yaw = -2.5;     %Positive is Head Right
    pitch = -3.1;    %Positive is Head Down
    roll = -3.2;    %Positive is Head Tilt right
    
    folderStr = 'DeepLearningOutput';
    
    %% Create generic json file with deformation information
    
    %Modify Json File
    jsonText = fileread(InputFile);
    
    % Convert JSON formatted text to MATLAB data types (3x1 cell array in this example)
    jsonData = jsondecode(jsonText);     
    
    jsonData.angles = [abs(yaw);abs(roll);abs(pitch)];
    jsonData.axes = abs(jsonData.axes);
    if yaw < 0
       jsonData.axes(1,1) = -1; 
    end
    if pitch < 0
       jsonData.axes(2,3) = -1; 
    end
    if roll < 0
       jsonData.axes(3,2) = -1; 
    end    
    
    % Convert to JSON text
    jsonText2 = jsonencode(jsonData);
    
    % Write to a json file
    fid = fopen(InputFile, 'w');
    fprintf(fid, '%s', jsonText2);
    fclose(fid);    
    
    %% loop through number of patients
    for i = 1:numel(Files)
        
        PatientDir = fullfile(BaseFile,Files{i});
        
        DeepLearningDir = fullfile(PatientDir,folderStr);
        
        if ~exist(DeepLearningDir,'dir')
            mkdir(DeepLearningDir)
        end
        
        FinalVols = fullfile(DeepLearningDir,'TrainingVols');
        
        if ~exist(FinalVols,'dir')
            mkdir(FinalVols)
        end        
        
        %% Create generic json file for each patient
        RotFile = fullfile(DeepLearningDir,'RotationInfo.json');        
        
        %Do python script to create scan specific rotation .json file
        [errVal,msg] = system(['python CreateRotationInfo.py "',InputFile,'" "',strrep(RotFile,'\','/'),'"']);
        if errVal ~= 0
            system(['python ',which('CreateRotationInfo.py'),' "',InputFile,'" "',strrep(RotFile,'\','/'),'"'])
        end        
        
        %Modify Json File
        jsonText = fileread(RotFile);

        % Convert JSON formatted text to MATLAB data types (3x1 cell array in this example)
        jsonData = jsondecode(jsonText);         
        
        %Add tumour shift to deformation
        
        if i < 6 
            StructShift = [2.0,1.7,0];    %Bruijnen et al. 2019 pg4
        elseif i >=6 && i < 11
            StructShift = [3.8,2.2,0];  %Bruijnen et al. 2019 pg4
        elseif i >= 11
            StructShift = [0,0,0];
        end
        
        if i == 13
            StructName = "PRE_CHEMO_GTV";
        elseif i == 15
            StructName = "POST_CHEMO_GTV";
        else
            StructName = "GTV";
        end
        
        jsonData.Structure_Shift = [StructShift];
        jsonData.Structure_Names = [StructName];        
        
        % Convert to JSON text
        jsonText2 = jsonencode(jsonData);

        % Write to a json file
        fid = fopen(RotFile, 'w');
        fprintf(fid, '%s', jsonText2);
        fclose(fid);         
        
        %% modify transformation file
        
        %Copy elastix parameter file to folder and modify to allow for
        %rigid/non-rigid registration
        NewParamFile = fullfile(DeepLearningDir,'Elastix_BSpline_OpenCL_RigidPenalty.txt');        
        
        if ~exist(NewParamFile,'file')
           
            copyfile(DeformationParamFile,NewParamFile)

            newstr = ['"',DeepLearningDir,'\nifti\BoneMaskInit.nii.gz','"'];

            ModifyElastixParameterFile(NewParamFile,'MovingRigidityImageName',newstr)        
        end        
        
        RotFileLoop = fullfile(DeepLearningDir,'RotationInfoLoop.json');   
        
        if ~exist(fullfile(PatientDir,'CT.mha'),'file')
            dicomList = lscell([fullfile(PatientDir,'dicom','ct'),'/*.dcm']);
            if isempty(dicomList)
                disp('Something')
            end
            DICOM2IEC(dicomList,fullfile(PatientDir,'CT.mha'))
        end
        
        [EgHeader,~] = MhaRead(fullfile(PatientDir,'CT.mha'));
        
        %% Loop through and create volumes
        
        for j = 1:VolumeNum 

            CTVolFile = fullfile(FinalVols,num2str(j,'CTVol%02d.mha'));
            
            if ~exist(CTVolFile,'file')
            
                %% Modify json file to only include one rotation/transformation
                jsonText = fileread(RotFile);
                jsonData = jsondecode(jsonText); 

                NewAngles = [0;0;0];
                NewStructShift = [0;0;0];

                if j == 2
                   NewAngles(1) = jsonData.angles(1)*2;
                end
                if j == 3
                   NewAngles(2) = jsonData.angles(2)*2;
                end
                if j == 4 
                   NewAngles(3) = jsonData.angles(3)*2;
                end
                if j == 5 %|| j == VolumeNum 
                   NewAngles(1) = jsonData.angles(1)*2;
                   jsonData.axes = abs(jsonData.axes);
                end
                if j == 6 %|| j == VolumeNum 
                   NewAngles(2) = jsonData.angles(2)*2;
                   jsonData.axes = abs(jsonData.axes);
                end
                if j == 7 %|| j == VolumeNum 
                   NewAngles(3) = jsonData.angles(3)*2;
                   jsonData.axes = abs(jsonData.axes);
                end
                if j == 8
                   NewStructShift(1) = jsonData.Structure_Shift(1)*2;
                end
                if j == 9
                   NewStructShift(2) = jsonData.Structure_Shift(2)*2;
                end
                if j == 10
                   NewStructShift(1) = -jsonData.Structure_Shift(1)*2;
                end  
                if j == 11
                   NewStructShift(2) = -jsonData.Structure_Shift(2)*2;
                end              
                if j < 12
                    jsonData.angles = NewAngles;
                    jsonData.Structure_Shift = NewStructShift;
                end
                if j == 13
                    jsonData.angles = jsonData.angles;
                    jsonData.Structure_Shift = -jsonData.Structure_Shift;
                    %jsonData.axes = abs(jsonData.axes);
                end
                if j == 14
                    % 25% deformation magnitude
                    jsonData.angles = jsonData.angles/2;
                    jsonData.Structure_Shift = jsonData.Structure_Shift/2;
                elseif j == 15
                    % 75% deformation magnitude
                    jsonData.angles = 3*(jsonData.angles/2);
                    jsonData.Structure_Shift = 3*(jsonData.Structure_Shift/2);
                end
                % Convert to JSON text
                jsonText2 = jsonencode(jsonData);

                % Write to a json file
                fid = fopen(RotFileLoop, 'w');
                fprintf(fid, '%s', jsonText2);
                fclose(fid);             

                %% Do Rotation
                [errVal,msg] = system(['python VirtualPhantom/MultipleHeadRotationsShift.py "',RotFileLoop,'"']);
                if errVal ~= 0
                    errVal = system(['python "',which('MultipleHeadRotationsShift.py'),'" "',RotFileLoop,'"']);
                    if errVal ~= 0
                        disp(['Error ',PatientDir]) 
                       continue
                    end
                end            

                data = ReadJsonData(RotFileLoop);

                [~,name,ext] = fileparts(RotFileLoop);

                %% Save created output volumes

                CTFile = fullfile(DeepLearningDir,'output_nifti',data.name,name,'ct','ct.nii.gz');
                GTVFile = fullfile(DeepLearningDir,'output_nifti',data.name,name,'rtstruct','GTV.nii.gz');
                if ~exist(GTVFile,'file')
                    GTVFile = fullfile(DeepLearningDir,'output_nifti',data.name,name,'rtstruct','POST_CHEMO_GTV.nii.gz');
                    if ~exist(GTVFile,'file')
                        GTVFile = fullfile(DeepLearningDir,'output_nifti',data.name,name,'rtstruct','PRE_CHEMO_GTV.nii.gz');
                        if ~exist(GTVFile,'file')
                           disp('Something went wrong. GTV file not found')
                           continue
                        end
                    end
                end

                %dicomFileList = lscell([fullfile(DeepLearningDir,data.name,name,'ct'),'/*dcm']); 
                [CTVol,CTHeader] = TransformVolumes(CTFile,EgHeader);
                [GTVVol,GTVHeader] = TransformVolumes(GTVFile,EgHeader);

                MhaWrite(CTHeader,CTVol,CTVolFile);
                MhaWrite(GTVHeader,GTVVol,fullfile(FinalVols,num2str(j,'GTVVol%02d.mha')));
            end
        end
        
    end
    
end

function[Vol2,Header2] = TransformVolumes(File,Header)

    waterAtt = 0.013;

    Vol = niftiread(File);
    Info = niftiinfo(File);
    
    IndOrder = [1 3 2];
    
    Header2 = Header;
    Header2.PixelDimensions = Info.PixelDimensions(IndOrder);
    Header2.Dimensions = Info.ImageSize(IndOrder);
    Header2.Offset = -0.5 * (Info.ImageSize(IndOrder) - 1) .* Info.PixelDimensions(IndOrder);

    if range(Vol(:)) > 2
        if min(Vol(:)) >= 0
            Vol = single(Vol) * waterAtt / 1000;
        else
            Vol = (single(Vol) + 1000) * waterAtt / 1000;
        end   
    else
        Vol = single(Vol);
    end
    Vol2 = permute(Vol(:,end:-1:1,end:-1:1),IndOrder);
    
end
