function GenerateDRRsUnetMultiVol(VolDir,OutputDir)
%% GenerateDRRsUnetMultiVol(BaseFile)
% ------------------------------------------
% FILE   : GenerateDRRsUnetMultiVol.m
% AUTHOR : Mark Gardner, The University of Sydney
% DATE   : 2022-09-12  Created.
% ------------------------------------------
% PURPOSE
%   Convert mulitple mha volumes into images that can be used as 
%   training/testing data for deep-learning networks.
% ------------------------------------------
% INPUT
%   VolumeDir:  Directory where the generated volumes are located. 
%   OutputDir:  Directory where the generated images and supporting data
%               files will be written to. 

%% Check inputs
    if nargin < 1
        VolDir = uigetdir(matlabroot,'Select the directory where the input .mha volumes are located');
    end

    if isempty(VolDir) || ~exist(VolDir,'dir')
        error(['Directory does not exist: ',VolDir])
    end

    if nargin < 2
        OutputDir = fullfile(VolDir,'TrainingData');
    end

%% Define initial variables
    close('all')

    VolumeNum = 15;     %Number of volumes created, 11 training plus 4 test volumes
    
    angles = 0:0.1:359.9;       %Projection angle of the DRRs created
       
    imageSize = [550 550 1];    %Size of images being created
    
    var_gauss = 0.00001;   %Variance of the gaussian noise added.

    ProjOffset = 0;     %X-offset of detector. Change to simulate half-fan geometry 
    
    ReconDir = fullfile(OutputDir,'Recon');
    
    if ~exist(ReconDir,'dir')
        mkdir(ReconDir)
    end    
    
    DRRDir = fullfile(OutputDir,'DRRs');
    
    if ~exist(DRRDir,'dir')
        mkdir(DRRDir)
    end
    
    %% Create geometry files
    geoFile = fullfile(DRRDir,'GeometryGAN1.xml');
    
    if ~exist(geoFile,'file')
        geoFile1 = fullfile(DRRDir,'GeometryGAN1.xml');
        system(['igtsimulatedgeometry -f 0 --arc 180 -n ',num2str(numel(angles)/2),' -o "',geoFile1,'"',...
            ' --proj_iso_x ',num2str(ProjOffset)])  
        geoFile2 = fullfile(DRRDir,'GeometryGAN2.xml');
        system(['igtsimulatedgeometry -f 180 --arc 180 -n ',num2str(numel(angles)/2),' -o "',geoFile2,'"',...
            ' --proj_iso_x ',num2str(ProjOffset)])      
    end
    
    %% Shift volumes so centre of tumour is centre pixel
    CTFile = fullfile(VolDir,'CTVol01.mha');
    
    CTShifted = fullfile(ReconDir,'CTShifted.mha'); 
     
    if ~exist(CTShifted,'file')
    
        MaskVolFile = fullfile(VolDir,'GTVVol01.mha');
        [MaskVolHeader,MaskVol] = MhaRead(MaskVolFile);

        [r,c,v] = ind2sub(size(MaskVol),find(MaskVol > 0));

        DiffR = (median(r) - MaskVolHeader.Dimensions(1)/2);%*MaskVolHeader.PixelDimensions(1);
        DiffC = (median(c) - MaskVolHeader.Dimensions(2)/2);%*MaskVolHeader.PixelDimensions(2);
        DiffV = (median(v) - MaskVolHeader.Dimensions(3)/2);%*MaskVolHeader.PixelDimensions(3);

        [CTHeader,CTVol] = MhaRead(CTFile);

        CTVol2 = NonCircularShift(CTVol,DiffR,DiffC,DiffV);
        
        MhaWrite(CTHeader,CTVol2,CTShifted)
    
    end
    
    MaskVolFile = fullfile(ReconDir,'GTVShifted.mha');
    if ~exist(MaskVolFile,'file')
        [GTVHeader,GTVVol] = MhaRead(fullfile(VolDir,'GTVVol01.mha'));
        GTVVol2 = NonCircularShift(GTVVol,DiffR,DiffC,DiffV);
        MhaWrite(GTVHeader,GTVVol2,MaskVolFile)
    end    
    
    %% Create DRRs and convert DRRs to pictures for first volume

    if ~exist(fullfile(DRRDir,'Vol01DRR1.mha'),'file')
        DRRFile = fullfile(DRRDir,'Vol01DRR.mha');
        DoForwardProjection(CTShifted,DRRFile,DRRDir)
        Projs2Pics(DRRFile,fullfile(OutputDir,'DRRImg01'))
    end        
    
    if ~exist(fullfile(DRRDir,'GTVVol01DRR1.mha'),'file')
        
        TumourDRRFile = fullfile(DRRDir,'GTVVol01DRR.mha');
        
        [CTHeader,CTVol] = MhaRead(CTShifted);
        %MaskVolFile = fullfile(BaseFile,'ContourMasks','GTVShifted.mha');
        MaskVolFile = fullfile(ReconDir,'GTVShifted.mha');
        [~,MaskVol] = MhaRead(MaskVolFile);
        
        CTVol2 = CTVol;
        CTVol2(~MaskVol) = 0;
        
        MhaWrite(CTHeader,CTVol2,fullfile(ReconDir,'GTVVol01.mha'));
        
        DoForwardProjection(fullfile(ReconDir,'GTVVol01.mha'),TumourDRRFile,DRRDir)
        
        Projs2Pics(TumourDRRFile,fullfile(OutputDir,'MaskImg01'))
        
    end
    
    %% Repeat for all other volumes

     for j = 2:VolumeNum

        CTFileLoop = fullfile(VolDir,num2str(j,'CTVol%02d.mha'));
        
        GTVFileLoop = fullfile(VolDir,num2str(j,'GTVVol%02d.mha'));
        
        CTShiftedLoop = fullfile(VolDir,num2str(j,'CTVol%02dShifted.mha'));
        
        GTVShiftedLoop = fullfile(VolDir,num2str(j,'GTVVol%02dShifted.mha'));
        
        if ~exist(CTShiftedLoop,'file')
        
            [CTHeader,CTVol] = MhaRead(CTFileLoop);
        
            CTVol2 = NonCircularShift(CTVol,DiffR,DiffC,DiffV);
            
            if j >= VolumeNum-1
                CTVol2 = imnoise(CTVol2,'gaussian',0,var_gauss);
            end
            MhaWrite(CTHeader,CTVol2,CTShiftedLoop)
            
        end
        
        if ~exist(GTVShiftedLoop,'file')
        
            [CTHeader,CTVol] = MhaRead(GTVFileLoop);
        
            CTVol2 = NonCircularShift(CTVol,DiffR,DiffC,DiffV);
            
            MhaWrite(CTHeader,CTVol2,GTVShiftedLoop)
            
        end        
        
        if ~exist(fullfile(DRRDir,num2str(j,'Vol%02dDRR1.mha')),'file')
            DRRFile = fullfile(DRRDir,num2str(j,'Vol%02dDRR.mha'));
            DoForwardProjection(CTShiftedLoop,DRRFile,DRRDir)   
            Projs2Pics(DRRFile,fullfile(OutputDir,num2str(j,'DRRImg%02d')))
        end
        
        if ~exist(fullfile(DRRDir,num2str(j,'GTVVol%02dDRR1.mha')),'file')

            TumourDRRFile = fullfile(DRRDir,num2str(j,'GTVVol%02dDRR.mha'));

            [CTHeader,CTVol] = MhaRead(CTShifted);
            [~,MaskVol] = MhaRead(GTVShiftedLoop);

            CTVol2 = CTVol;
            CTVol2(~MaskVol) = 0;

            MhaWrite(CTHeader,CTVol2,fullfile(ReconDir,num2str(j,'GTVVol%02d.mha')));

            DoForwardProjection(fullfile(ReconDir,num2str(j,'GTVVol%02d.mha')),TumourDRRFile,DRRDir)
            
            Projs2Pics(TumourDRRFile,fullfile(OutputDir,num2str(j,'MaskImg%02d')))
        end
        
    end
    
    %% Convert jp2 ims to final im
    
    FinalDir = fullfile(OutputDir,'DataSet','train');
    
    if ~exist(FinalDir,'dir')
       mkdir(FinalDir) 
    end
    
    filter_N = 3;
    a = 1; 
    b = ones(filter_N,1)./filter_N;     
    
    for j = 1:VolumeNum
        
        imdsDirectory = fullfile(OutputDir,num2str(j,'DRRImg%02d'));
        pxdsDirectory = fullfile(OutputDir,num2str(j,'MaskImg%02d'));
        
        [Images, ~] = Read_jp2(imdsDirectory);
        [Masks, ~] = Read_jp2(pxdsDirectory);
        
        if ~HalfFan       
            Images = Images(768/2 - imageSize(1)/2 +1: 768/2 + imageSize(1)/2 ,...
                1024/2 - imageSize(2)/2 +1: 1024/2 + imageSize(2)/2,:);
            Masks = Masks(768/2 - imageSize(1)/2 +1 : 768/2 + imageSize(1)/2,...
                1024/2 - imageSize(2)/2 +1 : 1024/2 + imageSize(2)/2 ,:);   
        else
            Images = Images(1:imageSize(1) ,...
                768/2 - imageSize(2)/2 +1: 768/2 + imageSize(2)/2,:);
            Masks = Masks(1:imageSize(1),...
                768/2 - imageSize(2)/2 +1 : 768/2 + imageSize(2)/2 ,:); 
        end
        
        for i = 1:numel(Images(1,1,:))
            Images(:,:,i) = filter(b, a, Images(:,:,i));
        end        
        
        for im = 1:size(Images,3)
            I = Images(:,:,im);
            M = Masks(:,:,im);
            
            I = double(I);
            I = (I - min(I(:)))/(max(I(:)) - min(I(:)));

            M = double(M);
            M = (M - min(M(:)))/(max(M(:)) - min(M(:)));

            I = im2uint16(I);
            M = im2uint16(M);            
           
            joined = [I M];
            
            fullFileName = fullfile(FinalDir,num2str([im,j],'CombineMaskDRR%04d_%02d.png'));
            imwrite(joined, fullFileName, 'BitDepth', 16);             
            
        end
        
    end
    
    %Automatically assign the data to training and testing
    AssignTrainingData(fullfile(OutputDir,'DataSet'))
    
end

function DoForwardProjection(InputVol,OutputFile,OutputDir)
% Setup the forward projections using rtkforwardprojections.exe. 
% Split into two seperate functions to save memory on the computer.
    geoFile1 = fullfile(OutputDir,'GeometryGAN1.xml');
    geoFile2 = fullfile(OutputDir,'GeometryGAN2.xml');

    [~,name,ext] = fileparts(OutputFile);
    
    OutputTemp1 = fullfile(OutputDir,[name,'1',ext]);
    OutputTemp2 = fullfile(OutputDir,[name,'2',ext]);
     
    system(['rtkforwardprojections -g "',geoFile1,'" -i "',InputVol,'" -o "', OutputTemp1,'"',...
        ' -f CudaRayCast --spacing 0.388 --dimension 768,1024'])
    system(['rtkforwardprojections -g "',geoFile2,'" -i "',InputVol,'" -o "', OutputTemp2,'"',...
        ' -f CudaRayCast --spacing 0.388 --dimension 768,1024'])    

end

function[Vol2] = NonCircularShift(Vol,DiffC,DiffR,DiffV)
% Shifts the volumes
    DiffCInt = round(DiffC);
    DiffRInt = round(DiffR);
    DiffVInt = round(DiffV);

    Vol2 = zeros(size(Vol),'like',Vol);
    
    if DiffRInt < 0
        RInds2 = abs(DiffRInt)+1:size(Vol,1);
        RInds = 1:size(Vol,2)-abs(DiffRInt);
    elseif DiffRInt > 0
        RInds2 = 1:size(Vol,2)-DiffRInt;
        RInds = DiffRInt+1:size(Vol,2);
    else
        RInds2 = 1:size(Vol,2);
        RInds = 1:size(Vol,2);
    end
    if DiffCInt < 0
        CInds2 = abs(DiffCInt)+1:size(Vol,1);
        CInds = 1:size(Vol,1)-abs(DiffCInt);
    elseif DiffCInt > 0
        CInds2 = 1:size(Vol,1)-DiffCInt;
        CInds = DiffCInt+1:size(Vol,1);
    else
        CInds2 = 1:size(Vol,1);
        CInds = 1:size(Vol,1);
    end
    if DiffVInt < 0
        VInds2 = abs(DiffVInt)+1:size(Vol,3);
        VInds = 1:size(Vol,3)-abs(DiffVInt);
    elseif DiffVInt > 0
        VInds2 = 1:size(Vol,3)-DiffVInt;
        VInds = DiffVInt+1:size(Vol,3);
    else
        VInds2 = 1:size(Vol,1);
        VInds = 1:size(Vol,1);
    end         
    for i = 1:numel(RInds)
        Vol2(CInds2,RInds2(i),VInds2) = Vol(CInds,RInds(i),VInds);
    end
    
end