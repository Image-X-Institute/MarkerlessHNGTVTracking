function DeepLearningRandEvaluation(ResultsFile,UnWarpedFile,OutputDir,videoFlag)
%% function DeepLearningRandEvaluation
% ------------------------------------------
% FILE   : DeepLearningRandEvaluation.m
% AUTHOR : Mark Gardner, The University of Sydney
% DATE   : 2022-09-12  Created.
% ------------------------------------------
% PURPOSE
%   Analyse the markerless tracking results for each patient.
% ------------------------------------------
% INPUT
%   ResultsFile:    The directory of images which is the output of the
%                   test.py command for each patient.
%   UnWarpedFile:   Directory to the training data (.pngs files) used to train the
%                   netowrk that is being tested.
%   OutputDir:      The directory where the output images and data will be
%                   saved.
%   VideoFlag:      Whether a video of the results will be created. 


%% Input check

    if nargin < 1
        ResultsDir = uigetdir(matlabroot,'Select the directory where the results images are stored');
    end

    if isempty(ResultsDir) || exist(ResultsDir,'dir')
        error(['Directory does not exist: ',ResultsDir])
    end

    if nargin < 2
        UnWarpedFile = uigetdir(matlabroot,'Select the directory where the training images are stored');
    end

    if isempty(UnWarpedFile) || exist(UnWarpedFile,'dir')
        error(['Directory does not exist: ',UnWarpedFile])
    end

    if nargin < 3 || isempty(OutputDir)
        OutputDir = fullfile(ResultsDir,'Results');        
    end

    if ~exist(OutputDir,'dir')
        mkdir(OutputDir)
    end

    if nargin < 4
        videoFlag = 1;
    end

    close('all')
   
    UnWarpedFile = fullfile(BaseFile,'DeepLearningMultiVol2','UnWarpedIms'); 

    str = 'CombineMaskDRR%04d_%02d';        %General String for input images

    StartingVal = 1;
    
    angles = 0:0.1:359.9;
    
    AreaThresh = 1000;
    
    BinaryThreshold = 0.05;

    pixelpermm = 0.256;
    
    DICE = [];
    PredictedCentroid = [];
    RealCentroid = [];
    MSDNorm = [];     
    
    %Results for if you assume no motion
    UnWarpedDICE = [];
    UnWarpedCentroid = [];
    UnWarpedMSD = [];
    
    %% Setup video and plotting

    figure('Menu','none','ToolBar','none')
    ax1 = gca;
    fig1 = gca;
    
    FrameRate = 5;
    
    if videoFlag
        video = VideoWriter(fullfile(OutputDir,'ResultsNormModel1'));  
        video2 = VideoWriter(fullfile(OutputDir,'ResultsNormModel2')); 
        video3 = VideoWriter(fullfile(OutputDir,'ResultsNormModel3')); 
        video4 = VideoWriter(fullfile(OutputDir,'ResultsNormModel4'));

        video.FrameRate = FrameRate;
        video2.FrameRate = FrameRate;
        video3.FrameRate = FrameRate;
        video4.FrameRate = FrameRate;

        open(video);
        open(video2);
        open(video3);
        open(video4);
    end
    
    % Get list of images
    ImageList = lscell([ResultsFile,'/*real_A.png']);
    
    VolumeFlag = zeros(1,numel(ImageList));

    %% Main analysis loop
    for i = 1:numel(ImageList)
        
        % Get image name for all files

        realstr2 = ImageList{i};
        
        [~,name,ext] = fileparts(ImageList{i});
        
        ImNum = str2double(name(15:18));
        ModelNum = str2double(name(20:21));
        
        str2  = [num2str([ImNum,ModelNum],str),'_'];
       
        realstr = fullfile(ResultsFile,[str2,'real_B.png']);
        fakestr = fullfile(ResultsFile,[str2,'fake_B.png']);
       
       if ModelNum == 12
            VolumeFlag(i) = 0;
       elseif ModelNum == 13
            VolumeFlag(i) = 1;
       elseif ModelNum == 14
            VolumeFlag(i) = 2;
       elseif ModelNum == 15
            VolumeFlag(i) = 3;
       end
       
       if exist(realstr,'file')
           % Read images
           try
               PredictionIm = imread(fakestr); 
           catch
               error(['File ',fakestr,' cannot be read'])
           end
           
           try
               TruthIm = imread(realstr); 
           catch
              error(['File ',realstr,' cannot be read'])
           end
           try
               ActualIm = imread(realstr2); 
           catch
               error(['File ',realstr2,' cannot be read'])
           end
           TruthImBin = TruthIm;
           
            %Binarise input images

           if numel(find(TruthIm == 0)) < 100
               TruthImBin = TruthIm - median(TruthIm(:));
               TruthIm = TruthIm - median(TruthIm(:));
           end           
           
           TruthImBin(TruthImBin > 0) = 1;
           
           PredictionImBin = im2double(PredictionIm - median(PredictionIm(:))); 
           
           PredictionImBin(PredictionImBin < BinaryThreshold) = 0;
           PredictionImBin(PredictionImBin > 0) = 1;           

           TruthBin2 = CreateBinMask(imfill(logical(TruthImBin),'holes'));
           PredictionBin2 = CreateBinMask(imfill(logical(PredictionImBin),'holes'));           
             
           DiceTemp = dice(TruthBin2,PredictionBin2);
           DICE = [DICE;DiceTemp];
           MSDNorm = [MSDNorm;MSD(TruthBin2, PredictionBin2)]; 
           
           PredictionProps = regionprops(PredictionBin2,'Centroid','Area','BoundingBox');
           
           BadProps = [];
           
           for j = 1:numel(PredictionProps)
               if PredictionProps(j).Area < AreaThresh
                   BadProps = [BadProps;j];
               end
           end
           
           PredictionProps(BadProps) = [];
           
           %Open real Unwarped File;
           try
              RealUnwarpedIm = imread(fullfile(UnWarpedFile,[num2str([ImNum,1],str),'.png']));
           catch
              error(['File ',fullfile(UnWarpedFile,[num2str([ImNum,1],str),'.png']),' cannot be read'])
           end

           [s1,s2] = size(RealUnwarpedIm);

           GTVCropped = RealUnwarpedIm(:,round(s2/2)+1:end);

           W = round(s1/2 - 512/2);

           GTVCropped = GTVCropped(W+1:end-W,W+1:end-W);

           RealPropsXCorr = regionprops(imbinarize(GTVCropped,0),'Centroid','Area','BoundingBox');
           
           area = zeros(1,numel(RealPropsXCorr));
           for j = 1:numel(RealPropsXCorr)
               area(j) = RealPropsXCorr(j).Area;
           end         

           if numel(PredictionProps) > 1
                PredImNew = zeros(size(PredictionBin2),'like',PredictionBin2);
                for j = 1:numel(PredictionProps)
                    bboxLoop = round(cat(1,PredictionProps(j).BoundingBox));
                    PredImNew(bboxLoop(2):bboxLoop(2)+bboxLoop(4)-1,bboxLoop(1):bboxLoop(1)+bboxLoop(3)-1) = PredictionBin2(bboxLoop(2):bboxLoop(2)+bboxLoop(4)-1,bboxLoop(1):bboxLoop(1)+bboxLoop(3)-1);
                end
                
                [r,c] = find(PredImNew);
               PredictedCentroid = [PredictedCentroid;[mean(c) mean(r)]];
           else
               ind = 1;
               PredictedCentroid = [PredictedCentroid;PredictionProps.Centroid];
           end

           PredictionImNew = zeros(size(PredictionIm),'like',PredictionIm);

           bbox = round(cat(1,PredictionProps(ind).BoundingBox));

           PredictionImNew(bbox(2):bbox(2)+bbox(4)-1,bbox(1):bbox(1)+bbox(3)-1) = PredictionIm(bbox(2):bbox(2)+bbox(4)-1,bbox(1):bbox(1)+bbox(3)-1);

           RealProps = regionprops(TruthBin2,'Centroid','Area','BoundingBox');

           if numel(RealProps) > 1
               Area = zeros(numel(RealProps),1);
               for j = 1:numel(Area)
                   Area(j) = RealProps(j).Area;
               end
               [~,RealInd] = max(Area);
               RealCentroid = [RealCentroid;RealProps(RealInd).Centroid];
           else
               RealCentroid = [RealCentroid;RealProps.Centroid];
               RealInd = 1;
           end  
            
            %% Unwarped tumour location/segmentation
            
            UnwarpedImBin = imbinarize(GTVCropped,BinaryThreshold);
            
            UnWarpedDICE = [UnWarpedDICE;dice(TruthBin2,UnwarpedImBin)];
            
            UnwarpedProps = regionprops(UnwarpedImBin,'Centroid','Area');
            
            areasUnwarped = cat(1, UnwarpedProps.Area);
            CentroidUnwarped = cat(1, UnwarpedProps.Centroid);
            
            [~, idxUnwarped] = max(areasUnwarped);
            UnWarpedCentroid = [UnWarpedCentroid;CentroidUnwarped(idxUnwarped,:)]; 
            
            UnWarpedMSD = [UnWarpedMSD;MSD(TruthBin2,UnwarpedImBin)];
            
            %% Plotting

            imshow(ActualIm,[],'parent',ax1,'border','tight')
            hold(ax1,'on')
            BTruth = bwboundaries(TruthBin2);
            BPredict = bwboundaries(PredictionBin2);
            BUnWarped = bwboundaries(UnwarpedImBin);
            
            %Plot ground truth
            visboundaries(ax1,BTruth,'Color','g')            
            plot(ax1,RealProps(RealInd).Centroid(1),RealProps(RealInd).Centroid(2),'g+')
            %Plots predicted segmentation
            visboundaries(ax1,BPredict,'Color','r')
            plot(ax1,PredictedCentroid(end,1),PredictedCentroid(end,2),'r+')
            %Plots unwarped data
            visboundaries(ax1,BUnWarped,'Color','b')
            plot(ax7,UnWarpedCentroid(end,1),UnWarpedCentroid(end,2),'b+')
            
            %Write frames to video file
            if videoFlag
                frame = getframe(fig1);
                sFrame = size(frame.cdata);
                if sFrame(1) ~= 512 || sFrame(2) ~= 512
                    frame.cdata = imresize(frame.cdata,[512 512]);
                end
                if ModelNum == 12
                    writeVideo(video,frame); 
                elseif ModelNum == 13
                    writeVideo(video3,frame); 
                elseif ModelNum == 14
                    writeVideo(video3,frame); 
                elseif ModelNum == 15 
                    writeVideo(video3,frame); 
                end               
            end
            
            hold(ax1,'off')
           
       else
        disp(['File does not exist: ',realstr])
       end
    end
    
    %Close video files
    if videoFlag
        close(video)
        close(video2)
        close(video3)
        close(video4)
    end
    
    angles = 0:359;
    
    %% Plot results

    %Results for normal centroid matching
    figure()
    boxplot(DICE,'Whisker',5)
    title('DICE - Norm')
    
    figure()
    boxplot((RealCentroid(:,1) - PredictedCentroid(:,1))*pixelpermm,'Whisker',5)
    title('X Centroid Error - Norm')    
    ylabel('Error (mm)')

    figure()
    boxplot((RealCentroid(:,2) - PredictedCentroid(:,2))*pixelpermm,'Whisker',5)
    title('Y Centroid Error - Norm')   
    ylabel('Error (mm)')    
    
    Dist = sqrt((RealCentroid(:,1) - PredictedCentroid(:,1)).^2+(RealCentroid(:,2) - PredictedCentroid(:,2)).^2)*pixelpermm;
    
    disp('No X Correlation:')
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICE)),'(',num2str(std(DICE)),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,1) - PredictedCentroid(:,1))*pixelpermm)),'(',num2str(std((RealCentroid(:,1) - PredictedCentroid(:,1))*pixelpermm)),')'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,2) - PredictedCentroid(:,2))*pixelpermm)),'(',num2str(std((RealCentroid(:,2) - PredictedCentroid(:,2))*pixelpermm)),')'])
    disp(['total distance [mean(SD)] = ',num2str(mean(Dist)),'(',num2str(std(Dist)),')'])
    disp(['Mean Surface Distance(mm) [mean(SD)] = ',num2str(mean(MSDNorm)*pixelpermm),'(',num2str(std(MSDNorm)*pixelpermm),')'])

    DistUnwarped = sqrt((RealCentroid(:,1) - UnWarpedCentroid(:,1)).^2+(RealCentroid(:,2) - UnWarpedCentroid(:,2)).^2)*pixelpermm;
    
    disp('With UnWarped Results:')
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(UnWarpedDICE)),'(',num2str(std(UnWarpedDICE)),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,1) - UnWarpedCentroid(:,1))*pixelpermm)),'(',num2str(std((RealCentroid(:,1) - UnWarpedCentroid(:,1))*pixelpermm)),')'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,2) - UnWarpedCentroid(:,2))*pixelpermm)),'(',num2str(std((RealCentroid(:,2) - UnWarpedCentroid(:,2))*pixelpermm)),')'])
    disp(['total distance [mean(SD)] = ',num2str(mean(DistUnwarped)),'(',num2str(std(DistUnwarped)),')'])     
    disp(['Mean Surface Distance(mm) [mean(SD)] = ',num2str(mean(UnWarpedMSD)*pixelpermm),'(',num2str(std(UnWarpedMSD)*pixelpermm),')'])
    
    figure()
    subplot(1,4,1)
    polarplot(deg2rad(TestAngles),abs(RealCentroid(:,1) - PredictedCentroid(:,1))*pixelpermm,...
        deg2rad(TestAngles),abs(RealCentroid(:,1) - UnWarpedCentroid(:,1))*pixelpermm)
    subplot(1,4,2)
    polarplot(deg2rad(TestAngles),abs(RealCentroid(:,2) - PredictedCentroid(:,2))*pixelpermm,...
        deg2rad(TestAngles),abs(RealCentroid(:,2) - UnWarpedCentroid(:,2))*pixelpermm)
    subplot(1,4,3)
    polarplot(deg2rad(TestAngles),Dist,deg2rad(TestAngles),DistUnwarped)
    subplot(1,4,4)
    polarplot(deg2rad(TestAngles),DICE,deg2rad(TestAngles),UnWarpedDICE)
    
    %% Save results as .mat file

    save(fullfile(OutputDir,'results.mat'),'RealCentroid','PredictedCentroid','UnWarpedCentroid',...
        'DICE','UnWarpedDICE','Dist','DistUnwarped','angles','MSDNorm','UnWarpedMSD','VolumeFlag')
    
end

function[Im2] = CreateBinMask(Im)

    Props = regionprops(Im);

    Im2 = Im;
    
    if numel(Props) > 1
       Area = zeros(numel(Props),1);
       for j = 1:numel(Area)
           Area(j) = Props(j).Area;
       end
       [~,ind] = max(Area);
       
       BBox = Props(ind).BoundingBox;
       
       xi = [BBox(1) BBox(1) BBox(1) + BBox(3) BBox(1) + BBox(3)];

       yi = [BBox(2) BBox(2) + BBox(4) BBox(2) + BBox(4) BBox(2)];

       [s1,s2] = size(Im);
       
       M = poly2mask(xi,yi,s1,s2);
       
       Im2(M == 0) = 0;

    end
    
end