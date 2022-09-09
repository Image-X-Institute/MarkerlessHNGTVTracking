function DeepLearningRandEvaluation(BaseFile,ResultsFile,UnWarpedFile,str,StartingVal)

    %TODO
    %Add MSD calculation

    close('all')

    %ResultsFile = 'Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\HNSCC-01-0001\12-05-1998-RT SIMULATION-43582\DeepLearning\Checkpoints\HNCMarkerless_Patient003\web\images';
    %ResultsFile = 'Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\HNSCC-01-0001\12-05-1998-RT SIMULATION-43582\DeepLearning\Results\HNCMarkerless_Patient001_RandomCrop\WarpedData\images';
    %ResultsFile = 'Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\HNSCC-01-0001\12-05-1998-RT SIMULATION-43582\DeepLearning\Results\HNCMarkerless_Patient001_MultiVol\Warpedimages';
    %ResultsFile = 'Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\HNSCC-01-0001\12-05-1998-RT SIMULATION-43582\DeepLearning\Results\HNCMarkerless_Patient001_MultiVol\UnWarpedimages\images';
    
    %ResultsFile = fullfile(BaseFile,'DeepLearningNew/results/images');
    %ResultsFile = fullfile(BaseFile,'DeepLearningNew2/results/images');
    %ResultsFile = fullfile(BaseFile,'DeepLearningShift/results/images');
    %ResultsFile = fullfile(BaseFile,'DeepLearningNew3/results/images');
    %ResultsFile = fullfile(BaseFile,'DeepLearningNew3/results_2/images');
    %ResultsFile = fullfile(BaseFile,'DeepLearningNew/Stuff/images');
    
    %UnWarpedFile = fullfile(BaseFile,'DeepLearningNew4/DataSet/train');
    %WarpedFile = fullfile(BaseFile,'DeepLearningNew4/DataSet/test');
    
    %JP2UnwarpedDir = fullfile(BaseFile,'DeepLearningNew4/DRRImg');
    %JP2WarpedDir = fullfile(BaseFile,'DeepLearningNew4/WarpedDRRImg');
    
    %UnWarpedFile = fullfile(BaseFile,'DeepLearningNew/Stuff/DataSet/train');
    %UnWarpedFile = fullfile(BaseFile,'DeepLearningNew/DataSet/train');
    %UnWarpedFile = fullfile(BaseFile,'DeepLearningNew/DataSet/train');
    UnWarpedFile = fullfile(BaseFile,'DeepLearningMultiVol2','UnWarpedIms');
    %UnWarpedFile = fullfile(BaseFile,'DeepLearningMultiVol/DataSet/train');
    %UnWarpedFile = fullfile(BaseFile,'DeepLearningAlt3/DataSet/train');
    WarpedFile = fullfile(BaseFile,'DeepLearningNew/DataSet/test');
    
    JP2UnwarpedDir = fullfile(BaseFile,'DeepLearningNew/DRRImg');
    %JP2WarpedDir = fullfile(BaseFile,'DeepLearningNew2/WarpedDRRImg');
    %JP2WarpedDir = fullfile(BaseFile,'DeepLearningAlt/WarpedDRRImg');
    JP2WarpedDir = fullfile(BaseFile,'DeepLearningMultiVol2/DRRImg13');
    
    %[UnWarpedWholeImages, ~] = Read_jp2(JP2UnwarpedDir);
%     if exist(JP2WarpedDir,'dir')
%         [WarpedWholeImages, ~] = Read_jp2(JP2WarpedDir);
%         WarpedFlag = 1;
%     else
%         WarpedFlag = 0;
%     end
    WarpedFlag = 0;
    
    %str = 'epoch%03d';
    str = 'CombineMaskDRR%04d_%02d';
    %StartingVal = 2;
    StartingVal = 1;
    
    %angleFile = 'Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\GeometryFiles\AnglesOverSample15degwidth.csv';
    
    %angles = ReadAnglesFile(angleFile);
    
    TestAngles = [];
    
    angles = 0:0.1:359.9;
    
    AreaThresh = 1000;
    
    %angleFile = 'Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\VirtualAngles.csv';
    
    %Create video of results?
    videoFlag = 1;
    
    AddMargin = 0;
    
    %angles = csvread(angleFile);
    
    i = StartingVal;  
    
    flag = 1;
    
    DICE = [];
    PredictedCentroid = [];
    PredictedCentroidOther = [];
    RealCentroid = [];
    MSDNorm = [];    

    PredictedCentroidXCorr = [];
    DICEXCorr = [];
    MSDXCorr = [];
    
    
    PredictedCentroidRigid = [];
    DICERigid = [];    
    
    %Results for if you assume no motion
    UnWarpedDICE = [];
    UnWarpedCentroid = [];
    UnWarpedMSD = [];
    
    t = [];
    
%     figure()
%     ax1 = gca;
%     
%     figure()
%     ax2 = gca;
%     
%     figure()
%     ax3 = gca;    
    
    figure()
    ax4 = gca;      
    fig4 = gca;
    
%     figure()
%     ax5 = gca; 
%     
%     figure()
%     ax6 = gca;
    
    figure('Menu','none','ToolBar','none')
    ax7 = gca;
    fig7 = gca;
 
    figure('Menu','none','ToolBar','none')
    ax8 = gca;
    fig8 = gca; 
    
    figure()
    ax9 = gca;
    fig9 = gca; 
    
%     figure()
%     ax10 = gca;
%     
    %figure('Menu','none','ToolBar','none')
    figure()
    ax11 = gca;
    fig11 = gca; 
    
    pixelpermm = 0.256;
    
    FrameRate = 5;
    
    if videoFlag
        video = VideoWriter(fullfile(ResultsFile,'ResultsNormModel1'));   
        video2 = VideoWriter(fullfile(ResultsFile,'ResultsXCorrModel1'));
        video3 = VideoWriter(fullfile(ResultsFile,'ResultsNormModel2'));
        video4 = VideoWriter(fullfile(ResultsFile,'ResultsXCorrModel2'));
        %video5 = VideoWriter(fullfile(ResultsFile,'ResultsRigidReg'));
        video6 = VideoWriter(fullfile(ResultsFile,'ResultsWholeFrameModel1'));
        video.FrameRate = FrameRate;
        video2.FrameRate = FrameRate;
        video3.FrameRate = FrameRate;
        video4.FrameRate = FrameRate;
        video6.FrameRate = FrameRate;
        %video5.FrameRate = FrameRate;
        
        open(video);
        open(video2);
        open(video3);
        open(video4);
        %open(video5);
        open(video6);
    end
    
    [optimizer, metric] = imregconfig('monomodal');
    
    ImageList = lscell([ResultsFile,'/*real_A.png']);
    
    VolumeFlag = zeros(1,numel(ImageList));

    %while(flag)
    for i = 1:numel(ImageList)
        
        realstr2 = ImageList{i};
        
        [~,name,ext] = fileparts(ImageList{i});
        
        ImNum = str2double(name(15:18));
        ModelNum = str2double(name(20:21));
        
        TestAngles = [TestAngles;angles(ImNum)];
        
        str2  = [num2str([ImNum,ModelNum],str),'_'];
       
        realstr = fullfile(ResultsFile,[str2,'real_B.png']);
        fakestr = fullfile(ResultsFile,[str2,'fake_B.png']);
       
       %realstr2 = fullfile(ResultsFile,[str2,'real_A.png']);
       
       if ModelNum == 12
            VolumeFlag(i) = 0;
       elseif ModelNum == 13
            VolumeFlag(i) = 1;
       elseif ModelNum == 14
            VolumeFlag(i) = 2;
       elseif ModelNum == 15
            VolumeFlag(i) = 3;
       end
       
       %str3  = num2str(i*10,str);
       
       %realstr2 = fullfile(UnWarpedFile,[str3,'.png']);
       %Fix
       
       if exist(realstr,'file')
           %% Normal stuff 
           try
               PredictionIm = imread(fakestr); 
           catch
               disp('Error')
           end
           
           try
               TruthIm = imread(realstr); 
           catch
              disp('Error') 
           end
           try
               ActualIm = imread(realstr2); 
           catch
               disp('Error')
           end
           TruthImBin = TruthIm;
           
           if numel(find(TruthIm == 0)) < 100
               TruthImBin = TruthIm - median(TruthIm(:));
               TruthIm = TruthIm - median(TruthIm(:));
           end           
           
           TruthImBin(TruthImBin > 0) = 1;
           
           PredictionImBin = im2double(PredictionIm - median(PredictionIm(:))); 
           
           %PredictionImBin(PredictionImBin < 0.05) = 0;
           PredictionImBin(PredictionImBin < 0.1) = 0; 
           PredictionImBin(PredictionImBin > 0) = 1;           
           
           %imshow(TruthImBin,[0 1],'parent',ax1)
           %imshow(PredictionImBin,[0 1],'parent',ax2)
           
           TruthBin2 = CreateBinMask(imfill(logical(TruthImBin),'holes'));
           PredictionBin2 = CreateBinMask(imfill(logical(PredictionImBin),'holes'));
           %PredictionBin2 = CreateBinMask(logical(PredictionImBin));
           %PredictionBin2 = logical(PredictionImBin);
           
           %imshow(PredictionBin2,[0 1],'parent',ax3)
           
           %imshow(imgradient(PredictionBin2) + im2double(TruthIm)*10,[],'parent',ax4)
           %imshow(imgradient(PredictionBin2) + im2double(ActualIm)*10,[],'parent',ax4)
           
           
           if videoFlag
               %frame = getframe(fig4);
               %writeVideo(video,frame);             
           end
%            if i > 180
%                disp('Something')
%            end
           
           

           imshow(imgradient(PredictionBin2) + im2double(ActualIm)*10,[],'parent',ax4)
           if videoFlag
               %frame = getframe(fig4);
               %writeVideo(video2,frame);
           end
             
           DiceTemp = dice(TruthBin2,PredictionBin2);
           if isempty(DiceTemp)
               disp('error')
           end
           DICE = [DICE;DiceTemp];
           MSDNorm = [MSDNorm;MSD(TruthBin2, PredictionBin2)]; 
%            if MSD(TruthBin2, PredictionBin2) > 10
%               disp('something') 
%            end
           
           PredictionProps = regionprops(PredictionBin2,'Centroid','Area','BoundingBox');
           
           BadProps = [];
           
           for j = 1:numel(PredictionProps)
               if PredictionProps(j).Area < AreaThresh
                   BadProps = [BadProps;j];
               end
           end
           
           PredictionProps(BadProps) = [];
           
           %PredictionProps = regionprops(PredictionImBin,'Centroid','Area','BoundingBox');
           
           %Double check the angles
           %[~,angleInd] = min(abs(angles - ImNum-1));
           
           %Open real Unwarped File;
           try
              RealUnwarpedIm = imread(fullfile(UnWarpedFile,[num2str([ImNum,1],str),'.png']));
           catch
              disp('error') 
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
           [~,RealXCorrInd] = max(area);          

           if numel(PredictionProps) > 1
%                Area = zeros(numel(PredictionProps),1);
                PredImNew = zeros(size(PredictionBin2),'like',PredictionBin2);
                for j = 1:numel(PredictionProps)
%                    Area(j) = PredictionProps(j).Area;
                    bboxLoop = round(cat(1,PredictionProps(j).BoundingBox));
                    PredImNew(bboxLoop(2):bboxLoop(2)+bboxLoop(4)-1,bboxLoop(1):bboxLoop(1)+bboxLoop(3)-1) = PredictionBin2(bboxLoop(2):bboxLoop(2)+bboxLoop(4)-1,bboxLoop(1):bboxLoop(1)+bboxLoop(3)-1);
                end
                
                [r,c] = find(PredImNew);
                
%                [~,ind] = max(Area);
%                PredictedCentroid = [PredictedCentroid;PredictionProps(ind).Centroid];
               PredictedCentroid = [PredictedCentroid;[mean(c) mean(r)]];
               PredictedCentroidOther = [PredictedCentroidOther;DoWeirdCentroidStuff(RealPropsXCorr(RealXCorrInd),GTVCropped,PredictionProps(ind))];
               %PredictedCentroidOther = [PredictedCentroidOther;PredicitionProps(ind).Centroid];
           else
               ind = 1;
               PredictedCentroid = [PredictedCentroid;PredictionProps.Centroid];
               PredictedCentroidOther = [PredictedCentroidOther;DoWeirdCentroidStuff(RealPropsXCorr(RealXCorrInd),GTVCropped,PredictionProps)];
               %PredictedCentroidOther = [PredictedCentroidOther;PredicitionProps.Centroid];
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

           %% Template matching
           %Edit to check multiple projections and see if that increase
           %accuracy

           %Open real Warped File; 
           %RealwarpedIm = imread(fullfile(WarpedFile,[num2str(i,str),'.png']));
           
           %Open real Warped Cropped File;
           %ActualIm
           
%            disp(fullfile(UnWarpedFile,[num2str(i*10,str),'.png']))
%            disp(fullfile(WarpedFile,[num2str(i*10,str),'.png']))
%            disp(realstr)
%            disp(realstr2)
           
           %[DRRCropped,GTVCropped] = CropUnWarpedIm(RealUnwarpedIm,RealwarpedIm,ActualIm);
           
           tic()
           RealPropsXCorr = regionprops(imbinarize(GTVCropped,0),'Centroid','Area','BoundingBox');

           area = zeros(1,numel(RealPropsXCorr));
           for j = 1:numel(RealPropsXCorr)
               area(j) = RealPropsXCorr(j).Area;
           end
           [~,RealXCorrInd] = max(area);
           
           
%            if numel(RealPropsXCorr) > 0
%            end
%            else 
%                
%            end
           
            bbox = round(cat(1, RealPropsXCorr(RealXCorrInd).BoundingBox));
            try
                template = GTVCropped(bbox(2):bbox(2)+bbox(4)-1,bbox(1):bbox(1)+bbox(3)-1);
            catch
                disp('Something')
            end
            if bbox(2)+bbox(4)-1 > 512
                disp('Something')
            end
            if bbox(1)+bbox(3)-1 > 512
                disp('Something')
            end            
            %bbox = round(cat(1, PredicitionProps(1).BoundingBox));
            %template = PredictionIm(bbox(2):bbox(2)+bbox(4)-1,bbox(1):bbox(1)+bbox(3)-1);
            
            %bx = size(PredictionIm, 2); 
            %by = size(PredictionIm, 1);            
            %tx = size(template, 2); % used for bbox placement
            %ty = size(template, 1);

            %bx = size(TruthIm, 2); 
            %by = size(TruthIm, 1);
            
            %Ga = fft2(PredictionIm);

            %Gb = fft2(template,bx,by);
            
            %c = normxcorr2(template,PredictionIm);
            c = normxcorr2(template,PredictionImNew);
            %c = real(ifft2((Ga.*conj(Gb))./abs(Ga.*conj(Gb))));
            
            [ypeak,xpeak] = find(c==max(c(:)));

            %yoffSet = ypeak;
            %xoffSet = xpeak;            
            %if isempty(xpeak)
                %c = normxcorr2(template,PredictionIm);
                %[ypeak,xpeak] = find(c==max(c(:)));
            yoffSet = ypeak-size(template,1);
            xoffSet = xpeak-size(template,2); 
            if yoffSet < 1
                yoffSet = 1;
            end
            if xoffSet < 1
                xoffSet = 1;
            end                
            %end
                        
            % Create the new segmentation
            newFake = zeros(size(TruthIm));
            xInds = xoffSet:xoffSet+size(template,2)-1;
            yInds = yoffSet:yoffSet+size(template,1)-1;
            if yInds(end) > 512
                yInds = yoffSet:512;
            end
            if xInds(end) > 512
                xInds = xoffSet:512;
            end            
            newFake(yInds,xInds) = template(1:numel(yInds),1:numel(xInds));
            
            %imshowpair(newFake,TruthIm,'montage','parent',ax5)
            newFake = imbinarize(newFake,0.1); 
            %newFake = imbinarize(newFake,0.05); 
            
            t = [t;toc];
            
            try
            %imshow(newFake + TruthBin2,[0 2],'parent',ax6)
            catch
                disp('Something')
            end
            % Determine the fake centroid
            s = regionprops(newFake,'Centroid','Area');    
            areas = cat(1, s.Area);
            centroids = cat(1, s.Centroid);
            [m, idx] = max(areas);
            PredictedCentroidXCorr = [PredictedCentroidXCorr;centroids(idx,:)];

            % Calculate the DICE
            DICEXCorr = [DICEXCorr;dice(TruthBin2, newFake)]; 
            
            MSDXCorr = [MSDXCorr;MSD(TruthBin2, newFake)]; 
            
            %% Registration of GTV
            
            moving_reg = imregister(GTVCropped,PredictionImNew,'Rigid',optimizer,metric);
            
            %imshowpair(moving_reg,TruthIm,'montage','parent',ax10)
            
            RegBin = imbinarize(moving_reg,0);
            
            RigidProps = regionprops(RegBin,'Centroid','Area'); 
            
            %Get Centroid Location
            areasRigid = cat(1, RigidProps.Area);
            centroidsRigid = cat(1, RigidProps.Centroid);
            
            [~, idxRigid] = max(areasRigid);
            PredictedCentroidRigid = [PredictedCentroidRigid;centroidsRigid(idxRigid,:)];            
            
            % Calculate the DICE
            DICERigid = [DICERigid;dice(TruthBin2, RegBin)];             
            
            %% Unwarped tumour location/segmentation
            
            UnwarpedImBin = imbinarize(GTVCropped,0.1);
            
            if AddMargin
                Margin = 5/pixelpermm;
                UnwarpedImBin = AddMarginToSegmentation(UnwarpedImBin,Margin);
            end
            
            UnWarpedDICE = [UnWarpedDICE;dice(TruthBin2,UnwarpedImBin)];
            
            UnwarpedProps = regionprops(UnwarpedImBin,'Centroid','Area');
            
            areasUnwarped = cat(1, UnwarpedProps.Area);
            CentroidUnwarped = cat(1, UnwarpedProps.Centroid);
            
            [~, idxUnwarped] = max(areasUnwarped);
            UnWarpedCentroid = [UnWarpedCentroid;CentroidUnwarped(idxUnwarped,:)]; 
            
            UnWarpedMSD = [UnWarpedMSD;MSD(TruthBin2,UnwarpedImBin)];
            
            %% Plotting
            
            imshow(ActualIm,[],'parent',ax7,'border','tight')
            hold(ax7,'on')
            BTruth = bwboundaries(TruthBin2);
            BPredict = bwboundaries(PredictionBin2);
            BUnWarped = bwboundaries(UnwarpedImBin);
            
            visboundaries(ax7,BTruth,'Color','g')
            plot(ax7,RealProps(RealInd).Centroid(1),RealProps(RealInd).Centroid(2),'g+')
            %plot(ax7,PCentroid(1),RealProps(RealInd).Centroid(2),'g+')
            try
            visboundaries(ax7,BPredict,'Color','r')
            %plot(ax7,PredictionProps(ind).Centroid(1),PredictionProps(ind).Centroid(2),'r+')
            plot(ax7,PredictedCentroid(end,1),PredictedCentroid(end,2),'r+')
            catch
               disp('Something') 
            end
            %Plots unwarped data
            visboundaries(ax7,BUnWarped,'Color','b')
            plot(ax7,UnWarpedCentroid(end,1),UnWarpedCentroid(end,2),'b+')

            imshow(ActualIm,[],'parent',ax8,'border','tight')
            hold(ax8,'on')
            BTruth = bwboundaries(TruthBin2);
            BPredict = bwboundaries(newFake);
            BUnWarped = bwboundaries(UnwarpedImBin);
           
            visboundaries(ax8,BTruth,'Color','g')
            plot(ax8,RealProps(RealInd).Centroid(1),RealProps(RealInd).Centroid(2),'g+')
            try
            visboundaries(ax8,BPredict,'Color','r')
            plot(ax8,centroids(idx,1),centroids(idx,2),'r+')
            catch
               disp('Something') 
            end    
            %Plots unwarped data
            %visboundaries(ax8,BUnWarped,'Color','b')
            %plot(ax8,UnWarpedCentroid(end,1),UnWarpedCentroid(end,2),'b+')
            
            if WarpedFlag
                %Add unwarped contour to image
                C2 = normxcorr2(ActualIm,WarpedWholeImages(:,:,ImNum));
                
                [ypeak2,xpeak2] = find(C2==max(C2(:)));
                
                yoffSet2 = ypeak2 - size(ActualIm,1);
                xoffSet2 = xpeak2 - size(ActualIm,2);
                
                %Need to flip so head is at the top of image because people are
                %stupid. 
                %imshow(imnoise(flipud(WarpedWholeImages(:,:,i)'),'Poisson'),[],'parent',ax11,'border','tight')
                imshow(imnoise(fliplr(WarpedWholeImages(:,:,ImNum)'),'Poisson'),[],'parent',ax11,'border','tight')
                hold(ax11,'on')
                try               
                    BTruth2 = BTruth;
                    if numel(BTruth) > 1
                        BTruth2(2:numel(BTruth)) = [];
                    end
                    BTruth2{:,1} = BTruth2{:,1} + [yoffSet2,xoffSet2];
                    %BTruth2{:,1} = [1024-BTruth2{:,1}(:,2) 768-BTruth2{:,1}(:,1)];
                    BTruth2{:,1} = [BTruth2{:,1}(:,2) 1024-BTruth2{:,1}(:,1)];                    
                    visboundaries(ax11,BTruth2,'Color','g')
                catch
                    visboundaries(ax11,BTruth,'Color','g')
                end
                %plot(ax11,[RealProps(RealInd).Centroid(2) + yoffSet2],1024-[RealProps(RealInd).Centroid(1) + xoffSet2],'g+')
                %plot(ax11,768-[RealProps(RealInd).Centroid(2) + yoffSet2],[RealProps(RealInd).Centroid(1) + xoffSet2],'g+')
                plot(ax11,1024-[RealProps(RealInd).Centroid(2) + yoffSet2],[RealProps(RealInd).Centroid(1) + xoffSet2],'g+')
                try
                    BPredict2 = BPredict;
                    BPredict2{:,1} = BPredict2{:,1} + [yoffSet2,xoffSet2];
                    %BPredict2{:,1} = [1024-BPredict2{:,1}(:,2) BPredict2{:,1}(:,1)];
                    %BPredict2{:,1} = [BPredict2{:,1}(:,2) 768-BPredict2{:,1}(:,1)];
                    BPredict2{:,1} = [BPredict2{:,1}(:,2) 1024-BPredict2{:,1}(:,1)];
                    visboundaries(ax11,BPredict2,'Color','r','LineStyle','--','EnhanceVisibility',false)
                catch
                    visboundaries(ax11,BPredict,'Color','r')
                end
                %plot(ax11,[centroids(idx,2) + yoffSet2],1024-[centroids(idx,1)+ xoffSet2],'r+')
                %plot(ax11,768-[centroids(idx,2) + yoffSet2],[centroids(idx,1)+ xoffSet2],'r+')
                plot(ax11,1024-[centroids(idx,2) + yoffSet2],[centroids(idx,1)+ xoffSet2],'r+')
                
                %Plots unwarped data
                try
                    BUnWarped2 = BUnWarped;
                    BUnWarped2{:,1} = BUnWarped2{:,1} + [yoffSet2,xoffSet2];
                    %BPredict2{:,1} = [1024-BPredict2{:,1}(:,2) BPredict2{:,1}(:,1)];
                    %BUnWarped2{:,1} = [BUnWarped2{:,1}(:,2) 768-BUnWarped2{:,1}(:,1)];
                    BUnWarped2{:,1} = [BUnWarped2{:,1}(:,2) 1024-BUnWarped2{:,1}(:,1)];
                    visboundaries(ax11,BUnWarped2,'Color','b','LineStyle','--','EnhanceVisibility',false)
                catch
                    visboundaries(ax11,BUnWarped2,'Color','b')
                end
                %plot(ax11,[centroids(idx,2) + yoffSet2],1024-[centroids(idx,1)+ xoffSet2],'r+')
                %plot(ax11,768-[UnWarpedCentroid(end,2) + yoffSet2],[UnWarpedCentroid(end,1)+ xoffSet2],'b+')      
                plot(ax11,1024-[UnWarpedCentroid(end,2) + yoffSet2],[UnWarpedCentroid(end,1)+ xoffSet2],'b+')  
                
                hold(ax11,'off')
                
            end
            
            imshow(ActualIm,[],'parent',ax9)
            hold(ax9,'on')
            BTruth = bwboundaries(TruthBin2);
            BPredict = bwboundaries(RegBin);
           
            visboundaries(ax9,BTruth,'Color','g')
            plot(ax9,RealProps(RealInd).Centroid(1),RealProps(RealInd).Centroid(2),'g+')
            try
            visboundaries(ax9,BPredict,'Color')
            plot(ax9,centroidsRigid(idxRigid,1),centroidsRigid(idxRigid,2),'r+')
            catch
               disp('Something') 
            end 
            
            if videoFlag
                frame = getframe(fig7);
                sFrame = size(frame.cdata);
                if sFrame(1) ~= 512 || sFrame(2) ~= 512
                    frame.cdata = imresize(frame.cdata,[512 512]);
                end
                try
                    if ModelNum == 12
                        writeVideo(video,frame); 
                    else
                        writeVideo(video3,frame); 
                    end
                catch
                    disp('Error')
                end
                frame = getframe(fig8);
                try
                    if ModelNum == 12
                        writeVideo(video2,frame);
                    else
                        writeVideo(video4,frame);
                    end
                catch
                   disp('Error') 
                end
                frame = getframe(fig11);
                try
                    if ModelNum == 12
                        writeVideo(video6,frame);
                    end
                catch
                    disp('Error')
                end
                %frame = getframe(fig9);
                %writeVideo(video5,frame);                 
            end
            
            hold(ax7,'off')
            hold(ax8,'off')
            hold(ax9,'off')
            

            if norm(PredictedCentroid(end,:)-RealCentroid(end,:))*pixelpermm > 3
                disp('Something')
            end
%             if(ImNum == 1 || ImNum == 601 || ImNum == 1201 || ImNum == 2401) && ModelNum == 12
%                 disp('Something')
%             end
            
%             if (i == 1) || (i == 60*2-1) || (i == 120*2-1) || (i == 240*2-1)
%                disp('Something') 
%             end
%             if i == 76 || i == 235
%                disp('Something') 
%             end
            
%             if i == 87 || i == 197 || i == 273
%                 disp('Something')
%             end

           %i = i+1;   
           
       else
           flag = 0;
       end
    end
    
    if videoFlag
        close(video)
        close(video2)
        close(video3)
        try
            close(video4)
        catch
            disp('Error')
        end
        %close(video5)
        close(video6)
    end
    %angles = linspace(0,360,numel(RealCentroid(:,1)));
    
%     VolumeFlag = zeros(1,numel(ImageList));
%     VolumeFlag(1:2:end) = 1;
    
    t = t + 0.005;
    
    disp(mean(t))
    
    %angles = csvread(angleFile);
    
    angles = 0:359;
    
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
    
    Dist = sqrt((RealCentroid(:,1) - PredictedCentroidOther(:,1)).^2+(RealCentroid(:,2) - PredictedCentroidOther(:,2)).^2)*pixelpermm;  

    disp('No X Correlation - other centroid metric:')
    %disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICE)),'(',num2str(std(DICE)),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,1) - PredictedCentroidOther(:,1))*pixelpermm)),'(',num2str(std((RealCentroid(:,1) - PredictedCentroidOther(:,1))*pixelpermm)),')'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,2) - PredictedCentroidOther(:,2))*pixelpermm)),'(',num2str(std((RealCentroid(:,2) - PredictedCentroidOther(:,2))*pixelpermm)),')'])
    disp(['total distance [mean(SD)] = ',num2str(mean(Dist)),'(',num2str(std(Dist)),')'])

    DistXCorr = sqrt((RealCentroid(:,1) - PredictedCentroidXCorr(:,1)).^2+(RealCentroid(:,2) - PredictedCentroidXCorr(:,2)).^2)*pixelpermm;
    
    disp('With X Correlation:')
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICEXCorr)),'(',num2str(std(DICEXCorr)),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,1) - PredictedCentroidXCorr(:,1))*pixelpermm)),'(',num2str(std((RealCentroid(:,1) - PredictedCentroidXCorr(:,1))*pixelpermm)),')'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,2) - PredictedCentroidXCorr(:,2))*pixelpermm)),'(',num2str(std((RealCentroid(:,2) - PredictedCentroidXCorr(:,2))*pixelpermm)),')'])
    disp(['total distance [mean(SD)] = ',num2str(mean(DistXCorr)),'(',num2str(std(DistXCorr)),')'])    
    disp(['Mean Surface Distance(mm) [mean(SD)] = ',num2str(mean(MSDXCorr)*pixelpermm),'(',num2str(std(MSDXCorr)*pixelpermm),')'])
    
    DistRigid = sqrt((RealCentroid(:,1) - PredictedCentroidRigid(:,1)).^2+(RealCentroid(:,2) - PredictedCentroidRigid(:,2)).^2)*pixelpermm;
    
    disp('With Rigid Registration:')
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICERigid)),'(',num2str(std(DICERigid)),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,1) - PredictedCentroidRigid(:,1))*pixelpermm)),'(',num2str(std((RealCentroid(:,1) - PredictedCentroidRigid(:,1))*pixelpermm)),')'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,2) - PredictedCentroidRigid(:,2))*pixelpermm)),'(',num2str(std((RealCentroid(:,2) - PredictedCentroidRigid(:,2))*pixelpermm)),')'])
    disp(['total distance [mean(SD)] = ',num2str(mean(DistRigid)),'(',num2str(std(DistRigid)),')'])        
    
    DistUnwarped = sqrt((RealCentroid(:,1) - UnWarpedCentroid(:,1)).^2+(RealCentroid(:,2) - UnWarpedCentroid(:,2)).^2)*pixelpermm;
    
    disp('With UnWarped Results:')
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(UnWarpedDICE)),'(',num2str(std(UnWarpedDICE)),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,1) - UnWarpedCentroid(:,1))*pixelpermm)),'(',num2str(std((RealCentroid(:,1) - UnWarpedCentroid(:,1))*pixelpermm)),')'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(:,2) - UnWarpedCentroid(:,2))*pixelpermm)),'(',num2str(std((RealCentroid(:,2) - UnWarpedCentroid(:,2))*pixelpermm)),')'])
    disp(['total distance [mean(SD)] = ',num2str(mean(DistUnwarped)),'(',num2str(std(DistUnwarped)),')'])     
    disp(['Mean Surface Distance(mm) [mean(SD)] = ',num2str(mean(UnWarpedMSD)*pixelpermm),'(',num2str(std(UnWarpedMSD)*pixelpermm),')'])
    
    %Results using Xcorr matching
    figure()
    boxplot(DICEXCorr,'Whisker',5)
    title('DICE - xcorr')
    
    figure()
    boxplot((RealCentroid(:,1) - PredictedCentroidXCorr(:,1))*pixelpermm,'Whisker',5)
    title('X Centroid Error - xcorr')    
    ylabel('Error (mm)')

    figure()
    boxplot((RealCentroid(:,2) - PredictedCentroidXCorr(:,2))*pixelpermm,'Whisker',5)
    title('Y Centroid Error - xcorr')   
    ylabel('Error (mm)')      
       
    
    
    Dist2 = sqrt((RealCentroid(:,1) - PredictedCentroidXCorr(:,1)).^2+(RealCentroid(:,2) - PredictedCentroidXCorr(:,2)).^2)*pixelpermm;
    
    
    figure()
    boxplot(DICERigid,'Whisker',5)
    title('DICE - RigidReg')
    
    figure()
    boxplot((RealCentroid(:,1) - PredictedCentroidRigid(:,1))*pixelpermm,'Whisker',5)
    title('X Centroid Error - RigidReg')    
    ylabel('Error (mm)')

    figure()
    boxplot((RealCentroid(:,2) - PredictedCentroidRigid(:,2))*pixelpermm,'Whisker',5)
    title('Y Centroid Error - RigidReg')   
    ylabel('Error (mm)')      
       
    Dist3 = sqrt((RealCentroid(:,1) - PredictedCentroidRigid(:,1)).^2+(RealCentroid(:,2) - PredictedCentroidRigid(:,2)).^2)*pixelpermm;
    
    angles2 = angles;
    %angles2(angles2 > 200) = angles2(angles2 > 200) - 360;
    
%     figure()
%     plot(angles2,Dist,'*')
%     xlabel('Projection Angle (deg)')
%     ylabel('Centroid Distance Error (mm)')
%     
%     figure()
%     plot(angles2,Dist2,'*')
%     xlabel('Projection Angle (deg)')
%     ylabel('Centroid Distance Error (mm)')    
%     
%     figure()
%     plot(angles2,Dist3,'*')
%     xlabel('Projection Angle (deg)')
%     ylabel('Centroid Distance Error (mm)') 
    
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

    figure()
    subplot(1,4,1)
    polarplot(deg2rad(TestAngles),abs(RealCentroid(:,1) - PredictedCentroidXCorr(:,1))*pixelpermm,...
        deg2rad(TestAngles),abs(RealCentroid(:,1) - UnWarpedCentroid(:,1))*pixelpermm)
    subplot(1,4,2)
    polarplot(deg2rad(TestAngles),abs(RealCentroid(:,2) - PredictedCentroidXCorr(:,2))*pixelpermm,...
        deg2rad(TestAngles),abs(RealCentroid(:,2) - UnWarpedCentroid(:,2))*pixelpermm)
    subplot(1,4,3)
    polarplot(deg2rad(TestAngles),Dist2,deg2rad(TestAngles),DistUnwarped)
    subplot(1,4,4)
    polarplot(deg2rad(TestAngles),DICEXCorr,deg2rad(TestAngles),UnWarpedDICE)
    
    [ResultsDir,name] = fileparts(ResultsFile);
    
    save(fullfile(ResultsDir,'results.mat'),'RealCentroid','PredictedCentroid','PredictedCentroidOther','PredictedCentroidXCorr','PredictedCentroidRigid','UnWarpedCentroid',...
        'DICE','DICEXCorr','DICERigid','UnWarpedDICE','Dist','Dist2','Dist3','DistUnwarped','angles2','MSDNorm','MSDXCorr','UnWarpedMSD','VolumeFlag')
    
end

function[OtherCentroid] = DoWeirdCentroidStuff(XCorrProps,GTVCroppedIm,PredictionProps)
%Resize segmentation to be the same size as DRR segmentation

    GTVImOnly = GTVCroppedIm(round(XCorrProps.BoundingBox(2):XCorrProps.BoundingBox(2)+XCorrProps.BoundingBox(4)),...
        round(XCorrProps.BoundingBox(1):XCorrProps.BoundingBox(1)+XCorrProps.BoundingBox(3)));

    GTVScaled = imresize(GTVImOnly,[PredictionProps.BoundingBox(4) PredictionProps.BoundingBox(3)]);

    NewIm = zeros(size(GTVCroppedIm),'like',GTVCroppedIm);

    NewIm(round(PredictionProps.BoundingBox(2):PredictionProps.BoundingBox(2)+PredictionProps.BoundingBox(4)-1),...
        round(PredictionProps.BoundingBox(1):PredictionProps.BoundingBox(1)+PredictionProps.BoundingBox(3)-1)) = GTVScaled;

    NewProps = regionprops(imbinarize(NewIm,0));

    MaxArea = 0;

    for i = 1:numel(NewProps)
        if NewProps(i).Area > MaxArea
            MaxInd = i;
            MaxArea = NewProps(i).Area;
        end
    end

    OtherCentroid = NewProps(MaxInd).Centroid;

end

function[DRRCropped,GTVCropped] = CropUnWarpedIm(RealUnwarpedIm,RealwarpedIm,ActualIm)

    %Get Cropping between warped files
    RealwarpedDRR = RealwarpedIm(:,1:550);
    RealwarpedGTV = RealwarpedIm(:,551:end);
    
    RealUnwarpedDRR = RealUnwarpedIm(:,1:550);
    RealUnwarpedGTV = RealUnwarpedIm(:,551:end);    
    
    c = normxcorr2(ActualIm,RealwarpedDRR);
    [ypeak,xpeak] = find(c==max(c(:)));
    
    yoffSet = ypeak-size(ActualIm,1);
    xoffSet = xpeak-size(ActualIm,2);     
    
    try
        xInds = xoffSet:xoffSet+size(ActualIm,2)-1;
        yInds = yoffSet:yoffSet+size(ActualIm,1)-1;   
        DRRCropped = RealUnwarpedDRR(yInds,xInds);
        GTVCropped = RealUnwarpedGTV(yInds,xInds);        
    catch
       disp('Something') 
    end
    %Apply cropping to unwarped files
    %DRRCropped = zeros(size(ActualIm),'like',ActualIm);
    %GTVCropped = zeros(size(ActualIm),'like',RealUnwarpedGTV);
    
    %DRRCropped = RealUnwarpedDRR(yInds,xInds);
    %GTVCropped = RealUnwarpedGTV(yInds,xInds);
    
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