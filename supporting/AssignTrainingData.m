function AssignTrainingData(BaseFile)
%% AssignTrainingData(BaseFile)
%Assign all the training data to two testing data folders, one is made of
%random images from all volumes, and the other is all the images from the
%final volume

    %BaseFile = 'J:\CTData\HNSCC-01-0130\DeepLearningMultiVol2\DataSet';

    OriginalFile = fullfile(BaseFile,'train');

    TrainTestThresh = 0.99;   %Percentage of images that will be assinged to training vs testing
    
    RandomTrainFile = fullfile(BaseFile,'Random','train');

    if ~exist(RandomTrainFile,'dir')
        mkdir(RandomTrainFile)
    end    
    
    RandomTestFile = fullfile(BaseFile,'Random','test');
    
    if ~exist(RandomTestFile,'dir')
        mkdir(RandomTestFile)
    end

    NonRandomTrainFile = fullfile(BaseFile,'Other','train');

    if ~exist(NonRandomTrainFile,'dir')
        mkdir(NonRandomTrainFile)
    end    
    
    NonRandomTestFile = fullfile(BaseFile,'Other','test');
    
    if ~exist(NonRandomTestFile,'dir')
        mkdir(NonRandomTestFile)
    end    

    NonRandomTrainFile2 = fullfile(BaseFile,'Other2','train');

    if ~exist(NonRandomTrainFile2,'dir')
        mkdir(NonRandomTrainFile2)
    end    
    
    DataSetFile = fullfile(BaseFile,'UnWarpedIms','train');
    
    if ~exist(DataSetFile,'dir')
        mkdir(DataSetFile)
    end            
    
    VolNum = 15;
    
    imNum = 3600;
    
    cropThing = round((550-512)/2);
    
    for i = 1:VolNum
    %for i = VolNum:VolNum
       
        TestInds = randi(imNum,round(imNum*(1-TrainTestThresh)),1);     %Get indices to move from train to testing data
        
        for j = 1:imNum
            
            strName = num2str([j,i],'CombineMaskDRR%04d_%02d.png');
            if i < VolNum
                if ismember(j,TestInds)
                    im = imread(fullfile(OriginalFile,strName));
                    imCropped = CropBothImages(im,cropThing);
                    imwrite(imCropped,fullfile(RandomTestFile,strName))
                else
                    try
                        copyfile(fullfile(OriginalFile,strName),fullfile(RandomTrainFile,strName)) 
                    catch
                        disp(fullfile(OriginalFile,strName))
                    end
                end
            end
            %if (i == VolNum || i == (VolNum-1))
            if i > 11
                if mod(j,10) == 1
                    im = imread(fullfile(OriginalFile,strName));
                    imCropped = CropBothImages(im,cropThing);
                    imwrite(imCropped,fullfile(NonRandomTestFile,strName))
                end
            else
               copyfile(fullfile(OriginalFile,strName),fullfile(NonRandomTrainFile,strName))
               if i < 8
                   copyfile(fullfile(OriginalFile,strName),fullfile(NonRandomTrainFile2,strName))
               end
            end
            
            if i == 1 
                copyfile(fullfile(OriginalFile,strName),fullfile(DataSetFile,strName))
            end
            
        end
        
    end
    
end

function[imCropped] = CropBothImages(im,cropThing)

    imA = im(:,1:550);
    imB = im(:,551:end);
    
    imCropped = [imA(cropThing+1:end-cropThing,cropThing+1:end-cropThing) imB(cropThing+1:end-cropThing,cropThing+1:end-cropThing)];
    
end