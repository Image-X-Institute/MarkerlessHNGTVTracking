function DeepLearningSummaryMultiVol(Files,OutputFile)
%% DeepLearningSummaryMultiVol(Files,OutputDir)
% ------------------------------------------
% FILE   : DeepLearningSummaryMultiVol.m
% AUTHOR : Mark Gardner, The University of Sydney
% DATE   : 2022-09-12  Created.
% ------------------------------------------
% PURPOSE
%   Get the results data for each patient and combine the results together for analysis.
% ------------------------------------------
% INPUT
%   Files:          Cell Array where each index in the cell array is the
%                   directory where the 'results.mat' file was saved
%   OutputFile:     Location or filename where the output file will be
%                   saved. 
%% Input check

    for i = 1:numel(Files)
        if ~exist(Files{i},'dir')
            error(['Directory ',Files{i},' does not exist'])
        end
    end

    if nargin < 2
        OutputFile = uigetdir(matlabroot,OutputFile);
    end

    if isfolder(OutputFile)
        OutputFile = fullfile(OutputFile,'Results.csv');
    end

    close('all')

%     BaseFile = 'Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC';
% 
%     %OutputFile = 'C:\Users\mgar5380\Documents\Writing\MarkerlessHNTracking\Results.csv';
%     %OutputFile ='Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\NewResults2.csv';
%     OutputFile ='Z:\2RESEARCH\2_ProjectData\RemoveTheMask\CTData\HNSCC\NewResults3.csv';
% %     Files = {'HNSCC-01-0001\12-05-1998-RT SIMULATION-43582\DeepLearning\Results\images_Patient001',...
% %         'HNSCC-01-0002\11-20-2001-RT SIMULATION-57976\DeepLearning\results',...
% %         'HNSCC-01-0003\11-20-2001-RT SIMULATION-30560\DeepLearning\results'
% %         };
%     
%     Files = {'HNSCC-01-0001/12-05-1998-RT SIMULATION-43582',...
%     'HNSCC-01-0003/11-20-2001-RT SIMULATION-30560',...
%     'HNSCC-01-0004/08-24-1996-RT SIMULATION-72882',...
%     'HNSCC-01-0005/01-19-1998-NA-RT SIMULATION-39358',...
%     'HNSCC-01-0007/04-29-1997-RT SIMULATION-32176',...
%     'HNSCC-01-0040/12-12-1998-RT SIMULATION-35323',...
%     'HNSCC-01-0039/12-14-1998-RT SIMULATION-26503',...
%     'HNSCC-01-0070/05-10-1999-RT SIMULATION-21211',...
%     'HNSCC-01-0101/09-25-1999-NA-RT SIMULATION-93304',...
%     'HNSCC-01-0094/09-08-1999-NA-RT SIMULATION-68058',...    
%     'HNSCC-01-0002/11-20-2001-RT SIMULATION-57976',...
%     'HNSCC-01-0018/03-01-2009-NA-RT SIMULATION-16942',...    
%     'HNSCC-01-0122/04-18-2000-RT SIMULATION-20757',...
%     'HNSCC-01-0123/03-24-2000-RT SIMULATION-47852',...    
%     'HNSCC-01-0130/05-29-2000-NA-RT SIMULATION-72920'   
%     };
%     
%     BaseResultsOtherDir = 'DeepLearningMultiVol2\resultsAll';
%     %BaseResultsOtherDir = 'DeepLearningMultiVol2\resultsOther';
%     BaseResultsOneDir = 'DeepLearningMultiVol2\resultsOneVol2';
%     
%     %Generate Results
% %     for i = 1:numel(Files)
% %         slashes = strfind(Files{i},'\');
% %         try
% %             DeepLearningRandEvaluation(fullfile(BaseFile,Files{i}(1:slashes(2)-1)))
% %         catch
% %            disp(['Something went wrong for ',Files{i}]) 
% %         end
% %     end
    
    %% Define inputs

    pixelpermm = 0.256;

    DICE = [];
    RealCentroid = [];
    PredictedCentroid = [];
    UnwarpedCentroid = [];
    DICEUnwarped = [];
    Dist = [];
    DistUnwarped = [];
    angles2 = [];
    MSDNorm = [];
    MSDUnWarped = [];
    VolumeFlag = [];
    
    meanDICE = zeros(1,numel(Files));
    
    meanMSD = zeros(1,numel(Files));
    
    %for i = 1:numel(Files)

    %% Loop through files

    for i = 1:numel(Files)
        
        ResultsFile = fullfile(Files{i},'results.mat');
       
        dataOther = load(ResultsFile);

        DICE = [DICE ;[dataOther.DICE zeros(size(dataOneVol.DICE))+i]];

        VolumeFlag = [VolumeFlag; [dataOther.VolumeFlag' zeros(size(dataOther.VolumeFlag'))+i]];

        DICEUnwarped = [DICEUnwarped ;[dataOther.UnWarpedDICE zeros(size(dataOther.UnWarpedDICE))+i]];
        RealCentroid = [RealCentroid ;[dataOther.RealCentroid zeros(size(dataOther.RealCentroid(:,1)))+i]];
        PredictedCentroid = [PredictedCentroid ;[dataOther.PredictedCentroid zeros(size(dataOther.PredictedCentroid(:,1)))+i]];

        UnwarpedCentroid = [UnwarpedCentroid ;[dataOther.UnWarpedCentroid zeros(size(dataOneVol.UnWarpedCentroid(:,1)))+i]];
        Dist = [Dist;[dataOther.Dist zeros(size(dataOther.Dist))+i]];

        DistUnwarped = [DistUnwarped;[dataOther.DistUnwarped zeros(size(dataOther.DistUnwarped))+i]];
        angles2 = [angles2;[dataOther.angles2 zeros(size(dataOther.angles2))+i]];
        MSDUnWarped = [MSDUnWarped;[dataOther.UnWarpedMSD zeros(size(dataOther.UnWarpedMSD))+i]];
        MSDNorm = [MSDNorm;[dataOther.MSDNorm zeros(size(dataOther.MSDNorm))+i]];

        clear dataOther dataOneVol
    end
    
    disp('Data Loaded')
    
    Inds1 = 1:4:numel(DICE(:,1));
    Inds2 = 2:4:numel(DICE(:,1));
    Inds3 = 3:4:numel(DICE(:,1));
    Inds4 = 4:4:numel(DICE(:,1));

    %% Display results

    disp('Volume 1')
    
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICE(Inds1,1))),'(',num2str(std(DICE(Inds1,1))),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds1,1) - PredictedCentroid(Inds1,1))*pixelpermm)),'(',num2str(std((RealCentroid(Inds1,1) - PredictedCentroid(Inds1,1))*pixelpermm)),')'])
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds1,1) - PredictedCentroid(Inds1,1))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds1,1) - PredictedCentroid(Inds1,1))*pixelpermm,95)),'mm'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds1,2) - PredictedCentroid(Inds1,2))*pixelpermm)),'(',num2str(std((RealCentroid(Inds1,2) - PredictedCentroid(Inds1,2))*pixelpermm)),')'])    
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds1,2) - PredictedCentroid(Inds1,2))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds1,2) - PredictedCentroid(Inds1,2))*pixelpermm,95)),'mm'])    
    disp(['total distance [mean(SD)] = ',num2str(mean(Dist(Inds1,1))),'(',num2str(std(Dist(Inds1,1))),')'])
    disp(['MSD [mean(SD)] = ',num2str(mean(MSDNorm(Inds1,1))*pixelpermm),'(',num2str(std(MSDNorm(Inds1,1))*pixelpermm),')'])  
    disp(['Max distance = ',num2str(max(Dist(Inds1,1)))])

    disp('UnWarped')
    
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICEUnwarped(Inds1,1))),'(',num2str(std(DICEUnwarped(Inds1,1))),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds1,1) - UnwarpedCentroid(Inds1,1))*pixelpermm)),'(',num2str(std((RealCentroid(Inds1,1) - UnwarpedCentroid(Inds1,1))*pixelpermm)),')'])
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds1,1) - UnwarpedCentroid(Inds1,1))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds1,1) - UnwarpedCentroid(Inds1,1))*pixelpermm,95)),'mm'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds1,2) - UnwarpedCentroid(Inds1,2))*pixelpermm)),'(',num2str(std((RealCentroid(Inds1,2) - UnwarpedCentroid(Inds1,2))*pixelpermm)),')'])    
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds1,2) - UnwarpedCentroid(Inds1,2))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds1,2) - UnwarpedCentroid(Inds1,2))*pixelpermm,95)),'mm'])    
    %disp(['total distance [mean(SD)] = ',num2str(mean(Dist2(:,1))),'(',num2str(std(Dist2(:,1))),')'])    
    disp(['MSD [mean(SD)] = ',num2str(mean(MSDUnWarped(Inds1,1)*pixelpermm)),'(',num2str(std(MSDUnWarped(Inds1,1)*pixelpermm)),')'])   
    %disp(['Max distance = ',num2str(max(Dist2(:,1)))])

    disp('Volume 2')
    
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICE(Inds2,1))),'(',num2str(std(DICE(Inds2,1))),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds2,1) - PredictedCentroid(Inds2,1))*pixelpermm)),'(',num2str(std((RealCentroid(Inds2,1) - PredictedCentroid(Inds2,1))*pixelpermm)),')'])
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds2,1) - PredictedCentroid(Inds2,1))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds2,1) - PredictedCentroid(Inds2,1))*pixelpermm,95)),'mm'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds2,2) - PredictedCentroid(Inds2,2))*pixelpermm)),'(',num2str(std((RealCentroid(Inds2,2) - PredictedCentroid(Inds2,2))*pixelpermm)),')'])    
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds2,2) - PredictedCentroid(Inds2,2))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds2,2) - PredictedCentroid(Inds2,2))*pixelpermm,95)),'mm'])    
    disp(['total distance [mean(SD)] = ',num2str(mean(Dist(Inds2,1))),'(',num2str(std(Dist(Inds2,1))),')'])
    disp(['MSD [mean(SD)] = ',num2str(mean(MSDNorm(Inds2,1))*pixelpermm),'(',num2str(std(MSDNorm(Inds2,1))*pixelpermm),')'])  
    disp(['Max distance = ',num2str(max(Dist(Inds2,1)))])

    disp('UnWarped')
    
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICEUnwarped(Inds2,1))),'(',num2str(std(DICEUnwarped(Inds2,1))),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds2,1) - UnwarpedCentroid(Inds2,1))*pixelpermm)),'(',num2str(std((RealCentroid(Inds2,1) - UnwarpedCentroid(Inds2,1))*pixelpermm)),')'])
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds2,1) - UnwarpedCentroid(Inds2,1))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds2,1) - UnwarpedCentroid(Inds2,1))*pixelpermm,95)),'mm'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds2,2) - UnwarpedCentroid(Inds2,2))*pixelpermm)),'(',num2str(std((RealCentroid(Inds2,2) - UnwarpedCentroid(Inds2,2))*pixelpermm)),')'])    
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds2,2) - UnwarpedCentroid(Inds2,2))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds2,2) - UnwarpedCentroid(Inds2,2))*pixelpermm,95)),'mm'])    
    %disp(['total distance [mean(SD)] = ',num2str(mean(Dist2(:,1))),'(',num2str(std(Dist2(:,1))),')'])    
    disp(['MSD [mean(SD)] = ',num2str(mean(MSDUnWarped(Inds2,1)*pixelpermm)),'(',num2str(std(MSDUnWarped(Inds2,1)*pixelpermm)),')'])   
    %disp(['Max distance = ',num2str(max(Dist2(:,1)))])

    disp('Volume 3')
    
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICE(Inds3,1))),'(',num2str(std(DICE(Inds3,1))),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds3,1) - PredictedCentroid(Inds3,1))*pixelpermm)),'(',num2str(std((RealCentroid(Inds3,1) - PredictedCentroid(Inds3,1))*pixelpermm)),')'])
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds3,1) - PredictedCentroid(Inds3,1))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds3,1) - PredictedCentroid(Inds3,1))*pixelpermm,95)),'mm'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds3,2) - PredictedCentroid(Inds3,2))*pixelpermm)),'(',num2str(std((RealCentroid(Inds3,2) - PredictedCentroid(Inds3,2))*pixelpermm)),')'])    
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds3,2) - PredictedCentroid(Inds3,2))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds3,2) - PredictedCentroid(Inds3,2))*pixelpermm,95)),'mm'])    
    disp(['total distance [mean(SD)] = ',num2str(mean(Dist(Inds3,1))),'(',num2str(std(Dist(Inds3,1))),')'])
    disp(['MSD [mean(SD)] = ',num2str(mean(MSDNorm(Inds3,1))*pixelpermm),'(',num2str(std(MSDNorm(Inds3,1))*pixelpermm),')'])  
    disp(['Max distance = ',num2str(max(Dist(Inds3,1)))])

    disp('UnWarped')
    
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICEUnwarped(Inds3,1))),'(',num2str(std(DICEUnwarped(Inds3,1))),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds3,1) - UnwarpedCentroid(Inds3,1))*pixelpermm)),'(',num2str(std((RealCentroid(Inds3,1) - UnwarpedCentroid(Inds3,1))*pixelpermm)),')'])
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds3,1) - UnwarpedCentroid(Inds3,1))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds3,1) - UnwarpedCentroid(Inds3,1))*pixelpermm,95)),'mm'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds3,2) - UnwarpedCentroid(Inds3,2))*pixelpermm)),'(',num2str(std((RealCentroid(Inds3,2) - UnwarpedCentroid(Inds3,2))*pixelpermm)),')'])    
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds3,2) - UnwarpedCentroid(Inds3,2))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds3,2) - UnwarpedCentroid(Inds3,2))*pixelpermm,95)),'mm'])    
    %disp(['total distance [mean(SD)] = ',num2str(mean(Dist2(:,1))),'(',num2str(std(Dist2(:,1))),')'])    
    disp(['MSD [mean(SD)] = ',num2str(mean(MSDUnWarped(Inds3,1)*pixelpermm)),'(',num2str(std(MSDUnWarped(Inds3,1)*pixelpermm)),')'])   
    %disp(['Max distance = ',num2str(max(Dist2(:,1)))])

    disp('Volume 4')
    
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICE(Inds4,1))),'(',num2str(std(DICE(Inds4,1))),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds4,1) - PredictedCentroid(Inds4,1))*pixelpermm)),'(',num2str(std((RealCentroid(Inds4,1) - PredictedCentroid(Inds4,1))*pixelpermm)),')'])
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds4,1) - PredictedCentroid(Inds4,1))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds4,1) - PredictedCentroid(Inds4,1))*pixelpermm,95)),'mm'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds4,2) - PredictedCentroid(Inds4,2))*pixelpermm)),'(',num2str(std((RealCentroid(Inds4,2) - PredictedCentroid(Inds4,2))*pixelpermm)),')'])    
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds4,2) - PredictedCentroid(Inds4,2))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds4,2) - PredictedCentroid(Inds4,2))*pixelpermm,95)),'mm'])    
    disp(['total distance [mean(SD)] = ',num2str(mean(Dist(Inds4,1))),'(',num2str(std(Dist(Inds4,1))),')'])
    disp(['MSD [mean(SD)] = ',num2str(mean(MSDNorm(Inds4,1))*pixelpermm),'(',num2str(std(MSDNorm(Inds4,1))*pixelpermm),')'])  
    disp(['Max distance = ',num2str(max(Dist(Inds4,1)))])

    disp('UnWarped')
    
    disp(['Dice Coefficient [mean(SD)] = ',num2str(mean(DICEUnwarped(Inds4,1))),'(',num2str(std(DICEUnwarped(Inds4,1))),')'])
    disp(['x distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds4,1) - UnwarpedCentroid(Inds4,1))*pixelpermm)),'(',num2str(std((RealCentroid(Inds4,1) - UnwarpedCentroid(Inds4,1))*pixelpermm)),')'])
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds4,1) - UnwarpedCentroid(Inds4,1))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds4,1) - UnwarpedCentroid(Inds4,1))*pixelpermm,95)),'mm'])
    disp(['y distance(mm) [mean(SD)] = ',num2str(mean((RealCentroid(Inds4,2) - UnwarpedCentroid(Inds4,2))*pixelpermm)),'(',num2str(std((RealCentroid(Inds4,2) - UnwarpedCentroid(Inds4,2))*pixelpermm)),')'])    
    disp(['5% & 95% error = ',num2str(prctile((RealCentroid(Inds4,2) - UnwarpedCentroid(Inds4,2))*pixelpermm,5)),'mm ',num2str(prctile((RealCentroid(Inds4,2) - UnwarpedCentroid(Inds4,2))*pixelpermm,95)),'mm'])    
    %disp(['total distance [mean(SD)] = ',num2str(mean(Dist2(:,1))),'(',num2str(std(Dist2(:,1))),')'])    
    disp(['MSD [mean(SD)] = ',num2str(mean(MSDUnWarped(Inds4,1)*pixelpermm)),'(',num2str(std(MSDUnWarped(Inds4,1)*pixelpermm)),')'])   
    %disp(['Max distance = ',num2str(max(Dist2(:,1)))])


    
    PatientLabels = {''};
    for i = 1:numel(Files)
        if i <= 5
            PatientLabels{numel(PatientLabels)+1} = ['Oropharynx Patient',num2str(i)];
        elseif i > 5 && i <=10
            PatientLabels{numel(PatientLabels)+1} = ['Glottis Patient',num2str(i)];
        else
            PatientLabels{numel(PatientLabels)+1} = ['Nasopharynx',num2str(i)];
        end
    end
    PatientLabels{numel(PatientLabels)+1} = 'Total';
    PatientLabels{numel(PatientLabels)+1} = '';      
    
    TumourLocation = cell(size(PredictedCentroid(:,3)));
    
    for i = 1:numel(TumourLocation)
        if PredictedCentroid(i,3) < 5
            TumourLocation{i} = 'Oropharynx';
        elseif PredictedCentroid(i,3) >= 5 && PredictedCentroid(i,3) < 10
            TumourLocation{i} = 'Glottis';
        else
            TumourLocation{i} = 'Nasopharynx';
        end
    end
    
%     T = table(XDist(1:numel(RealCentroid(:,1)),1),YDist(1:numel(RealCentroid(:,1)),1),DICEXCorr(:,1),MSDXCorr(:,1),...
%         XDistunwarped(1:numel(RealCentroid(:,1)),1),YDistunwarped(1:numel(RealCentroid(:,1)),1),DICEUnwarped(:,1),MSDUnWarped(:,1),...
%         PredictedCentroidXCorr(:,3),TumourLocation,'VariableNames',varNames);

    %% Save results

    XDist = [(RealCentroid(:,1) - PredictedCentroid(:,1))*pixelpermm (RealCentroid(:,3) - PredictedCentroid(:,3))*pixelpermm];
    YDist = [(RealCentroid(:,2) - PredictedCentroid(:,2))*pixelpermm (RealCentroid(:,4) - PredictedCentroid(:,4))*pixelpermm];
    
    XDistUnWarped = [(RealCentroid(:,1) - UnwarpedCentroid(:,1))*pixelpermm (RealCentroid(:,3) - UnwarpedCentroid(:,3))*pixelpermm];
    YDistUnWarped = [(RealCentroid(:,2) - UnwarpedCentroid(:,2))*pixelpermm (RealCentroid(:,4) - UnwarpedCentroid(:,4))*pixelpermm];
    
    varNames = {'XDistMultiVol',...
        'YDistMultiVol',...
        'DICEMultiVol',...
        'MSDMultiVol',...
        'XDistUnWarpedMultiVol',...
        'YDistUnWarpedMultiVol',... 
        'DICEUnwarpedMultiVol',...
        'MSDUnwarpedMultiVol',......
        'PatientNum','TumourLocation','VolumeFlag'};      
    
    T = table(XDist(:,1),...
        YDist(:,1),...
        DICE(:,1),...
        MSDNorm(:,1)*pixelpermm,...
        XDistUnWarped(:,1),...
        YDistUnWarped(:,1),...
        DICEUnwarped(:,1),...
        MSDUnWarped(:,1)*pixelpermm,...
        PredictedCentroid(:,end),TumourLocation,VolumeFlag(:,1),...
        'VariableNames',varNames);

    writetable(T,OutputFile)    
    
end