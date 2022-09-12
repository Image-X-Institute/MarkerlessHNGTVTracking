function[MSDVal] = MSD(Im1,Im2)
%% function[MSDVal] = MSD(Im1,Im2)
%Calculates the Mean Surface Distance which is the average value of the
%minimum distance between the two surfaces
%Mark Gardner 11/05/2021


%Im1 & Im2 - binarised reference image
%MSDVal - Mean Surface Distance value

    %% Convert Images to outlines
  
    ImOutline1 = bwboundaries(Im1);
    ImOutline2 = bwboundaries(Im2);
    
    if (numel(ImOutline1) == 0) || (numel(ImOutline2) == 0)
       error('Outline not able to be detected in input image(s)') 
    end
    
    if numel(ImOutline1) > 1
        boundaryNums = zeros(1,numel(numel(ImOutline1)));
        for i = 1:numel(ImOutline1)
            boundaryNums(i) = numel(ImOutline1{i}(:,1));
        end
        [~,Inds1] = max(boundaryNums);
    else
        Inds1 = 1;
    end
    
    if numel(ImOutline2) > 1
        boundaryNums = zeros(1,numel(numel(ImOutline2)));
        for i = 1:numel(ImOutline2)
            boundaryNums(i) = numel(ImOutline2{i}(:,1));
        end
        [~,Inds2] = max(boundaryNums);        
    else
        Inds2 = 1;
    end    
    
    %% Calculate MSD    
    dist12 = zeros(size(ImOutline1{Inds1}(:,1)));
    for i = 1:numel(ImOutline1{Inds1}(:,1))
        distLoop = sqrt((ImOutline1{Inds1}(i,1) - ImOutline2{Inds2}(:,1)).^2 + (ImOutline1{Inds1}(i,2) - ImOutline2{Inds2}(:,2)).^2);
        dist12(i) = min(distLoop);
    end
    
    dist21 = zeros(size(ImOutline2{Inds2}(:,1)));
    for i = 1:numel(ImOutline2{Inds2}(:,1))
        distLoop = sqrt((ImOutline2{Inds2}(i,1) - ImOutline1{Inds1}(:,1)).^2 + (ImOutline2{Inds2}(i,2) - ImOutline1{Inds1}(:,2)).^2);
        dist21(i) = min(distLoop);
    end    
    
    MSDVal = (sum(dist12) + sum(dist21))/(numel(dist12) + numel(dist21));
    
end