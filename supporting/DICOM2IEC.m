function DICOM2IEC(dicomList,outputFile,waterAtt)
%% DICOM2IEC(dicomList,outputFile,waterAtt)
% ------------------------------------------
% FILE   : DICOM2IEC.m
% AUTHOR : Andy Shieh, ACRF Image X Institute, The University of Sydney
% DATE   : 2018-02-14  Created.
% ------------------------------------------
% PURPOSE
%   Convert a series of DICOM files to a .mha file in IEC geometry.
% ------------------------------------------
% INPUT
%   dicomList:      A list of dicom files in cell format (use "lscell")
%   outputFile:     Output file path (.mha format)
%   waterAtt:       (Optional) The intended attenuation value for water.
%                   Default: 0.013
%                   If left as "NaN", then no transformation will be
%                   applied.

%% Input check
if nargin < 2
    error('ERROR: dicomList and outputFile must be input.');
end

if nargin < 3
    waterAtt = 0.013;
end

%% Read DICOM files
for k = 1:length(dicomList)
    dcmHeader{k} = dicominfo(dicomList{k});
    M(:,:,k) = dicomread(dicomList{k});
    try
        sl(k) = dcmHeader{k}.SliceLocation;
    catch
        sl(k) = dcmHeader{k}.ImagePositionPatient(3);
    end
end

%% Convert to IEC geometry and scale intensity value
M = permute(single(M),[1 3 2]);
signTransformMat = dcmHeader{1}.ImageOrientationPatient;
signTransformMat = signTransformMat(signTransformMat~=0);

% Inverting LR and AP if needed
if dcmHeader{1}.ImageOrientationPatient(1) ~= 0
    if signTransformMat(1) < 0
        M = M(:,:,end:-1:1);
    end
    if signTransformMat(2) > 0
        M = M(end:-1:1,:,:);
    end
else
    if signTransformMat(1) < 0
        M = M(end:-1:1,:,:);
    end
    if signTransformMat(2) > 0
        M = M(:,:,end:-1:1);
    end
end

% Axial slice ordering
[~,indSlice] = sort(sl);
[~,indSlice] = unique(sl);
M = M(:,indSlice(end:-1:1),:);

% Swap LR and AP if needed
if dcmHeader{1}.ImageOrientationPatient(1) ~= 0
    M = permute(M,[3 2 1]);
end

if ~isnan(waterAtt)
    if min(M(:)) >= 0
        M = M * waterAtt / 1000;
    else
        M = (M + 1000) * waterAtt / 1000;
    end
end

%% Construct header
load(which('mhaHeaderTemplate.mat'));
mhaHeader = mhaHeaderTemplate;
mhaHeader.Dimensions = size(M);
%Thickness = dcmHeader{1}.SliceThickness;
%if isempty(Thickness)
try
    Thickness = abs(dcmHeader{2}.SliceLocation-dcmHeader{1}.SliceLocation);
catch
    Thickness = dcmHeader{1}.SliceThickness;
end
%end
mhaHeader.PixelDimensions = [dcmHeader{1}.PixelSpacing(1),Thickness,dcmHeader{1}.PixelSpacing(2)];
mhaHeader.Offset = -0.5 * (mhaHeader.Dimensions - 1) .* mhaHeader.PixelDimensions;

%% Write to file
MhaWrite(mhaHeader,M,outputFile);

end