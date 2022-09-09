function [info M] = MhaRead(filename)
%% [info M] = MhaRead(filename)
% ------------------------------------------
% FILE   : MhaRead.m
% AUTHOR : Andy Shieh, School of Physics, The University of Sydney
% DATE   : 2012-11-08  Created.
% ------------------------------------------
% PURPOSE
%   Read metaimage mha header and image.
%   If the file is a archived file in the format of .7z, .zip, or .rar, the
%   program will unzip and read.
% ------------------------------------------
% INPUT
%   filename : The full path to the file.
% ------------------------------------------
% OUTPUT
%   info  : header information in a struct.
%   M     : the image stored in a 3D matrix.
% ------------------------------------------

%% Checking input arguments & Opening the file

M = [];

% If no input filename => Open file-open-dialog
if nargin < 1
    
    % go into a default directory
    DefaultDir = pwd;
    
    % get the input file (hnc) & extract path & base name
    [FileName,PathName] = uigetfile( ...
        {'*.mha;*.7z;*.zip;*.rar','MetaImage (*.mha) or Archive file (*.7z, *.zip, *.rar)';...
        '*.*','All files (*.*)'},...
        'Select an image or image archive file', ...
        DefaultDir);
    
    % catch error if no file selected
    if isnumeric(FileName)
        info = struct([]);
        error('ERROR: No file selected.');
    end
    
    % make same format as input
    filename = fullfile(PathName, FileName);
    
end
filename = strtrim(filename);

%% Unzip if needed
[~,~,ext] = fileparts(filename);
if strcmpi(ext,'.7z') || strcmpi(ext,'.zip') || strcmpi(ext,'.rar')
    tempdir = tempname;
    mkdir(tempdir);
    [~,~] = system(['7z x "',filename,'" -o"',tempdir,'"']);
    fl = lscell(fullfile(tempdir,'*.mha'));
    filename = fl{1};
end

%% Reading the header and body using external module
info = mha_read_header(filename);

if nargout == 1
    if exist('tempdir','var')
        rmdir(tempdir,'s');
    end
    return;
else
    M = mha_read_volume(info);
    if exist('tempdir','var')
        rmdir(tempdir,'s');
    end
end

return
