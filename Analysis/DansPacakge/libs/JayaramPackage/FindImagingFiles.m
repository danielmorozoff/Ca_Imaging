% CalciumTraceFilter Performs basic filtering on a 1-D calcium imaging
% trace.
%
% calciumtracefilter receivea a time-series vector and performs filtering
% (either 'median' or 'integrated').
%
% type: function
%
% inputs:
%   folder: string of absolute folder path containing data
%    
% outputs:
%   outputmeanfile: filename of averageimage
%   outputfiles:    cell array of filenames for each trial
%   celldatafile:   filename of celldata.mat
%   celltracesfile: filename of celltraces.mat
%   isRobert:       binary indicator if dataset from Robert or Jayaram
%
% dependencies on custom functions:
%   FindFiles
%
% Robert Barretto, robertb@gmail.com
% 03/24/13 7:50pm initial commit
% 03/27/13 4:10pm add tag to identify type of dataset

function [outputmeanfile, outputfiles, celldatafile, celltracesfile, celltransientsfile, isRobert] = FindImagingFiles(folder)
    % check if folder exists
    if ~exist(folder)
        disp('Folder does not exist')
        outputfiles = '';
        outputmeanfile = '';
        celldatafile = '';
        celltracesfile = '';
        isRobert = '';
        celltransientsfile = '';
        return
    end
    % see what kind of files there are
    tic
    s = FindFiles(folder,'Image_Registration',3);
    toc
    r = FindFiles(folder,'aligned',3);
    toc
    if isempty(r) && isempty(s)
        disp('Aligned files not found')
        outputfiles = '';
        outputmeanfile = '';
        celldatafile = '';
        celltracesfile = '';
        celltransientsfile = '';
        isRobert = '';
        return
    end

    if ~isempty(s) % case of Jayaram's dataset
        isRobert = 0;
        for i=1:length(s)
            outputfiles{i} = fullfile(s(i).path,s(i).name);
        end
        % identify given avgimage
        ss = FindFiles(folder,'\S+sdmean\S+',3);
    end
    if ~isempty(r); % case of Robert's dataset
        isRobert = 1;
        for i=1:length(r)
            outputfiles{i} = fullfile(r(i).path,r(i).name);
        end
        % identify given avgimage    
        ss = FindFiles(folder,'\S+avgimg\S+',3);
        if isempty(ss)
            ss = FindFiles(folder,'AVG_\S+',3);
        end
    end
    if ~isempty(ss)
        outputmeanfile = fullfile(ss(end).path,ss(end).name);    
    else
        outputmeanfile = '';
        disp('SD Files not found');
    end 
    
    % identify celldata and celltraces files   
    celldatafile = fullfile(folder,'celldata.mat');
    celltracesfile = fullfile(folder,'celltraces.mat');
    celltransientsfile = fullfile(folder,'celltransients.mat');
    if ~exist(celldatafile)
        celldatafile = '';
    end
    if ~exist(celltracesfile)
        celltracesfile = '';
    end
    if ~exist(celltransientsfile)
        celltransientsfile = '';
    end
end
