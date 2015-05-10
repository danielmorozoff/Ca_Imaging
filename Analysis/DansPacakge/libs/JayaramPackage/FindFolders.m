% FINDFOLDERS Searches a folder (and its subfolders) for folders matching a
% search string.
%
% type: function
% 
% inputs:
%   dirpath:  string identifying main directory path
%   search:  string regular expression identifying the type of folders to look for
%   depth:  integer representing the level of subfolders to examine.  1
%     represents only searching the main folder.
%
% outputs:
%   outputs:  structure array containing directory information
%
% Robert Barretto, robertb@gmail.com
% 05/03/2011 12:29pm

function [outputs] = findfolders(dirpath, search, depth )
    %% check inputs
    % ensure filesep is missing from end of dirpath
    if strcmp(dirpath(length(dirpath)),filesep)
        dirpath = dirpath(1:-1+length(dirpath));
    end
  
    %% get all the subfolders
    [subfolders] = getsubfolder(dirpath);
    if depth>1
        startindex = 1;
        for j=2:depth
            jsubfolders = [];
            for m=startindex:length(subfolders)
                temp = getsubfolder([subfolders(m).path filesep subfolders(m).name] );
                if isempty(temp) == 0
                    jsubfolders = [jsubfolders ; temp];
                end
            end
            startindex = length(subfolders) + 1;
            subfolders = [subfolders; jsubfolders];
        end
    end
    
    %% validate the subfolders
    outputs = [];
    for i=1:length(subfolders)
        if validatefolder(subfolders(i).name,search)
            outputs = [outputs; subfolders(i)];
        end
    end
        
        
    %% subfunctions
    % scans a given folder for subfolders
    function [out] = getsubfolder(folder)    
        foldercontents = dir(folder);
        % skip the first two folders because of '.' and '..'
        foldercontents = foldercontents([3:length(foldercontents)]);
        out = foldercontents([foldercontents.isdir] == 1);
        % appends absolute path 
        for i=1:length(out)
            out(i).path = folder;
        end
    end
    % validates folder against a search string
    function [output] = validatefolder(folder,search)
        result = regexp(folder,search);
        if isempty(result)
            output = 0;
        elseif result ~=1
            output = 0;
        elseif result == 1
            output = 1;
        else 
            output = 0;
        end
    end
end    
    
    

