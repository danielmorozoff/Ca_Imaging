% FINDFILES Searches a folder (and its subfolders) for files matching a
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
%   outputs:  structure array containing file information
%
% Robert Barretto, robertb@gmail.com
% 05/03/2011 12:29pm

function [outputs] = FindFiles(dirpath, search, depth )
    %% check inputs
    % ensure filesep is missing from end of dirpath
    if strcmp(dirpath(length(dirpath)),filesep)
        dirpath = dirpath(1:-1+length(dirpath));
    end
  
    %% get the files at depth 1
    outputs = [];
    files = getsubfile(dirpath);
    for j=1:length(files)
        if validatefile(files(j).name,search)
            outputs = [outputs; files(j)];
        end
    end
    
    
    %% get all the subfolders
    if depth>1
        % at depth 2
        [subfolders] = getsubfolder(dirpath);
            if depth>2
                startindex = 1;
                for j=3:depth
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
        %% get the files of all subfolders         
    
        for i=1:length(subfolders)
%             disp([subfolders(i).path, filesep, subfolders(i).name])
            files = getsubfile([subfolders(i).path, filesep, subfolders(i).name]);
            for j=1:length(files)
                if validatefile(files(j).name,search)
                    outputs = [outputs; files(j)];
                end
            end
        end        
    end
    
        
    %% subfunctions
    % scans a given folder for files
    function [out] = getsubfile(folder)
        foldercontents = dir(folder);
        % skip the first two folders because of '.' and '..'
        foldercontents = foldercontents([3:length(foldercontents)]);
        out = foldercontents([foldercontents.isdir] == 0);
        % appends absolute path 
        for i=1:length(out)
            out(i).path = folder;
        end
    end
        
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
    % validates file against a search string
    function [output] = validatefile(file,search)
        result = regexp(file,search);
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
    
    

