% FILEMACPC determines the correct path depending on mac or pc
% installation.
%
% filemacpc gives the correct filepath for a given file or folder. it will
% correct the back(front) slash conventions for the operating
% system. it will also remove the trailing slash, if provided in the input.
%
% example mac output
% /Users/rpjb/Dropbox/project - geniculate originaldata/large fov
% example pc outptu
% 'C:\Users\rpjb\My Dropbox\project - geniculate originaldata\large fov'
%
% type: function
%
% inputs: 
%   folder:  string identifying folder to be called
%
% outputs:
%   output:  string identifying folder, correctly formatted for machine
%
% Robert Barretto, robertb@gmail.com
% 3/4/13 added dropbox autocorrect functionality

function [output] = Filemacpc(folder)
    % for debugging purposes
    %folder = 'C:/Users/rpjb/Dropbox/project - hippocampus/imaging';
    %folder = '\Users\rpjb\Dropbox\';
    
    % set system base folder depending on machine type
    % deconstruct folder into subfolders
    % the first case deals with USB drives that have other letter volumes
    if ispc && strcmp(folder(2:3),':\')
        base = folder(1:3);
    elseif ispc && ~strcmp(folder(2:3),':\')
        base = 'C:\';
    else
        base = '/';
    end
    
    % checks the functionality of the backslash/frontslash.
    if ispc
        folder = strrep(folder,'/','\');
    else
        folder = strrep(folder,'\','/'); 
    end
    % strip the last file separator by convention, if it was presented
    if strcmpi(folder(end),filesep)
        folder(end) = [];
    end               
    
    % remove string if it has a colon
    [firstfolder, remain] = strtok(folder,filesep); 
    if ~isempty(strfind(firstfolder,':'))
        folder = remain;
    end
    
    % attach the proper base
    output = [base,folder(2:end)];
    % if output folder doesn't exist check if there is a dropbox issue
    if ~exist(output)
        if findstr('Dropbox',output)
            if findstr('My Dropbox',output)
                output = strrep(output,'My Dropbox','Dropbox');
            else
                output = strrep(output,'Dropbox','My Dropbox');
            end
        end
    end
    
     % create output, and present warnings if either base or folder are not
     % found
    if ~isdir(base)
        disp(['Warning: base folder, ' base ' was not found.'])
        output = [];
    else
        if exist(output) == 0
            disp(['Warning: ' output ' was not found.'])
            output = [];
        end
    end    
end
