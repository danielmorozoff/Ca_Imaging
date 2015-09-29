% IMREADTIFFSTACK Opens a single file tiff stack.
%
% type: function
%
% inputs:
%   file: string of absolute filepath
%    
% outputs:
%   data: 3d image array
%
% dependencies on custom functions:
% 
% Robert Barretto, robertb@gmail.com
% 03/24/13 7:50pm initial commit

function [data] = imreadtiffstack(file)
    infoimg = imfinfo(file);
    w = infoimg(1).Width;
    h = infoimg(1).Height;
    numframes = length(infoimg);
    data = zeros(h,w,numframes,'uint16');
    for i=1:numframes
        data(:,:,i) = imread(file,i);
    end
end

