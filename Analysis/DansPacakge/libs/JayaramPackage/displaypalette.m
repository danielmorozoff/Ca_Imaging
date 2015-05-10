% displaypalette Visualizes a colormap.
%
% displaypalette calculates and displays an RGB colormap
%
% type: script
%
% inputs:
%    
% outputs:
%
% dependencies on custom functions:
%   none
%
% Robert Barretto, robertb@gmail.com
% 03/24/13 11:19am fix edge case scenario with median filter

% create a basic colormap
numbits = 256;
cmap = gray(numbits);
cmap(:,1) = 1;
cmap = flipud(cmap);
% create an image
cimg = zeros([numbits numbits 3]);

for i = 1:numbits % for each row
    for j=1:3 % for each color
        cimg(i,:,j) = cmap(i,j);
    end
end

% display image
figure
imshow(cimg)