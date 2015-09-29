% STACK2MONTAGE Flattens a stack into sequential images
%
% type: function
%   
% inputs:
%   originaldata: 3d volume stack
%    
% outputs:
%   none
%   creates a figure display
%
% dependencies on custom functions:
%   none
%
% Robert Barretto, robertb@gmail.com
% 03/24/13 9:24pm initial commit

function [] = stack2montage(originaldata)
datadim = size(originaldata);
figure;
data = reshape(originaldata,datadim(1),datadim(2),1,datadim(3));
montage(data,'DisplayRange',[])
end