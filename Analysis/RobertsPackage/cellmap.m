% CELLMAP Identifies a box subregion of an image when given box parameters.
%
% calciumtracefilter receives x,y and the radii of a rectangular box, and
% returns the data within the box.
%
% type: function
%
% inputs:
%   xpt: the horizontal center coordinate 0 at left
%   ypt: the vertical center coordinate 0 at top
%   xr: horizontal radius
%   yr: vertical radius
%   imgin: a 2D or 3D image
%
% outputs:
%   imgout: calculated subset of imgin
%   xrange: coordinates used
%   yrange: coordinates used
% dependencies on custom functions:
%   none
%
% Robert Barretto, robertb@gmail.com
% 03/24/13 7:29am initial commit

function [imgout,xrange,yrange] = cellmap(xpt,ypt,xr,yr,imgin )
    % calculate x,y coordinate ranges
    imgindim = size(imgin);
    xrange = max(xpt-xr,1):min(xpt+xr,imgindim(2));
    yrange = max(ypt-yr,1):min(ypt+yr,imgindim(1));    
    % process image or imagestack
    if length(imgindim) == 2
        imgout = imgin(yrange,xrange);
    elseif length(imgindim) == 3
        imgout = imgin(yrange,xrange,:);
    end    
end

