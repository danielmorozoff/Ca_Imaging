% NULLIMAGEBORDERS Sets to zero various edges of an image.
%
% nullimageborders receives the borders dimensions and pads a 2D/3D image.
%
% type: function
%
% inputs:
%   left: scalar nulling leftmost columns
%   right: scalar nulling rightmost columns
%   top: scalar nulling topmost rows
%   bottom: scalar nulling bottommost rows
%   imgin: a 2D or 3D image
%
% outputs:
%   imgout: calculated subset of imgin
%
% dependencies on custom functions:
%   none
%
% Robert Barretto, robertb@gmail.com
% 03/26/13 1:18pm initial commit

function [ imgout ] = NullImageBorders(left,right,top,bottom,imgin)
    imgindim = size(imgin);
    imgout = imgin;
   
    % process image or imagestack
    if length(imgindim) == 2
        if left>0
            imgout(:,1:left) = 0;
        end
        if right>0
            imgout(:,end+1-(1:right)) = 0;
        end
        if top>0
            imgout(1:top,:) = 0;
        end
        if bottom>0
            imgout(end+1-(1:bottom),:) = 0;
        end        
    elseif length(imgindim) == 3
        if left>0
            imgout(:,1:left,:) = 0;
        end
        if right>0
            imgout(:,end+1-(1:right),:) = 0;
        end
        if top>0
            imgout(1:top,:,:) = 0;
        end
        if bottom>0
            imgout(end+1-(1:bottom),:,:) = 0;
        end        
    end    
end
