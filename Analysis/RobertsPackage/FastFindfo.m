% FASTFINDFO Performs skewness analysis to determine baseline of trace.
%
% fastfindfo receives a time-series vector and performs various operations
% on the distribution of vector values. it uses this information to
% identify the largest subvector with a near symmetric normal distribution,
% using skewness as a measure.
% 
% transients are identified based on the following assumptions. in
% legitmate calcium imaging data, transients are assumed to only be
% increasing intensity events. thus, a firing neuron ought to hava a
% positive skewness, relative to the underlying distribution of intensities
% when the cell is silent. 
%
% type: function
% 
% inputs:
%   in: 1-D array containing calcium trace
%
% outputs:
%   fo: scalar representing the baseline fluorescence level
%
% dependencies on custom functions:
%   none
%
% Robert Barretto, robertb@gmail.com
% 03/21/13 06:22am

function [fo] = FastFindfo(in)

% initialize variables
numpoints = length(in);
% sort data into rising area
sortinput = sort(in);

% calculate largest subset of data with near zero skewness
bigskew = skewness(sortinput);
if bigskew >= 0
    i = numpoints-1;
    hit = 0;
    while i>1 && hit == 0
        littleskew = skewness(sortinput(1:i));
        if bigskew > 0 && littleskew < 0
            hit = 1;
        end
        bigskew = littleskew;
        i = i-1;        
    end
    fo = median(sortinput(1:i));
else
    i = numpoints;
    fo = median(sortinput(1:i));
end

