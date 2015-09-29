% CalciumTraceFilter Performs basic filtering on a 1-D calcium imaging
% trace.
%
% calciumtracefilter receivea a time-series vector and performs filtering
% (either 'median' or 'integrated').
%
% type: function
%
% inputs:
%   signal:  1-D vector containing calcium trace
%   filtertype: string specifying filter either 'median' or 'integrated'
%   window: integer specifying size of filter window
%    
% outputs:
%   output:  1-D array containing filtered calcium trace
%   filtertype: passthrough string specifying filter either 'median' or 'integrated'
%
% dependencies on custom functions:
%   none
%
% Robert Barretto, robertb@gmail.com
% 05/17/11 06:31pm
% 03/23/13 11:19am fix edge case scenario with median filter

function [output, filtertype] = CalciumTraceFilter(signal,filtertype,window)
    if strcmp(filtertype,'median')
        % performs median filter over a sliding window       
        % pad the start and ends to avoid edge issues
        numpad = ceil(window / 2);        
        longsignal = [signal(1)*ones(numpad,1); signal; signal(end)*ones(numpad,1)];        
        output = medfilt1(longsignal,window);
        output = output((numpad+1):(end-numpad));        
    elseif strcmp(filtertype,'integrated')
        % performs integrated response filter over a sliding window
        output = filter(ones(1,window),1,signal);
    else
        % no operation performed
        output = signal;
    end
end
