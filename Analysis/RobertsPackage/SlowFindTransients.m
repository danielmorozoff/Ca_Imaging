% SLOWFINDTRANSIENTS Performs transient detection on a 1-D calcium imaging trace.
%
% slowfindtransients receives a time-series vector and performs various operations
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
% this is the slow version of FASTFINDTRANSIENTS. it allows one to look at
% the underlying distributions and performs a cumuluative skewness.
%
% type: function
% 
% inputs:
%   in: 1-D array containing calcium trace
%   xtrace: 1-D array containing timestamps of each datapoint of in
%
% outputs:
%   out: stucture for identified transients
%       out.start:  first index of transient
%       out.end:    last index of transient
%       out.x:      1-D vector of timepoints
%       out.y:      1-D vector of amplitude
%       out.sn:     peak amplitude
%   sigma: dispersion of background distribution
%   fo: scalar representing the baseline fluorescence level
%
% dependencies on custom functions:
%   none
%
% Robert Barretto, robertb@gmail.com
% 03/21/13 06:22am

function [out,sigma,fo] = SlowFindTransients(in,xtrace)

% initialize variables
numpoints = length(in);
% sort data into rising area
sortinput = sort(in);

%% for graphical purposes
% initialize variables
kurtvalues = zeros(1,numpoints);
skewvalues = kurtvalues;
medvalues = kurtvalues;
% sort data into rising area
for i=1:numpoints
    kurtvalues(i) = kurtosis(sortinput(1:i));
    skewvalues(i) = skewness(sortinput(1:i));
    medvalues(i) = median(sortinput(1:i));
end

% calculate largest subset of data with near zero skewness
bigskew = skewness(sortinput);

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

% calculate standard deviation of the underlying distribution with a 2
% sigma cutoff
sigma = 1.4826*mad(sortinput(1:i),1);
cutoff = fo+2*sigma;

% find putative transients
eventvector = in >= cutoff;
diffvector = diff([0; eventvector; 0]);

out.start = find(diffvector == 1);
out.end = find(diffvector == -1) - 1;
% extend starts to onset of rising edge and ends to offset of falling edge
diffvector = diff([in; 0]);
eventvector = diffvector>0;
numtransients = length(out.start);
for i=1:numtransients
    n = out.start(i);
    while eventvector(n)
        n = n-1;
        if n <= 0
            n = 0;
            break
        end
    end
    out.start(i) = n+1;       
end

for i=1:numtransients
    n = out.end(i);
    while ~eventvector(n)
        n = n+1;
        if n >= numpoints
            n = numpoints;
            break
        end
    end
    out.end(i) = n;       
end

for i=1:numtransients
    out.x{i} = xtrace(out.start(i):out.end(i));
    out.y{i} = in(out.start(i):out.end(i));
end

% prune transients that are less than temporal cutoff
ind = find((out.end - out.start)>=10);
out.start = out.start(ind);
out.end = out.end(ind);
out.x = out.x(ind);
out.y = out.y(ind);

%% display results
figure
plot(in)
hold on
plot([0 numpoints],fo*[1 1],'Color',[0 0 0])
plot([0 numpoints],(fo+2*sigma)*[1 1],'Color',[1 0 0])
ylabel('data')

figure
plot(skewvalues)
ylabel('skewness')