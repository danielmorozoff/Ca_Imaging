function [outputtrace,outputimage] = Midtemporalica(data,plotfigure)
debug = 0;
if debug
    load('testdata1.mat')
    data = volumestack;
end
datasize = size(data);
%% optional filters
% data = KalmanPadStackFilter(data);
% for i=1:datasize(3)
%     data(:,:,i) = medfilt2(data(:,:,i),[5 5]);
% end

X = reshape(double(data),(datasize(1)*datasize(2)),datasize(3));
A =[]; c = 0; numICA = 4; numPCA = 4;
iterlimit = 4;
while isempty(A) && c<iterlimit % perform temporal ICA with PCA to reduce time to 4 dimensions
    c = c+1;
    [~,A,W] = fastica(X,'approach','symm','lastEig',numPCA,'verbose','off','numOfIC',numICA);   
end
% if convergence could not be attained
if c==iterlimit
    outputimage = [];
    outputtrace = [];
    return
end

% use the skewness to correct the directionality of traces (as firing
% transients have positive skewness)

%T = X'*A;  % basically saying the calcium transients are a mixture of independent temporal signals
%icaskew = skewness(T);
%T = T.*repmat(sign(icaskew),[size(T,1),1]);
icaAskew = skewness(A,0,1);
A = A.*repmat(sign(icaAskew),[size(A,1),1]);

% a way to prioritize the actual cell transient (this is not good however)
%[~,icaorder] = sort(icaskew,'descend');
% another way to prioritize actual cell transient based on size of 2sigma
% transients (this prioritizes highest responder -- but sometimes adjacent
% cells fuck it up)
[~,icaorder] = sort(sum(zscore(A)>2),'descend');

% another method is to look at each image, throw out border images, and
% then take the highest zscore image
excludeimg = zeros(1,numICA);
for i=1:numICA
    % normalize, and then threshold
    tempimg = A(:,i);
    tempimg = (tempimg-min(tempimg))/(max(tempimg) - min(tempimg));    
    tempimg(tempimg<graythresh(tempimg)) = 0;

    % extract largest blob
    cellmap = bwlabel(logical(reshape(tempimg,[size(data,1),size(data,2)])));
    numblobs = max(cellmap(:));
    blobareas = zeros(1,numblobs);
    for blob = 1:numblobs
        blobareas(blob) = sum(sum(cellmap==blob));
    end
    [icabinaryarea{i}, bigblobind] = max(blobareas);
    icabinaryimage{i} = cellmap==bigblobind;    
    % exclude image if main blob is on border
    if (sum(icabinaryimage{i}(1,:))>0) || (sum(icabinaryimage{i}(end,:))>0) || (sum(icabinaryimage{i}(:,1))>0) || (sum(icabinaryimage{i}(:,end))>0 )
        excludeimg(i) = 1;
    elseif (icabinaryarea{i} < .125 * datasize(1)*datasize(2))
        excludeimg(i) = 2;
    end
    
    % reconstruct ICA images (full weightings for all pixels)
    icabaseimage{i} = reshape(A(:,i),[size(data,1),size(data,2)]); 
    % some smoothing
    icareducedimage{i} = icabaseimage{i} .* bwmorph(icabinaryimage{i},'close');  
    % calculate trace (reduced weightings for binary pixels of cell)
    icareducedtrace{i} = X'*icareducedimage{i}(:);    
end

% visualize the process
visualize = 0;
if visualize
%   show the raw frames in montage
%   stack2montage(data);
    % prepare figure
    ICAFigure = figure(plotfigure);
    clf(ICAFigure);
    hold on
    for tracer = 1:numICA
       subplot(numICA,4,4*tracer-2)
       plot(T(:,tracer),'Color',[0 0 0]);
       axis off
       subplot(numICA,4,4*tracer)
       plot(icareducedtrace{tracer},'Color',[0 0 1]);
       axis off
       
       subplot(numICA,4,4*tracer-3)
       imshow(icabaseimage{tracer},[])
       ylabel(num2str(tracer));
 
       subplot(numICA,4,4*tracer-1)
       imshow(icabinaryimage{tracer},[])
       if excludeimg(tracer) == 1
           ylabel('border')
       elseif excludeimg(tracer) == 2
           ylabel('small')
       end
    end
end

% select one index for output
if sum(~excludeimg) == 0
    outputimage = [];
    outputtrace = [];
    return
else
    selectedindex = (sum(zscore(A)>2) .* ~excludeimg);
    [~,selectedindex] = max(selectedindex);
end

visualize = 0;
if visualize
    subplot(numICA,4,4*selectedindex)
    xrect = get(gca,'XLim'); yrect = get(gca,'YLim');
    rectangle('Position',[xrect(1),yrect(1),xrect(2),yrect(2)],'EdgeColor',[0 1 0])
end

% format image
outputimage = icareducedimage{selectedindex};
outputtrace = icareducedtrace{selectedindex};
    
% visualize = 0;
% if visualize
%     figure;
%     subplot(2,1,1)
%     plot(outputtrace);
%     subplot(2,1,2)
%     imshow(outputimage,[],'InitialMagnification',400);
% end
% 
% if sum(~excludeimg)>0
%     figure;
%     imshow(bwmorph(icabinaryimage{selectedindex},'close'),'InitialMagnification',1000)
%     keyboard
% end
% 
