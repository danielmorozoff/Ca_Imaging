function [outputtrace,outputimage] = spatialica(data)
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

X = reshape(double(data),(datasize(1)*datasize(2)),datasize(3))';
A =[]; c = 0; numICA = 4; numPCA = 4;
while isempty(A) && c<10 % perform spatial ICA with PCA to reduce space to 4 dimensions.
    c = c+1;
    [S,A,W] = fastica(X,'approach','symm','lastEig',numPCA,'verbose','off','numOfIC',numICA);
end
% if convergence could not be attained
if c==10
    outputimage = [];
    outputtrace = [];
    return
end
    
% use the skewness to correct the directionality of traces (as firing
% transients have positive skewness), and images
icaskew = skewness(A);
A = A.*repmat(sign(icaskew),[size(A,1),1]);
icaSskew = skewness(S,0,2);
S = S.*repmat(sign(icaSskew),[1,size(S,2)]);

% a way to prioritize the actual cell transient (this is no good however)
% [~,icaorder] = sort(icaskew,'descend');
% another way to prioritize actual cell transient based on size of 2sigma
% transients
[~,icaorder] = sort(sum(zscore(A)>2),'descend');
% if necessary to visualize
visualize = 0;
if visualize
   stack2montage(data);
   stack2montage(reshape(S',[size(data,1),size(data,2),size(S,1)]));
   figure
   hold on
   for tracer = 1:numICA
       subplot(numICA,1,tracer)
       plot(A(:,tracer));
       ylabel(num2str(tracer));
   end
end   
    
% format image
icaimage = reshape(S(icaorder(1),:),[datasize(1),datasize(2)]);
minimg = min(icaimage(:));
maximg = max(icaimage(:));
normicaimage = (icaimage-minimg)/(maximg - minimg);
bwicaimage = im2bw(normicaimage,graythresh(normicaimage));
bwicaimage = bwmorph(bwicaimage,'open');    
cellmap = bwlabel(bwicaimage);
outputimage = icaimage .* (cellmap==1);
outputtrace = A(:,icaorder(1));

visualize = 0;
if visualize
    figure;
    subplot(2,1,1)
    plot(outputtrace);
    subplot(2,1,2)
    imshow(outputimage,[],'InitialMagnification',400);
end


