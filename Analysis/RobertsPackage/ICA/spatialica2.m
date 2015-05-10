function [outputtrace,outputimage] = temporalicafinal(data)
load('testdata1.mat')
data = volumestack;
datasize = size(data);
%% optional filters
% data = KalmanPadStackFilter(data);
% for i=1:datasize(3)
%     data(:,:,i) = medfilt2(data(:,:,i),[5 5]);
% end

X = reshape(double(data),(datasize(1)*datasize(2)),datasize(3));
% dimensional reduction on time
cutoff = 4;
reducedX = score(:,1:cutoff);
reducedcoeff = coeff(:,1:cutoff);
A =[]; c = 0;
while isempty(A)
    [S,A,W] = fastica(reducedX','approach','symm','verbose','off','numOfIC',4);
end
A = (reducedcoeff*A);
% use the skewness to correct the directionality of traces (as firing
% transients have positive skewness)
icaskew = skewness(A);
A = (A).*repmat(sign(icaskew),[size(A,1),1]);
icaSskew = skewness(S,0,2);
S = S.*repmat(sign(icaSskew),[1,size(S,2)]);

% a way to prioritize the actual cell transient (this is no good however)
[~,icaorder] = sort(icaskew,'descend');
% another way to prioritize actual cell transient

Simg = reshape(S',[datasize(1),datasize(2),cutoff]);

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
    
    visualize = 1;
    if visualize == 1;
       stack2montage(data);
       stack2montage(reshape(S',[size(data,1),size(data,2),size(S,1)]));
       figure
       hold on
       for tracer = 1:length(icaorder)
           plot(A(:,tracer)+2*(tracer-1))
       end
    end
    keyboard
end

