% debug ICA

%% create dataset
data = rand([25 25 100]);

cell1x = 15;
cell1y = 15;
cell1t = 55:75;

cell2x = 15;
cell2y = 20;
cell2t = 15:40;

rad = 5;
gaussianblob = 40*fspecial('gaussian',2*rad+1,2);

for i=1:length(cell1t)
    data( (cell1x-rad):(cell1x+rad),(cell1y-rad):(cell1y+rad),cell1t(i) ) = gaussianblob + rand(2*rad+1);
end    
for i=1:length(cell2t)
    data( (cell2x-rad):(cell2x+rad),(cell2y-rad):(cell2y+rad),cell2t(i) ) = gaussianblob + rand(2*rad+1);
end
% show original data movie
stack2montage(data);
datasize = size(data);


[outputtrace,outputimage] = temporalica(data);

[outputtrace,outputimage] = spatialica(data);



% %% following fastica example
% X = reshape(data,(datasize(1)*datasize(2)),datasize(3));
% normX = zscore(X);
% [coeff,score,variances,t2] = princomp(normX,'econ');
% % dimensional reduction on time
% cutoff = 4;
% reducedX = score(:,1:cutoff);
% reducedcoeff = coeff(:,1:cutoff);
% A =[]; c = 0;
% while isempty(A)
%     [S,A,W] = fastica(reducedX','approach','symm','verbose','off','numOfIC',4);
% end
% A = (reducedcoeff*A);
% icaskew = skewness(A);
% [~,icaorder] = sort(skewness(A),'descend');
% A = (A).*repmat(sign(skewness(A)),[size(A,1),1]);
%     % format image
%     icaimage = reshape(S(icaorder(1),:),[datasize(1),datasize(2)]);
%     minimg = min(icaimage(:));
%     maximg = max(icaimage(:));
%     normicaimage = (icaimage-minimg)/(maximg - minimg);
%     bwicaimage = im2bw(normicaimage,graythresh(normicaimage));
%     bwicaimage = bwmorph(bwicaimage,'open');    
%     cellmap = bwlabel(bwicaimage);
%     outputimage = icaimage .* (cellmap==1);
%     outputtrace = A(:,icaorder(1));
%     
%     
%     
%     
%     
% %% do principal components to reduce the time, then follow with temporal
% %% ica and look at spatial weightings - THIS WORKS PRETTY WELL
% % principal components
% X = reshape(data,(datasize(1)*datasize(2)),datasize(3));
% normX = zscore(X);
% [coeff,score,variances,t2] = princomp(X,'econ');
% % dimensional reduction on time
% cutoff = 3;
% reducedX = score(:,1:cutoff);
% % temporal ICA
% [pcaica_tA,pcaica_tS] = tpicabase(reducedX);
% % identify the skewness of the image and set correct
% pcaica_tA = pcaica_tA.*repmat(sign(skewness(pcaica_tA)),[size(pcaica_tA,1),1]);    
% k = reshape(pcaica_tA,[size(data,1),size(data,2),size(pcaica_tA,2)]);
% stack2montage(k);
% xlabel('temporal ICA weighted spatial maps')
% 
% % transform back to original space
% timetraces = (pcaica_tS*coeff(:,1:cutoff)');%.*repmat(std(X),[cutoff,1])+repmat(mean(X),[cutoff,1]);
% timetraces = timetraces .* repmat(sign(skewness(timetraces')),[size(timetraces,2),1])';
% figure
% imshow(timetraces,[],'initialmagnification',400)
% xlabel('temporal ICA traces')
% %% do principal components and reduce the noise in space, then follow with
% %% temporal ica and look at spatial weightings
% % principal components
% X = reshape(data,(datasize(1)*datasize(2)),datasize(3))';
% normX = zscore(X);
% [coeff,score,variances,t2] = princomp(normX,'econ');
% % dimension reduction on space
% cutoff = 80;
% reducedX = score(:,1:cutoff);
% filteredX = reducedX*coeff(:,1:cutoff)';
% % temporal ICA
% [pcaica_tA,pcaica_tS] = tpicabase(filteredX');
% % identify the skewness of the image and set correct
% pcaica_tA = pcaica_tA.*repmat(sign(skewness(pcaica_tA)),[size(pcaica_tA,1),1]);    
% k = reshape(pcaica_tA,[size(data,1),size(data,2),size(pcaica_tA,2)]);
% stack2montage(k);
% xlabel('temporal ICA weighted maps')
% 
% %% do principal components and reduce the noise in space, then follow with
% %% spatial ica and look at spatial indendent components
% % principal components
% X = reshape(data,(datasize(1)*datasize(2)),datasize(3))';
% normX = zscore(X);
% [coeff,score,variances,t2] = princomp(normX,'econ');
% % dimension reduction on space
% cutoff = 80;
% reducedX = score(:,1:cutoff);
% pcaspatial = reducedX*coeff(:,1:cutoff)';
% stack2montage(reshape(pcaspatial',[25 25 100]));
% xlabel('PCA weighted maps')
% % temporal ICA
% [pcaica_sA,pcaica_sS] = spicabase(reducedX');
% % transform back to full spatial representation
% icaspatial = pcaica_sS*coeff(:,1:cutoff)';
% stack2montage(reshape(icaspatial',[25 25 80]))
% xlabel('spatial ICA weighted maps')
% 
% %% spatial ICA
% [pcaica_sA,pcaica_sS] = spicabase(reducedX);
% % identify the skewness of the image and set correct
% pcaica_sS = pcaica_sS.*repmat(sign(skewness(pcaica_sS)),[size(pcaica_sS,1),1]);    
% % show a montage of independent spatial components
% j = reshape(pcaica_sS',[size(data,1),size(data,2),size(pcaica_sS,2)]);
% stack2montage(j);
% xlabel('spatial independent maps')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% do principal components on temporal domain and reduce variables
% X = reshape(data,(datasize(1)*datasize(2)),datasize(3))';
% normX = zscore(X);
% montage(normX)
% %%
% [pcaica_sA,pcaica_sS] = spicabase(reducedX);
% % identify the skewness of the image and set correct
% pcaica_sS = pcaica_sS.*repmat(sign(skewness(pcaica_sS)),[size(pcaica_sS,1),1]);    
% % show a montage of independent spatial components
% j = reshape(pcaica_sS',[size(data,1),size(data,2),size(ica_sS,2)]);
% stack2montage(j);
% xlabel('spatial independent maps')
% 
% 
% %% do principal components on spatial domain and perform dimensional reduction
% % http://matlabdatamining.blogspot.com/2010/02/principal-components-analysis.html
% close all 
% X = reshape(data,(datasize(1)*datasize(2)),datasize(3));
% normX = zscore(X);
% montage(normX)
% [coeff,score,variances,t2] = princomp(normX,'econ');
% cutoff = 5;
% reducedX = score(:,1:cutoff);
% normreducedX = reducedX*coeff(:,1:cutoff)';
% normreduceddata = reshape(normreducedX,[datasize(1),datasize(2),datasize(3)]);
% 
% [pcaica_tA,pcaica_tS] = tpica(normreduceddata);
% k = reshape(pcaica_tA,[size(data,1),size(data,2),size(pcaica_tA,2)]);
% stack2montage(k);
% xlabel('temporal mix maps')
% 
% [ica_sA,ica_sS] = spica(data);
% % show a montage of independent spatial components
% j = reshape(ica_sS',[size(data,1),size(data,2),size(ica_sS,2)]);
% stack2montage(j);
% xlabel('spatial independent maps')
% 
% %% do principal components on tempoeral domain
% close all
% X = reshape(data,(datasize(1)*datasize(2)),datasize(3))';
% normX = zscore(X);
% montage(normX)
% [coeff,score,variances,t2] = princomp(normX,'econ');
% cutoff = 30;
% reducedX = score(:,1:cutoff);
% normreducedX = reducedX*coeff(:,1:cutoff)';
% normreduceddata = reshape(normreducedX',[datasize(1),datasize(2),datasize(3)]);
% stack2montage(normreduceddata)
% 
% [pcaica_tA,pcaica_tS] = tpica(normreduceddata);
% k = reshape(pcaica_tA,[size(data,1),size(data,2),size(pcaica_tA,2)]);
% stack2montage(k);
% xlabel('temporal mix maps')
% 
% [ica_sA,ica_sS] = spica(data);
% % show a montage of independent spatial components
% j = reshape(ica_sS',[size(data,1),size(data,2),size(ica_sS,2)]);
% stack2montage(j);
% xlabel('spatial independent maps')
% 
% 
% 
% %% perofrm spatial ica -- should not work 
% [ica_sA,ica_sS] = spica(data);
% % show a montage of independent spatial components
% j = reshape(ica_sS',[size(data,1),size(data,2),size(ica_sS,2)]);
% stack2montage(j);
% xlabel('spatial independent maps')
% 
% %% perform temporal ica -- should work
% [ica_tA,ica_tS] = tpica(data);
% % show a montage of independent temporal components
% k = reshape(ica_tA,[size(data,1),size(data,2),size(ica_tA,2)]);
% stack2montage(k);
% xlabel('temporal mix maps')
% 
% %%
% figure
% imshow(ica_tS,[],'InitialMagnification',400)
% 
% %%
    [a b] = spatialica(data);
    figure
    plot(a);
    xlabel('spatial ica')
    figure
    imshow(b,[],'InitialMagnification',400)
    xlabel('spatial component')
    %% 
    [c d] = temporalica(data);
    figure
    plot(c);
    xlabel('temporal ica')
    figure
    imshow(d,[],'InitialMagnification',400)
    xlabel('temporal component')