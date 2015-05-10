function [centers,radii]=segmentSingleProjection(fName_segmentorImg)
disp('Processing Segmentor image...');
    I = imread(strcat(fName_segmentorImg,'.tif'));
    % tI=imadjust(I,[.5;.6],[0,1]);
    normI = double(I)./double(max(I(:)));
    tI=imadjust(normI);
    J = tI;
    
%     J = adapthisteq(I);

    %Wavelet Filtering of noise
    rng default;
    [thr,sorh,keepapp] = ddencmp('den','wv',J);
    xd = wdencmp('gbl',J,'sym4',2,thr,sorh,keepapp);

    % colormap('gray'); 
    % subplot(131);
    % imagesc(J);title('Base Image');
    % subplot(132);

    %Image filtering
    K = wiener2(xd,[3 3]);
    K = medfilt2(K,[4,4]); !This basically thresholds the amount of detected singla quite significally!
    % imagesc(K);title('Denoised x2 Image')
    % subplot(133);

    % Morphological modifications of the image
    Nhood = strel([0 1 0;1 1 1;0 1 0]);
    closedK = K;
    for i=1:3
        closedK = imclose(closedK,Nhood);
    end
    % imagesc(closedK); title('Closed denoised');
    % figure
    % subplot(121);
    se = strel('disk', 10, 0);
    tophatFilteredK = closedK;
    for i=1:3  
        tophatFilteredK = imtophat(tophatFilteredK,se);
    end
    % imagesc(tophatFilteredK); title('Tophat');
    % subplot(122);
    se =strel('disk',2,0);
    dilatedFilteredK =tophatFilteredK;
    for i=1:2
        dilatedFilteredK = imdilate(dilatedFilteredK,se);
    end
    % imagesc(dilatedFilteredK); title('Dilated');
    se =strel('disk',7,0);
    IM2 = imopen(dilatedFilteredK,se);


    %rescale 0 - 1 
    IM2 =  IM2./max(max(IM2));
    %threshold image
    IM2(IM2 < .15) = 0;

    % figure
    % subplot(121);
    % imagesc(J);title('Base Image');
    % subplot(122);
    % imagesc(IM2); title('Post processing - paper');

    %find circles
    radius = [1 100];
    [centers, radii, metric] = imfindcircles(IM2,radius,'ObjectPolarity','bright');

    % processing step

    cell_signals = {};
    sI = size(I);
    
    subplot(1,2,1)
    imagesc(I);
    viscircles(centers, radii,'EdgeColor','w');
    
    subplot(1,2,2)
    imagesc(IM2);
    viscircles(centers, radii,'EdgeColor','w');
    
    figure;
%     rawI = imread('Image_Registration_4_s10_2015_03_26_main_016.tif');
%     rI = imadjust(rawI);
    imshowpair(I,IM2,'montage');
    viscircles(centers, radii,'EdgeColor','w');
      colormap('jet')
    disp('Segmentor image processed');

% imagename is being generted in fullPathExecution script 'Image_Registration_4_s11_2015_03_26_main_001.tif';
% this loop will be exectured by the full path a number of times to
% basically map the cells to the trials.
% disp('Add stack data based on segmented image...');
% 
% info = imfinfo(imageName);
% num_images = numel(info);
%  
% for k = 1:num_images
%     data = imread(imageName,k);
%     for i = 1:length(centers)
%         tempI = data ;
%         cent_x = centers(i,1);
%         cent_y = centers(i,2);
%         r= radii(i);
% 
%         [x,y]=meshgrid(-(cent_x-1):(sI(1,1)-cent_x),-(cent_y-1):(sI(1,2)-cent_y));
%         c_mask=((x.^2+y.^2)<=r^2);
%         tempI(~c_mask) = 0;
% 
%         cell_signals{k,i,1} = [ centers(i,1),centers(i,2)];
%         cell_signals{k,i,2} =   radii(i);
%         cell_signals{k,i,3} = sum(sum(tempI));
% 
%     end
% end
% disp('Finished processing stack');
end