% handy tool to view a small area of an image and see individual cells
% within that area

%% load parameters
x = 216;
y = 317;
xr = 10;
yr = 10;

imgwidth = 512; imgheight = 512;
finalmask = zeros(imgwidth,imgheight);


% locate cells within the requested box
cellids = find(cellcentroidx > (x-xr) & cellcentroidx < (x+xr) & cellcentroidy > (y-yr) & cellcentroidy < (y+yr));

figure;
for i=1:length(cellids)
    % create a full image
    tempmask = zeros(imgwidth,imgheight);

    % convert from local spatial coordinates to global
    idx = cellids(i);
    xpt = cellx(idx); ypt = celly(idx);
    boxx = max(xpt-cellxrad(idx),1):min(xpt+cellxrad(idx),imgwidth);
    boxy = max(ypt-cellyrad(idx),1):min(ypt+cellyrad(idx),imgheight);
    tempmask(boxy,boxx) = tempmask(boxy,boxx) + 1*(cellspatialfilters{idx} > 0);
    finalmask = finalmask+tempmask;
    % show figure
    subplot(1,length(cellids)+1,i)
    imshow(tempmask(y-yr:y+yr,x-xr:x+xr),[],'InitialMagnification',200)
pause(1)
end
    subplot(1,length(cellids)+1,length(cellids)+1)
    imshow(finalmask(y-yr:y+yr,x-xr:x+xr),[],'InitialMagnification',200)
    % black background, white single cells, and colored overlap
%    colormap([[0 0 0];[1 1 1]; jet(-1+max(finalmask(:)))])
