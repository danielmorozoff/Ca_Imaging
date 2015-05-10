% script to discover which points are closest to each other and prunes

%% load data and parameters
%load('/Volumes/usb/JC_GO_NOGO/11202013/celldata.mat')
load([dropboxroot '/Share (from jayaram)/zoomed_out/celldata.mat'])
mindist = max([cellxrad cellyrad]);
imgwidth = 512; imgheight = 512;

%% calculate a distance matrix and cells of interest
distmat = squareform(pdist([cellcentroidx' cellcentroidy']));
pcells = find(distmat<mindist);
% S is a sorted list of cell pairs
% R gives the second partner of the pair
[R S] = ind2sub(size(distmat),pcells);

finalmask = zeros(imgwidth,imgheight);
savelist = zeros(1,length(cellcentroidx));

%% examine cell pairs of interest
for i=1:length(cellcentroidx)
    % locate indices of pairs
    overlappingidx = find(S==i);    
    numcells = length(overlappingidx);

%   % create a simple display of the spatial filters
%     for j=1:length(overlappingidx)
%         idx = R(overlappingidx(j))
%         subplot(1,length(overlappingidx),j);
%         imshow(cellspatialfilters{idx})
%     end

    % identify where on the image cell belongs    
    contourmask = [];
    for j=1:numcells
        idx = R(overlappingidx(j));
        contourmask(:,:,j) = zeros(imgwidth,imgheight);
        xpt = cellx(idx); ypt = celly(idx);
        % convert from local spatial coordinates to global
        boxx = max(xpt-cellxrad(idx),1):min(xpt+cellxrad(idx),imgwidth);
        boxy = max(ypt-cellyrad(idx),1):min(ypt+cellyrad(idx),imgheight);
        contourmask(boxy,boxx,j) = contourmask(boxy,boxx,j) + 1*(cellspatialfilters{idx} > 0);  
    end
    
%     % create an image
%     summask = sum(contourmask,3);
%     overlap = [];
%     for j=1:numcells
%         idx = R(overlappingidx(j));
%         overlap(j) = sum(sum(summask == numcells & contourmask(:,:,j) == 1))/cellarea(idx);
%     end    
%     [overlappct,overlapidx] = max(overlap);

%     % check cells that have smaller overlap
%     if overlappct<.75
%         for j=1:numcells
%             subplot(1,1+length(overlappingidx),j);
%             imshow(squeeze(contourmask(:,:,j)))    
%         end
%         subplot(1,j+1,j+1);
%         imshow(summask,[0 numcells])
%     end

    % if there is only one cell (the ith cell) save
    if numcells == 1
        savelist(i) = 1;
        continue
    end
    
    % find the position of the ith auto-comparison
    autoidx = find(S==i & R==i);    
    autojidx = find(overlappingidx==autoidx);
    ifraction = []; jfraction = [];
    overlap = []; discard = [];    
    for j=1:numcells
        idx = R(overlappingidx(j));
        if j == autojidx % skip the auto-comparison            
            ifraction(j) = 0;
            jfraction(j) = 0;
            discard(j) = 0;
        else % calculate the overlap between the jth cell and the cell of interest
            % calculate intersection
            overlap = sum(sum( contourmask(:,:,j) & contourmask(:,:,autojidx) ));
            ifraction(j) = overlap / cellarea(i);            
            jfraction(j) = overlap / cellarea(idx);
            % prune highly overlapping cells
            if ifraction(j) > .5 && jfraction(j) > .5
                if idx > i
                    discard(j) = 1;
                else 
                    discard(j) = 0;
                end
            % prune if j is a subset of i            
            elseif ifraction(j) <= .8 && jfraction(j) > .8
                discard(j) = 1;
            else
                discard(j) = 0;
            end
        end                
    end    
    % collect cell to save
    if sum(discard) == 0
        savelist(i) = 1;
    end
end


%% display all cells 
for i=1:length(cellcentroidx)
    xpt = cellx(i); ypt = celly(i);
    % convert from local spatial coordinates to global
    boxx = max(xpt-cellxrad(i),1):min(xpt+cellxrad(i),imgwidth);
    boxy = max(ypt-cellyrad(i),1):min(ypt+cellyrad(i),imgheight);
    finalmask(boxy,boxx) = finalmask(boxy,boxx) + 1*(cellspatialfilters{i} > 0);  
end
figure;
imshow(finalmask,[])
% black background, white single cells, and colored overlap
colormap([[0 0 0];[1 1 1]; jet(-1+max(finalmask(:)))])


figure;
compositemask = zeros(imgwidth,imgheight);
for i=1:length(cellcentroidx)
    if savelist(i)
        xpt = cellx(i); ypt = celly(i);
        % convert from local spatial coordinates to global
        boxx = max(xpt-cellxrad(i),1):min(xpt+cellxrad(i),imgwidth);
        boxy = max(ypt-cellyrad(i),1):min(ypt+cellyrad(i),imgheight);
        compositemask(boxy,boxx) = compositemask(boxy,boxx) + 1*(cellspatialfilters{i} > 0);    
        if max(max(compositemask(boxy,boxx))) > 2
            disp(i)
        end
    end
end
imshow(compositemask,[])
colormap([[0 0 0];[1 1 1]; jet(-1+max(compositemask(:)))])
