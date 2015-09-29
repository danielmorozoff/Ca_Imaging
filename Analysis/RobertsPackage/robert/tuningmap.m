basefolder = Filemacpc('C:\Users\rpjb\Dropbox\project - geniculate newdata\pkd2l1-tent\');
% load scoring data
csvfile = fullfile(basefolder,'pkd2l1-tent gcamp6.csv');
[~,data] = csv2data(csvfile);
cellset = data(:,1);
cellids = data(:,2);
celltunings = str2double(data(:,3:end));
celltypes = classifycells(celltunings);
tuningspecificity = sum(celltypes,2);

% select a dataset
currentcellid = cellset{1};
celllist = find(strcmp(cellset,currentcellid));
% load cell map data
folder = fullfile(basefolder,currentcellid);
[outputmeanfile, outputfiles, ~, celltracesfile,celltransientsfile,isRobert] = FindImagingFiles(folder);
temp = load(celltracesfile);
cellspatialfilters = temp.cellspatialfilters;
cellx = temp.cellx;
celly = temp.celly;
cellxrad = temp.cellxrad;
cellyrad = temp.cellyrad;

% create a map of tuning depth
imgdata = imread(outputmeanfile);
specificityimg = zeros(size(imgdata)); 
for q=1:length(celllist)
    cellnumber = celllist(q);
    xrad = cellxrad(cellnumber);  % the x and y labels need to be sorted out
    yrad = cellyrad(cellnumber);
    xmin = max(cellx(cellnumber)-xrad,1);
    xmax = min(cellx(cellnumber)+xrad,size(imgdata,1));
    ymin = max(celly(cellnumber)-yrad,1);
    ymax = min(celly(cellnumber)+yrad,size(imgdata,2));
    specificityimg(ymin:ymax,xmin:xmax) = max(overlayimg(ymin:ymax,xmin:xmax),tuningspecificity(q)*(cellspatialfilters{cellnumber}>0));
end

%% population a combination table
combinations = [];
for k = 1:5
    pos = nchoosek([1:5],k);
    for j = 1:size(pos,1)
        row = [0 0 0 0 0];
        row(pos(j,:)) = 1;
        combinations = [combinations;row];
    end
end
%% find occurences of each combination and classify each cell to an id type
occurences = zeros(1,31);
celltypenums = zeros(1,size(celltypes,1));
for i=1:length(combinations)
    matchingcells = ismember(celltypes,combinations(i,:),'rows');
    celltypenums(matchingcells) = i; % single number classifying the cell
    occurences(i) = sum(matchingcells);
end

%% make bars for the specificity of a given tastant
maintaste = 1;
% find combinations containing maintaste
comblist = find(combinations(:,maintaste));




