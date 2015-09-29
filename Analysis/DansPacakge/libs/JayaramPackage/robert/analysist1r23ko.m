basefolder = Filemacpc('C:\Users\rpjb\Dropbox\project - geniculate data\t1r2-3 ko\');
% load scoring data
csvfile = fullfile(basefolder,'t1r23ko.csv');
[~,data] = csv2data(csvfile);
cellset = data(:,1);
cellids = data(:,2);
celltunings = str2double(data(:,3:end));
celltypes = classifycells(celltunings,1);

% adjust for sour and umami swapping
temp = celltypes;
temp(:,4) = celltypes(:,5);
temp(:,5) = celltypes(:,4);


tuningspecificity = sum(celltypes,2);

% population a combination table
combinations = [];
for k = 1:5
    pos = nchoosek([1:5],k);
    for j = 1:size(pos,1)
        row = [0 0 0 0 0];
        row(pos(j,:)) = 1;
        combinations = [combinations;row];
    end
end
% find occurences of each combination and classify each cell to an id type
occurences = zeros(1,31);
celltypenums = zeros(1,size(celltypes,1));
for i=1:length(combinations)
    matchingcells = ismember(celltypes,combinations(i,:),'rows');
    celltypenums(matchingcells) = i; % single number classifying the cell
    occurences(i) = sum(matchingcells);
end

makecartoonfcn(occurences);