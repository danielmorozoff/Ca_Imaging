% load data from the first cell chunk of ANALYSISANIMALFREQDISTWITHSUCROSE

%% load raw cell data
a = load('/Users/rpjb/Dropbox/project - neuroresponse software/robertmanualscorethy1.mat');
b = load('/Users/rpjb/Dropbox/project - neuroresponse software/robertmanualscoreaav.mat');
c = load('/Users/rpjb/Dropbox/project - neuroresponse software/robertmanualscoreaav3.mat');
d = load('/Users/rpjb/Dropbox/project - neuroresponse software/robertmanualscorepharm.mat');

% group all sets together
cellstims = [a.cellstims b.cellstims c.cellstims d.cellstims];
cellnames = [a.cellname b.cellname c.cellname d.cellname];
cellscore = [a.cellscore(:,1:7); b.cellscore(:,1:7); c.cellscore ;d.cellscore];

%%
for i=1:length(cellstims)
%    if length(find(strcmp(cellstims{i},'sweet'))) == 2
        % assign sweet
        celltunings(i,1) = cellscore(i,find(strcmp(cellstims{i},'sweet'),1,'first'));
        % assign bitter
        celltunings(i,2) = cellscore(i,find(strcmp(cellstims{i},'bitter'),1,'first'));
        % assign low salt
        celltunings(i,3) = cellscore(i,find(strcmp(cellstims{i},'low salt'),1,'first'));
        % assign umami
        celltunings(i,4) = cellscore(i,find(strcmp(cellstims{i},'umami'),1,'first'));
        % assign sour
        celltunings(i,5) = cellscore(i,find(strcmp(cellstims{i},'sour'),1,'first'));
        % assign sweet (sucrose)
        celltunings(i,7) = cellscore(i,find(strcmp(cellstims{i},'sweet'),1,'last'));
%     else
%         celltunings(i,:) = 0;        
%     end
end
celltypes = classifycells(celltunings,1);

% identify how many animals fit our category
cellsubset = find(sum(celltypes,2)>0);
aa = cellnames(cellsubset);
for i=1:length(aa)
    aa{i} = aa{i}(1:end-3);
end
numanimals = length(unique(aa));


% adjust for sour and umami swapping
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

occurences(2) = occurences(2) + occurences(11);
occurences(11) = 0;
occurences(6) = occurences(6) + occurences(17);
occurences(17) = 0;

makecartoonfcn(occurences);

%% make a bar graph

figure
set(gca,'Units','centimeters');

for maintaste = 1:5
    types = find(combinations(:,maintaste));
    ybox = 2;
    xleft = 0;
    ybottom = maintaste*2.5;
    for i=1:length(types)
        colorset = combinations(types(i),:);
        xbox = occurences(types(i));
        if xbox>0
            makebox(xleft, ybottom, xbox, ybox, colorset)
            xleft = xleft+xbox;
        end
    end
end
textcolor = [1 1 1];
backgroundcolor = [0 0 0];
set(gcf,'Color',backgroundcolor)
set(gca,'Color',backgroundcolor)
%set(gca,'Color',[0 0 0])
set(gca,'YTick',1+[1:5]*2.5)
set(gca,'YTickLabel',{'sweet','bitter','salt','sour','umami'})
set(gca,'YColor',textcolor)
set(gca,'XColor',textcolor)



%%