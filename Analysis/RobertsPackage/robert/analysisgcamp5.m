basefolder = Filemacpc('C:\Users\rpjb\Dropbox\project - geniculate newdata\gcamp5\');
% load scoring data
csvfile = fullfile(basefolder,'gcamp5.csv');
[~,data] = csv2data(csvfile);
cellset = data(:,1);
cellids = data(:,2);
celltunings = str2double(data(:,3:end));

    % expects
    % 1: sweet
    % 2: bitter
    % 3: salt
    % 4: umami
    % 5: sour    
    
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