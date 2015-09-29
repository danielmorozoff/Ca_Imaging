basefolder = Filemacpc('C:\Users\rpjb\Dropbox\project - geniculate newdata\gcamp6s\');
% load scoring data
csvfile = fullfile(basefolder,'gcamp6s.csv');
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
    
celltypes = classifycells(celltunings,0);
for i=1:size(celltypes,1)
    if celltypes(i,2) && celltypes(i,4)
        celltypes(i,4) = 0;
    end
end

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

%% find occurences of each combination and classify each cell to an id type
occurences = zeros(1,31);
celltypenums = zeros(1,size(celltypes,1));
for i=1:length(combinations)
    matchingcells = ismember(celltypes,combinations(i,:),'rows');
    celltypenums(matchingcells) = i; % single number classifying the cell
    occurences(i) = sum(matchingcells);
end

makecartoonfcn(occurences);

%% make a bar graph of cell numbers

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

%% load cell transient data

% start points for each tastant delivery in seconds
tastantstart = [6.5 21.5 36.5 51.5 66.5];
tastantend = tastantstart+13;

% identify datasets used
datasets = unique(cellset);
for i=1:length(datasets)    
    % load celltransients
    setfolder = fullfile(basefolder,datasets{i});
    [~,~,~,~,transientsfile,~] = FindImagingFiles(setfolder);
    x = load(transientsfile);        
    datasize = size(x.normcellsigtransients);
    numcells = datasize(1);
    numtrials = datasize(2);
    
    % examine each cell
    for j=1:numcells
        tasteamplitude = nan(4,5);
        % examine the first four trials
        for k=1:4             
            celldata = x.normcellsigtransients{j,k};
            % examine each transient
            for m=1:length(celldata.x)
                onset = celldata.x{m}(1);
                tastant = find(onset>tastantstart & onset<tastantend);
                if ~isempty(tastant)
                    tasteamplitude(k,tastant) = max(tasteamplitude(k,tastant),max(celldata.y{m}));
                end
            end            
        end
        % calculate mean amplitude
        avgamplitude = nanmean(tasteamplitude);
        avgamplitude(isnan(avgamplitude)) = 0;
        
        % switch the umami and sour position
        avgamplitude(4:5) = avgamplitude(5:-1:4); 
        
        % find correct position by matching dataset and cellid
        correctcellid = find(strcmp(datasets{i},cellset) & strcmp(num2str(j),cellids));
        cellamplitudes(correctcellid,:) = celltypes(correctcellid,:).* avgamplitude;
    end   
end

%% assemble an "occurences vector" for amplitudes

meantastantamplitudes = sum(cellamplitudes)/sum(cellamplitudes>0);

figure
set(gca,'Units','centimeters');

for maintaste = 1:5    
    % find the types of combinations that contain a maintaste
    types = find(combinations(:,maintaste));
    ybox = 2;
    xleft = 0;
    ybottom = maintaste*2.5;
    for i=1:length(types)
        % determine the color to be used for box
        colorset = combinations(types(i),:);
        % determine the length of the box
            % identify the cells that have a given combination
            cellstoadd = ismember(celltypes,combinations(types(i),:),'rows');
            % sum the maintaste amplitudes
            xbox = sum(cellamplitudes(cellstoadd,maintaste));
        if i==1
            % isolate the meantastant amplitude of the singly tuned cell
            meantastantamplitudes(maintaste) = xbox/(sum(cellstoadd));
        end
        xbox = xbox/meantastantamplitudes(maintaste);        
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
set(gca,'YLim',[0 18])

%% assemble sweet mean amplitudes
tasteamplitudes = {};
for maintaste = 1:5;
    types = find(combinations(:,maintaste));
    for i=1:length(types)
        % identify the cells that have a given combination
        cellstoadd = ismember(celltypes,combinations(types(i),:),'rows');
        % save amplitudes into a struct
        tasteamplitudes{maintaste,types(i)} = cellamplitudes(cellstoadd,maintaste);
    end
end
%% make plot for broadly tuned cells
% that illustrates how amplitudes are with respect for amplitudes for each
% of the singly tuned cells

% setup figure
figure;
set(gcf,'Units','centimeters','Position',[5 5 14 14])
tastenames = {'sweet','bitter','low salt','sour','umami'};
xtaste = 1;
ytaste = 5;
doubletype = find(combinations(:,xtaste) & combinations(:,ytaste),1,'first');

binedges = [0:.25:5];
maxcells = 15;
xplot = axes;
[counts,bins]=histc(tasteamplitudes{xtaste,xtaste},binedges);
xbarplot = bar(binedges,counts,'histc');
set(xbarplot,'FaceColor',[.7 .7 .7],'EdgeColor',[.5 .5 .5]);
set(xplot,'XLim',[binedges(1),binedges(end)])
set(xplot,'YLim',[0 maxcells])
set(xplot,'XTick',[]);
set(xplot,'YTick',[]);
set(xplot,'XTickLabel',{})
set(xplot,'Units','centimeters')
set(xplot,'Position',[8 7 5 5])
xlabel(tastenames{xtaste})
set(get(gca,'XLabel'),'Position',[mean(binedges) maxcells+1.5])
ylabel(['# of singly tuned ' tastenames{xtaste} ' cells'])
box off
yplot = axes;
[counts,bins]=histc(tasteamplitudes{ytaste,ytaste},binedges);
ybarplot = barh(binedges,counts,'histc');
set(ybarplot,'FaceColor',[.7 .7 .7],'EdgeColor',[.5 .5 .5]);
set(yplot,'YLim',[binedges(1),binedges(end)])
set(yplot,'XLim',[0 maxcells])
set(yplot,'YTick',[]);
set(yplot,'XTick',[]);
set(yplot,'YTickLabel',{})
set(yplot,'Units','centimeters')
set(yplot,'Xdir','r')
set(yplot,'Position',[2 1 5 5])
set(yplot,'YDir','reverse')
ylabel(tastenames{ytaste})
set(get(gca,'YLabel'),'Position',[maxcells+.5  mean(binedges)])
xlabel(['# singly tuned ' tastenames{ytaste} ' cells'])
box off
zplot = axes;
plot(tasteamplitudes{xtaste,doubletype},tasteamplitudes{ytaste,doubletype},'o')
set(gca,'XLim',[binedges(1),binedges(end)])
set(gca,'YLim',[binedges(1),binedges(end)])
set(gca,'Units','centimeters')
set(zplot,'Position',[8 1 5 5])
set(zplot,'XTick',[binedges(1):1:binedges(end)]);
set(zplot,'YTick',[binedges(1):1:binedges(end)]);
set(zplot,'YDir','reverse')
set(zplot,'YAxisLocation','right')
xlabel('peak amplitude')
ylabel('peak amplitude')
box off