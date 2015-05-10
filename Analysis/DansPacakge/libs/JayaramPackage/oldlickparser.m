% select filename
fn = 'WS27_QUI_GO_Suc_NOGO_12_20131120190811_summary.txt';
fn = fullfile('/USB/JC_GO_NOGO/11202013',fn);

% gather data
fid = fopen(fn);
data = textscan(fid, '%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d','delimiter',',','BufSize',4096*16,'Commentstyle','trial');
fclose(fid);


%% display stuff

DisplayFigure = figure;
hold on

% data{1} gives millisecond timestamp
data{1} = data{1};

% data{2} gives framecount

% data{3} gives trial number
numtrials = max(data{3});



% 4              step num  (can ignore this)  
% 5-12       SV   (Can ignore these)
% 13-20     PV (#8 typically has water and 1-7 has tastants)
% 21           lick state (the 0s are the licks)
% 22           lick port state (1 means IN and 0 is OUT)
% 23           debugging info (Ignore this)


for trial = 1:numtrials
    % find indices corresponding to trial
    tstamps = find(data{3} == trial);
    
    % find start of trial in absolute time
    starttime = data{1}(tstamps(1));
    
    % find licks
    lickidx = find(data{3} == trial & data{21} == 0);    
    lickstarts = find(diff([1;lickidx])>1);
    lickstarts = lickstarts(1:end-1);
    lickends = find(diff(lickidx)>1);
    % plot licks
    for licknum = 1:length(lickstarts)
        x = data{1}(lickidx(lickstarts(licknum)))-starttime;
        y = 3*trial-2;
        wth = double(data{1}(lickidx(lickends(licknum))) - data{1}(lickidx(lickstarts(licknum))));
        if  wth <= 0 % in the instance that only one time point was lick
            wth = 1;
        end
        h = 1;
        rectangle('Position',[x y wth h],'EdgeColor',[0 1 0],'FaceColor',[0 1 0])
    end
    
    % find lickport state 
    lickidx = find(data{3} == trial & data{22} == 1);    
    lickstarts = find(diff([1;lickidx])>1);
    lickstarts = lickstarts(1:end-1);
    lickends = find(diff(lickidx)>1);
    % plot licks
    for licknum = 1:length(lickstarts)
        x = data{1}(lickidx(lickstarts(licknum)))-starttime;
        y = 3*trial-1;
        w = double(data{1}(lickidx(lickends(licknum))) - data{1}(lickidx(lickstarts(licknum))));
        if  w <= 0 % in the instance that only one time point was lick
            w = 1;
        end
        h = 1;
        rectangle('Position',[x y w h],'EdgeColor',[0.4 0.4 0.4],'FaceColor',[0.4 0.4 0.4])
    end
    
    % find frame period
    % find indices corresponding to trial
    fstamps = find(data{3} == trial);
    minframe = 1+min(unique(data{2}(fstamps)));
    maxframe = max(unique(data{2}(fstamps)));
    
    mintidx = min(find(data{3} == trial & data{2} == minframe));
    maxtidx = min(find(data{3} == trial & data{2} == maxframe));
    
    x = data{1}(mintidx)-starttime;
    y = 3*trial;
    w = double(data{1}(maxtidx)-data{1}(mintidx));
    h = 1;
    rectangle('Position',[x y w h],'EdgeColor',[0 0 1],'FaceColor',[0 0 1])
    
    
end
    
    