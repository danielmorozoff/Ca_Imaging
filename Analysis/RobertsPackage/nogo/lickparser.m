function out = lickparser(fn,fdir)
% fn = ''
%% select filename
if isempty(fn)
    fn = 'WS27_QUI_GO_Suc_NOGO_12_20131120190811.txt';
    fdir = '/Volumes/usb/JC_GO_NOGO/11202013';
    fn = fullfile(fdir,fn);
end

%% gather data
fid = fopen(fn);
data = textscan(fid, '%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d','delimiter',',','BufSize',4096*16,'Commentstyle','trial');
fclose(fid);


% data{1} gives millisecond timestamp
data{1} = data{1};

% data{2} gives framecount

% data{3} gives trial number
numtrials = max(data{3});

for t = 0:numtrials    
    trial = t+1;
    out.trialnum{trial} = t;
    
    % find indices corresponding to trial
    tstamps = find(data{3} == t);
    
    % find start of trial in absolute time
    starttime = data{1}(tstamps(1));
    out.starttime{trial} = starttime;
    
    % find start of probe delivery in absolute time
    starttrial = find(data{3}== t & data{4} == 1);
    out.starttrial{trial} = data{1}(min(starttrial))-starttime;
    
    % licks
    % find licks
    lickidx = find(data{3} == t & data{21} == 0);    
    lickstarts = find(diff([1;lickidx])>1);
    lickstarts = lickstarts(1:end-1);
    lickends = find(diff(lickidx)>1);
    out.licktimes{trial} = data{1}(lickidx(lickstarts))-starttime;
    
    % lickports    
    lickidx = find(data{3} == t & data{22} == 1);    
    if length(lickidx) == (1+sum(diff(lickidx)))
        lickstarts = 1;
        lickends = length(lickidx);
    else
        lickstarts = find(diff([1;lickidx])>1);
        lickstarts = lickstarts(1:end-1);
        lickends = find(diff(lickidx)>1);
    end
    out.lickporttimes{trial} = [data{1}(lickidx(lickstarts)) data{1}(lickidx(lickends))]-starttime;

    % frametimes
    % find frame period
    % find indices corresponding to trial
    fstamps = find(data{3} == t);
    minframe = 1+min(unique(data{2}(fstamps)));
    maxframe = max(unique(data{2}(fstamps)));
    mintidx = min(find(data{3} == t & data{2} == minframe));
    maxtidx = min(find(data{3} == t & data{2} == maxframe));    
    out.frametimes{trial} = [data{1}(mintidx) data{1}(maxtidx)]-starttime;
end
