ldata = lickparser('','');
bdata = nogoparser('','');

% adjustment for the last trial
numtrials = length(ldata.trialnum);
if length(ldata.trialnum) < length(bdata.num)
    bdata.num = bdata.num(1:numtrials);
    bdata.valve = bdata.valve(1:numtrials);
    bdata.type = bdata.type(1:numtrials);
    bdata.result = bdata.result(1:numtrials);
end

%% display data by time

DisplayFigure = figure;
hold on

for trial=1:numtrials
    % plot licks
    for n=1:length(ldata.licktimes{trial})
        line(double(ldata.licktimes{trial}(n)).*[1 1],trial+[-.5 .5])
    end
    % plot GO vs NOGO
    if bdata.type(trial)
        x = -2000;
    else
        x = -1000;
    end

    % plot success or failure
    if bdata.result(trial)
        usecolor = [0 1 0];
    else
        usecolor = [1 0 0];
    end
       
    rectangle('Position',[x trial-.5 1000 1],'EdgeColor',usecolor,'FaceColor',usecolor)            
end

set(gca,'XLim',[-2500 12000])


%% display data by GO or NOGO
for type = 0:1

    DisplayFigure = figure;
    hold on

    % find GO trials
    trialstouse = find(bdata.type==type);


    numtrials = length(trialstouse);
    for j=1:numtrials
        trial = trialstouse(j);
        % plot licks
        
        for n=1:length(ldata.licktimes{trial})
            line(double(ldata.licktimes{trial}(n)).*[1 1],j+[-.5 .5])
        end
        % plot GO vs NOGO
        if bdata.type(trial)
            x = -2000;
        else
            x = -1000;
        end

        % plot success or failure
        if bdata.result(trial)
            usecolor = [0 1 0];
        else
            usecolor = [1 0 0];
        end

        rectangle('Position',[x j-.5 1000 1],'EdgeColor',usecolor,'FaceColor',usecolor)            
    end

    set(gca,'XLim',[-2500 12000])
end

%% display data by sligning on the second lick for NOGO
    DisplayFigure = figure;
    hold on

    % find NOGO trials
    trialstouse = find(bdata.type==0);


    numtrials = length(trialstouse);
    for j=1:numtrials
        trial = trialstouse(j);
        % plot licks
        
        for n=1:length(ldata.licktimes{trial})
            x = double([ldata.licktimes{trial}(n)]-ldata.starttrial{trial}); 
            disp(x)
            
            if bdata.result(trial)
                usecolor = [.4 .4 .4];
            else
                usecolor = [1 0 0];
            end
            line(x.*[1 1],2+j+[-.5 .5],'Color',usecolor);
        end
        % plot NOGO
        x = -3000;

        % plot success or failure
        if bdata.result(trial)
            usecolor = [0 1 0];
        else
            usecolor = [1 0 0];
        end

        rectangle('Position',[x  2+j-.5 1000 1],'EdgeColor',usecolor,'FaceColor',usecolor)            
    end

    % show 2s line
    line(0*[1 1],[0 numtrials+2],'Color',[0 0 0])
    line(2000*[1 1],[0 numtrials+2],'Color',[0 0 0])
    line(1500*[1 1],[0 numtrials+2],'Color',[0 0 0])
    set(gca,'XLim',[-3500 12000])
    set(gca,'YLim',[0 numtrials+2])
    