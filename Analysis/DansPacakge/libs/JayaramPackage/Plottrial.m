% PLOTTRIAL Analyzes data and creates a raster heatmap of neuronal data
% for a given imaging trial
%
% plottrials loads cell data and computes transients for each cell for a
% given trial. a plot of all cells is made, allowing the user to scroll
% through fluorescence traces for an individual trial.
% cells are then sorted based on the onset time of the
% largest transient. a raster is then made showing the sorted cells. images
% are then saved as pdfs
%
% type: function
%
% inputs:
%   folder: absolute path of folder containing imaging data
%    
% outputs:
%   none
%   creates pdfs for each trial in the imaging set
%
% dependencies on custom functions:
%   FastFindTransients
%   FindImagingFiles
%
% Robert Barretto, robertb@gmail.com
% 03/24/13 8:14pm initial commit

function [TrialFigure,RasterFigure] = Plottrial(folder,trialnumber)

[~, outputfiles, ~, celltracesfile,~,isRobert] = FindImagingFiles(folder);
% end program early if critical files not found
if isempty(celltracesfile)
    disp('celltrace file does not exist')
    return
end
if isempty(outputfiles{1})    
    disp('imagestack file does not exist')
    return
end

% load data
load(celltracesfile)

% scanimage framerate detection
if ~isRobert
    x =imfinfo(outputfiles{1});
    imgdesc = x(1).ImageDescription;
    hdr.fastZEnable = regexp(imgdesc,'scanimage.SI4.fastZEnable = ([0-9.]+)','tokens');
    hdr.fastZEnable = str2num(char(hdr.fastZEnable{1}));
    if ~hdr.fastZEnable
        hdr.scanFramePeriod = regexp(imgdesc,'scanimage.SI4.scanFramePeriod = ([0-9.]+)','tokens');
        frameperiod = str2num(char(hdr.scanFramePeriod{1}));
    else
        hdr.fastZPeriod = regexp(imgdesc,'scanimage.SI4.fastZPeriod = ([0-9.]+)','tokens');
        frameperiod = str2num(char(hdr.fastZPeriod{1}));
    end
else % for prairie files (note this assumes the first frameperiod is representative of all trials)
    [samplefolder,~,~] = fileparts(outputfiles{1});
    frameperiod = getframeperiod(samplefolder);
end

% identify number of cells
numcells = size(cellsig,1);
% create a nomrmalized cell signal for each set
for currentcell=1:numcells
    celltrace = cellsig{currentcell,trialnumber};        
    % calculate time series
    xtrace = frameperiod*(1:length(celltrace));
    % median filter
    [celltrace,~] = CalciumTraceFilter(celltrace,'median',4);
    % calculate fo, transients, and sigma (this is faster than calling
    % FastFindfo) and adjust scaling after.
    [transients,v,med] = FastFindTransients(celltrace,xtrace);
    % scale traces to normalized units
    celltrace = (celltrace - med)./med;
    for j = 1:length(transients.y)
        transients.y{j} = (transients.y{j} - med)./ med;
    end        
    % scale sigma: sigma_normalized = sigma_original / mean_original
    v = v/med;
    % store in cell array
    normcellsigy{currentcell} = celltrace; 
    normcellsigx{currentcell} = xtrace;
    normcellsigtransients{currentcell} = transients; 
    normcellsigvar{currentcell} = v;
end
traceperiod = length(celltrace);
traceperiod = ceil(traceperiod * frameperiod);


%% make plot of all activity for each trial, then make a raster of all
% cells for each trial

    % setup up figure for showing traces
    TrialFigure = figure;
    set(TrialFigure,'Name','Cell Transients','NumberTitle','off')
    set(TrialFigure,'units','centimeters','Color','w')      
    TrialAxis = axes;
    set(TrialAxis,'Units','centimeters','Box','off');    
    cmpers = .5;
    tracewidth = cmpers * traceperiod;

    % plot traces, detect and overlay transients
    cutoff = 1.75;
    for currentcell=1:numcells        
        xtrace = normcellsigx{currentcell};
        ytrace = 2*currentcell+normcellsigy{currentcell};
        ytrace(ytrace>(2*currentcell+cutoff)) = 2*currentcell+cutoff;
        vartrace = normcellsigvar{currentcell};
        plot([xtrace(1) xtrace(end)],(2*currentcell+2*vartrace)*ones(1,2),'Color',[.7 .7 .7])
        hold on
        plot(xtrace,ytrace,'Color',[0 0 0]);
        if ~isempty(normcellsigtransients{currentcell}.x)
            for i=1:length(normcellsigtransients{currentcell}.x)
                ytransient = 2*currentcell+normcellsigtransients{currentcell}.y{i};
                ytransient(ytransient>(2*currentcell+cutoff)) = 2*currentcell+cutoff;                                
                plot(normcellsigtransients{currentcell}.x{i},ytransient,'Color',[1 0 0]);
            end
        end                            
    end
    % reformat axis
%     set(TrialFigure,'Position',[2 2 tracewidth+4 10+4])
    set(TrialAxis,'Position',[2 2 tracewidth 10],'Box','off')
    set(TrialAxis,'Position',[2 2 10 15],'Box','off')
    set(TrialAxis,'XLim',[0 traceperiod],'YLim',[0 12])
    set(TrialAxis,'YTick',2*(1:numcells))
    set(TrialAxis,'YTickLabel',1:numcells)
    xlabel('Time (s)')
    title(['Trial ' num2str(trialnumber)]);

    % make a raster of stuff
    % initialize an image
    imgdata = zeros([numcells traceperiod]);
    for currentcell=1:numcells
        for i=1:length(normcellsigtransients{currentcell}.x)
            d = round(normcellsigtransients{currentcell}.x{i}/frameperiod);
            imgrow = normcellsigtransients{currentcell}.y{i};
            imgrow(imgrow>cutoff) = cutoff;
            imgdata(currentcell,d) = imgrow;
        end
    end
    [~, maxloc] = max(imgdata,[],2);
    [~,cellorder] = sort(maxloc,'ascend');
    % quick hack
    % cellorder = [1:numcells];
    
    sortedimgdata = imgdata(cellorder,:);
    % use max(1,__) to prevent zero when numcolumns < numrows
    imgsquared = imresize(sortedimgdata,[size(imgdata,1)*max(1,floor(size(imgdata,2)/size(imgdata,1))), size(imgdata,2)],'nearest');

    RasterFigure = figure;
    set(RasterFigure,'units','centimeters','Color','w','Position',[2 2 15 15])
    set(RasterFigure,'Name','Trial Raster','NumberTitle','off')
    cmap = jet(255);
    cmap(1,:) = [0 0 0];
    RasterAxis = axes;
    set(RasterAxis,'units','centimeters')
    RasterImage = imshow(255*imgsquared/cutoff,cmap,'Parent',RasterAxis); % set max at 200%
    title(['Trial ' num2str(trialnumber)]);

    set(RasterFigure,'PaperOrientation','portrait')
    set(RasterFigure,'PaperUnits','centimeters','PaperPosition',[2 2 15 15])

    printstuff =1 ;
    if printstuff == 1
        printfolder = fullfile(folder,'plots');
        if ~exist(printfolder,'dir')
            mkdir(printfolder);
        end
        savefile = fullfile(folder,'plots',['trial' num2str(trialnumber) '.pdf']);
        gsprinter(RasterFigure,savefile);
    end
    
% %% visualize the activity
% showmovie = 0;
% if showmovie == 1
% % load the current trial
% data = imreadtiffstack(file);
% 
% % adjust contrast
% cdata = imadjust(data(:));
% cdata = reshape(cdata,size(data));
% maxpx = max(cdata(:));
% 
% % for each cell with a transient
% for currentcell=1:numcells    
%     for currenttransient=1:length(celltransients{currentcell}.x)
%         activetimes = round(celltransients{currentcell}.x{currenttransient}/frameperiod);
%         cdata(cellcentroidy(currentcell)-yrad,cellcentroidx(currentcell)+[-xrad:xrad],activetimes) = maxpx;
%         cdata(cellcentroidy(currentcell)+yrad,cellcentroidx(currentcell)+[-xrad:xrad],activetimes) = maxpx;
%         cdata(cellcentroidy(currentcell)+[-yrad:yrad],cellcentroidx(currentcell)-xrad,activetimes) = maxpx;
%         cdata(cellcentroidy(currentcell)+[-yrad:yrad],cellcentroidx(currentcell)+xrad,activetimes) = maxpx;
%     end
% end
% 
% implay(cdata,4);
% end