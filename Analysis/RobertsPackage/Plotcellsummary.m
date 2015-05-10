% PLOTCELLSUMMARY Analyzes data and creates a summary plot of activity for
% a given cell
%
% plotcellsummary loads cell data and computes transients for each cell for a
% every trial. a raster is then made showing the activity for each trial. 
% an image of the cell location is generated and the entire figure is saved
% as a pdf
%
% type: function
%
% inputs:
%   folder: absolute path of folder containing imaging data
%    
% outputs:
%   none
%   creates pdf of cell activity and location
%
% dependencies on custom functions:
%   CalciumTraceFilter
%   FastFindTransients
%   FastFindfo
%
% Robert Barretto, robertb@gmail.com
% 03/24/13 8:14pm initial commit

%function [CellFigure] = Plotcellsummary(folder,cellnumber,whichtrials)
function [CellFigure] = Plotcellsummary(folder,cellnumber,whichtrials,whichframes,printplot,parentfigure,markerpos)
% folder
% cellnumber
% whichtrials
% whichframes
% printplot
% parentfigure
% markerpos

%% validate inputs
if ~exist('whichtrials','var')
    whichtrials = [];
end
if ~exist('whichframes','var')
    whichframes = [];
end
if ~exist('printplot','var') || isempty(printplot)
    printplot = 1;
end
if ~exist('parentfigure','var') || isempty(parentfigure)
    parentfigure = figure;
end
if ~exist('markerpos','var') || isempty(markerpos)
    markerpos = [];
end

%% load important data
% cellnumber = 3;
% folder = Filemacpc('F:\02132013\');

[outputmeanfile, outputfiles, ~, celltracesfile,~,isRobert] = FindImagingFiles(folder);
% end program early if critical files not found
if isempty(celltracesfile)
    disp('celltrace file does not exist')
    return
end
if isempty(outputfiles{1})    
    disp('imagestack file does not exist')
    return
end
if isempty(outputmeanfile)
    disp('mean file does not exist')
    return
end
% load data
load(celltracesfile);

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

imgdata = imread(outputmeanfile);

%% prep figure
CellFigure = figure(parentfigure);
clf(CellFigure)
currentposition = get(CellFigure,'position');
set(CellFigure,'units','centimeters','color','w')
if round(currentposition(3)) == 24 % if this plot has already been run once keep x y coordinate
    set(CellFigure,'position',currentposition)
else % make from scratch
    set(CellFigure,'position',[1 1 24 16])
end
set(CellFigure,'PaperOrientation','landscape')
set(CellFigure,'PaperUnits','centimeters','PaperPositionMode','manual','PaperPosition',[2 2 24 16])
CellImageAxes = axes;
set(CellImageAxes,'units','centimeters')
set(CellImageAxes,'position',[8 0 16 16])
DataImg = imshow(imgdata,[],'parent',CellImageAxes);
% hold on
% ContourImg = imshow(imgdata,[],'parent',CellImageAxes);
TrialAxis = axes;
set(TrialAxis,'Units','centimeters','Box','off')
set(TrialAxis,'position',[1 1 6 14])
% %% overlay image (for display only)
% overlayimg = zeros(size(imgdata));
% overlayimg(cellx(cellnumber)+(-xrad:xrad),celly(cellnumber)+(-yrad:yrad)) = cellspatialfilters{cellnumber};
% contourmask = ind2rgb(overlayimg>0,[[0 0 0];[1 0 0]]);
% hold on
% DataImg = imshow(imgdata,[],'parent',CellImageAxes);
% ContourImg = imshow(contourmask,'parent',CellImageAxes);
% set(ContourImg,'AlphaData',  .8*(max(contourmask,[],3) > 0 ));
% hold off

%% integrated color image (to print properly)
i = cellnumber;
imgwidth = size(imgdata,1);
imgheight = size(imgdata,2);
overlayimg = zeros(imgwidth,imgheight);
xpt = cellx(i); ypt = celly(i);
% convert from local spatial coordinates to global
boxx = max(xpt-cellxrad(i),1):min(xpt+cellxrad(i),imgwidth);
boxy = max(ypt-cellyrad(i),1):min(ypt+cellyrad(i),imgheight);
overlayimg(boxy,boxx) = (cellspatialfilters{cellnumber})>0;
% convert to an image of double type
imglowdata = double(imgdata)/(double(max(imgdata(:))));
rgbdata = repmat(imglowdata,[1,1,3]);
rgbdatar = rgbdata(:,:,1);
rgbdatar(overlayimg>0) = .6;
rgbdatag = rgbdata(:,:,2);
rgbdatag(overlayimg>0) = .3;
rgbdatab = rgbdata(:,:,3);
rgbdatab(overlayimg>0) = .3;
rgbdata(:,:,1) = rgbdatar;
rgbdata(:,:,2) = rgbdatag;
rgbdata(:,:,3) = rgbdatab;
DataImg = imshow(rgbdata,[],'parent',CellImageAxes);

%% calculate traces
% determine number of trials to present
if isempty(whichtrials) %if user didn't specify, select all trials
    celltrials = 1:size(cellsig,2);
else
    celltrials = find(whichtrials);
end

% create a nomrmalized cell signal for each set
traceperiod = 0;
for q=1:length(celltrials)
    currenttrial = celltrials(q);
    celltrace = cellsig{cellnumber,currenttrial};
    % calculate time series
    xtrace = frameperiod*(1:length(celltrace));
    % median filter
    [celltrace,~] = CalciumTraceFilter(celltrace,'median',4);
    % calculate fo, transients, and sigma (this is faster than calling
    % FastFindfo) and adjust scaling after.
    [transients,v,med] = FastFindTransients(celltrace,xtrace);
    % scale traces to normalized units
    celltrace = (celltrace-med)./(med); 
    for j = 1:length(transients.y)
        transients.y{j} = (transients.y{j} - med)./ med;
    end        
    % scale sigma: sigma_normalized = sigma_original / mean_original
    v = v/med;
    % store in cell array
    normcellsigy{currenttrial} = celltrace; 
    normcellsigx{currenttrial} = xtrace;
    normcellsigtransients{currenttrial} = transients; 
    normcellsigvar{currenttrial} = v;
    traceperiod = max(traceperiod,length(celltrace));
end
traceperiod = ceil(traceperiod * frameperiod);

%% make plot of all activity for each trial
cutoff = 1.75;
for q=1:length(celltrials)
    currenttrial = celltrials(q);
    % offset each trial
    xtrace = normcellsigx{currenttrial};    
    ytrace = 2*q+normcellsigy{currenttrial};
    ytrace(ytrace>(2*q+cutoff)) = 2*q+cutoff;        
    vartrace = normcellsigvar{currenttrial};
    plot([xtrace(1) xtrace(end)],(2*q+2*vartrace)*ones(1,2),'Color',[.7 .7 .7])
    hold on
    plot(xtrace,ytrace,'Color',[0 0 0]);
    if ~isempty(normcellsigtransients{currenttrial}.x)
        for i=1:length(normcellsigtransients{currenttrial}.x)
            ytransient = 2*q+normcellsigtransients{currenttrial}.y{i};
            ytransient(ytransient>(2*q+cutoff)) = 2*q+cutoff;
            if isempty(whichframes) % plot regular red for all transients
                plot(normcellsigtransients{currenttrial}.x{i},ytransient,'Color',[1 0 0]);
            else
                trialonset = normcellsigtransients{currenttrial}.x{i}(1);
                if trialonset<frameperiod*(whichframes(1)) || trialonset>frameperiod*(whichframes(2))
                    plot(normcellsigtransients{currenttrial}.x{i},ytransient,'Color',[1 .7 .7]);
                else
                    plot(normcellsigtransients{currenttrial}.x{i},ytransient,'Color',[1 0 0]);
                end
            end
        end
    end    
end

% make a marker for frames analyzed
if ~isempty(whichframes)
    % ensure that whichframes are within the plot range    
    plot([1 1]*frameperiod*whichframes(1),[0 1000],'Color',[.8 .8 .8])
    plot([1 1]*frameperiod*whichframes(2),[0 1000],'Color',[.8 .8 .8])
end
% make a marker for selected frames
if ~isempty(markerpos)
    for m=1:length(markerpos)
        plot([1 1]*frameperiod*markerpos(m),[0 1000],'Color',[.8 .8 1])
    end
end

% reformat axis
set(TrialAxis,'Box','off')
if length(celltrials)<=10
    set(TrialAxis,'XLim',[0 traceperiod],'YLim',[0 21])
else
    set(TrialAxis,'XLim',[0 traceperiod],'YLim',[0 2*length(celltrials)+1])
end
set(TrialAxis,'YTick',2*[1:length(celltrials)])
set(TrialAxis,'YTickLabel',celltrials)
xlabel('Time (s)')
ylabel('Trial Number')
title(['Cell Number: ' num2str(cellnumber)])

if printplot == 1
    printfolder = fullfile(folder,'plots');
    if ~exist(printfolder)
        mkdir(printfolder)
    end
    savefile = fullfile(folder,'plots',['cell' num2str(cellnumber) '.pdf']);
    gsprinter(CellFigure,savefile);
end
