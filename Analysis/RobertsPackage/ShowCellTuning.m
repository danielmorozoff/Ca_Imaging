% SHOWCELLSN Visualizes cell activity in a field of view as a function of
% response frequency and response amplitude
%
% ShowCellSN calculates transients for a given imaging window. Based on the
% transient amplitudes (relative to the baseline dispersion) and the
% response frequencies for a set of trials, a map of cell activity is made.
%
% type: function
%
% inputs:
%   folder: absolute path containing all imaging data
%    
% outputs:
%   none
%
% dependencies on custom functions:
%   FastFindTransients
%   FastFindfo
%   FindImagingFiles
%
% Robert Barretto, robertb@gmail.com
% 03/24/13 8:52pm   initial commit

function [] = ShowCellSN(folder)
%% load important data
[outputmeanfile, outputfiles, ~, celltracesfile,celltransientsfile,isRobert] = FindImagingFiles(folder);
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

% load data from cell traces file
temp = load(celltracesfile);
cellspatialfilters = temp.cellspatialfilters;
cellx = temp.cellx;
celly = temp.celly;
cellxrad = temp.cellxrad;
cellyrad = temp.cellyrad;
cellsig = temp.cellsig;

% framerate detection
if ~isRobert % for scanimage files
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

% setup image info
imgdata = imread(outputmeanfile);
datadim = size(imgdata);

% calculate sizes of data
numcells = size(cellsig,1);
numtrials = size(cellsig,2);
numframes = 0;
for i=1:numtrials
    numframes = max(numframes,length(cellsig{1,i}));
end
startframe = 1;
endframe = numframes;
markerframe = [];

% create a cell id map to locate cellnumber.
cellidmap =zeros(datadim(1),datadim(2));

% create an overlayimg map
overlayimg = zeros(datadim(1),datadim(2));

%% prep figure
DisplayFigure = figure;
set(DisplayFigure,'units','centimeters')
set(DisplayFigure,'position',[5 5 24.5 17])
set(DisplayFigure,'PaperOrientation','landscape')
set(DisplayFigure,'PaperUnits','centimeters','PaperPositionMode','manual','PaperPosition',[2 2 24 16])
set(DisplayFigure,'Toolbar','figure')

% image display
ImageAxes = axes;
set(ImageAxes,'units','centimeters','position',[8 .5 16 16])
set(ImageAxes,'XLim',[1 datadim(1)],'YLim',[1 datadim(2)])
hold on
DataImg = imshow(imgdata,[],'parent',ImageAxes);
MaskImg = imshow(zeros(size(imgdata)),[],'parent',ImageAxes);
set(MaskImg,'AlphaData',  ones(size(imgdata)))

% create status display
% status update
StatusText = uicontrol(DisplayFigure,'Style','text','HorizontalAlignment','left');
set(StatusText,'String','','Units','centimeters');
set(StatusText,'Position',[.5 7 7 0.5]);
% file loaded
FileLoadText = uicontrol(DisplayFigure,'Style','text','HorizontalAlignment','left');
set(FileLoadText,'String','','Units','centimeters');
set(FileLoadText,'Position',[.5 8 7 1])

%% manual keyboard for debugging
ManualButton = uicontrol(DisplayFigure,'Style','PushButton');
set(ManualButton,'String','keyboard','Units','centimeters');
set(ManualButton,'Callback',@ManualFcn);
set(ManualButton,'Position',[.5 3 2 .5]);
    % function for popping up keyboard
    function [] = ManualFcn(~,~)
        keyboard
    end

%% controls for cell selection criteria
StartingFrame = uicontrol(DisplayFigure,'Style','Edit');
set(StartingFrame,'String','1','Units','centimeters')
set(StartingFrame,'Callback',@StartingFrameCallback)
set(StartingFrame,'Position',[5.5 2 1 .5])
set(StartingFrame,'Enable','off')
    function [] = StartingFrameCallback(~,~)
        markerin = get(StartingFrame,'String');
        sloc = strfind(markerin,'s');
        if isempty(sloc) % frames were input
            markerframe = round(str2num(markerin));
        else % time was input
            markerin(sloc) = [];
            markerframe = round(str2num(markerin)/frameperiod);
        end
        startframe = max(1,markerframe);
        set(StartingFrame,'String',num2str(startframe));
        updateimage
    end
StartingFrameTxt = uicontrol(DisplayFigure,'Style','Text');
set(StartingFrameTxt,'String','start','Units','centimeters')
set(StartingFrameTxt,'Position',[3 2 2 .5])

EndingFrame = uicontrol(DisplayFigure,'Style','Edit');
set(EndingFrame,'String',num2str(numframes),'Units','centimeters')
set(EndingFrame,'Callback',@EndingFrameCallback)
set(EndingFrame,'Position',[5.5 2.5 1 .5])
set(EndingFrame,'Enable','off')
    function [] = EndingFrameCallback(~,~)
        markerin = get(EndingFrame,'String');
        sloc = strfind(markerin,'s');
        if isempty(sloc) % frames were input
            markerframe = round(str2num(markerin));
        else % time was input
            markerin(sloc) = [];
            markerframe = round(str2num(markerin)/frameperiod);
        end
        endframe = min(numframes,markerframe);
        set(EndingFrame,'String',num2str(endframe));
        updateimage
    end
EndingFrameTxt = uicontrol(DisplayFigure,'Style','Text');
set(EndingFrameTxt,'String','end','Units','centimeters')
set(EndingFrameTxt,'Position',[3 2.5 2 .5])

PercentCutoff = uicontrol(DisplayFigure,'Style','slider');
set(PercentCutoff,'Units','centimeters','Position',[3 3 4 .5])
set(PercentCutoff,'Min',0,'Max',100,'Value',50,'Callback',@PercentSliderFcn)
set(PercentCutoff,'Enable','off')
    function [] = PercentSliderFcn(~,~)
        n = round(get(PercentCutoff,'Value'));
        n = min(max(1,n),100);
        set(PercentEdit,'String',num2str(n));
        updateimage
    end
PercentTxt = uicontrol(DisplayFigure,'Style','Text');
set(PercentTxt,'String','% response','Units','centimeters')
set(PercentTxt,'Position',[3 3.5 2 .5])
PercentEdit = uicontrol(DisplayFigure,'Style','Edit');
set(PercentEdit,'String',50,'Units','centimeters');
set(PercentEdit,'Position',[5 3.5 2 .5],'Callback',@PercentFcn)
set(PercentEdit,'Enable','off')
    function [] = PercentFcn(~,~)
        n = round(str2num(get(PercentEdit,'String')));
        n = min(max(1,n),100);
        set(PercentEdit,'String',num2str(n));
        set(PercentCutoff,'Value',n);
        updateimage
    end
        
SigmaCutoff = uicontrol(DisplayFigure,'Style','slider');
set(SigmaCutoff,'Units','centimeters','Position',[3 4 4 .5])
set(SigmaCutoff,'Min',2,'Max',10,'Value',2,'Callback',@SigmaSliderFcn)
set(SigmaCutoff,'Enable','off')
    function [] = SigmaSliderFcn(~,~)
        n = round(10*get(SigmaCutoff,'Value'))/10;
        n = min(max(2,n),10);
        set(SigmaEdit,'String',num2str(n));
        updateimage
    end
SigmaTxt = uicontrol(DisplayFigure,'Style','Text');
set(SigmaTxt,'String','cutoff','Units','centimeters')
set(SigmaTxt,'Position',[3 4.5 2 .5])
SigmaEdit = uicontrol(DisplayFigure,'Style','Edit');
set(SigmaEdit,'String',2,'Units','centimeters');
set(SigmaEdit,'Position',[5 4.5 2 .5],'Callback',@SigmaFcn)
set(SigmaEdit,'Enable','off')
    function [] = SigmaFcn(~,~)
        n = round(10*str2num(get(SigmaEdit,'String')))/10;
        n = min(max(2,n),10);
        set(SigmaEdit,'String',num2str(n));
        set(SigmaCutoff,'Value',n);
        updateimage
    end

MarkerFrame = uicontrol(DisplayFigure,'Style','Edit');
set(MarkerFrame,'String','','Units','centimeters')
set(MarkerFrame,'Callback',@MarkerFrameCallback)
set(MarkerFrame,'Position',[4.5 1.5 1.5 .5])
set(MarkerFrame,'Enable','off')
    function [] = MarkerFrameCallback(~,~)
        markerin = get(MarkerFrame,'String');
        sloc = strfind(markerin,'s');
        if isempty(sloc) % frames were input
            markerframe = round(str2num(markerin));
        else % time was input
            markerin(sloc) = [];
            markerframe = round(str2num(markerin)/frameperiod);
        end
    end
MarkerFrameTxt = uicontrol(DisplayFigure,'Style','Text');
set(MarkerFrameTxt,'String','marker','Units','centimeters')
set(MarkerFrameTxt,'Position',[3 1.5 1.5 .5])

% select whether to include onsets within the range of start/end frame, or
% whether simple occurence within start/end frame
useonsetonly = 0;
FilterTxt = uicontrol(DisplayFigure,'Style','Checkbox','Value',useonsetonly);
set(FilterTxt,'Units','centimeters','String','Onset only','Position',[.5 2 2 .5])
set(FilterTxt,'Callback',@UseOnsetOnly)
    function [] = UseOnsetOnly(~,~)
        if useonsetonly
            useonsetonly = 0;
        else
            useonsetonly = 1;
        end
        updateimage;
    end

%% trial selection interface
TrialTxt = uicontrol(DisplayFigure,'Style','Text');
set(TrialTxt,'Units','centimeters','String','Selected Trials','Position',[.5 14.5 3 .5])
trialstate = ones(1,numtrials);
for v = 1:numtrials
    UseTrial(v) = uicontrol(DisplayFigure,'Style','Checkbox','String',num2str(v),'Value',1);
    set(UseTrial(v),'Units','centimeters','Position',[.5+floor((v-1)/10) 9.5+4.5-.5*rem(v-1,10) 1.0 .5]);
    set(UseTrial(v),'Callback',@UseTrialFcn);
end
    function [] = UseTrialFcn(hObject,~)
        trialchanged = str2num(get(hObject,'String'));
        trialchangedstate = get(hObject,'Value');
        trialstate(trialchanged) = trialchangedstate;
        updateimage
    end

% trialset selection interface
trialsetfile = fullfile(folder,'trialtypes.mat');
if exist(trialsetfile)
    r = load(trialsetfile);
    trialsetnames = r.trialsetnames;
    trialsetstates = r.trialsetstates;
else
    trialsetnames = {'all'};
    trialsetstates = trialstate;
end
numtrialsets = length(trialsetnames);

% Label for trialset radio button 
TrialTxt = uicontrol(DisplayFigure,'Style','Text');
set(TrialTxt,'Units','centimeters','String','Selected Set','Position',[4.5 16 2 .5])

% Create trialset radiobutton group
TrialSetGroup = uibuttongroup('Units','centimeters','Position',[4.5 10.5 3 5.5]);
for v = 1:numtrialsets
    UseTrialSet(v) = uicontrol(DisplayFigure,'Style','radio','String',trialsetnames{v},'Parent',TrialSetGroup);
    set(UseTrialSet(v),'Units','centimeters','Position',[0 4.5-.5*rem(v-1,10) 2.5 .5]);
    set(UseTrialSet(v),'Callback',@UseTrialSetFcn);
end

    function [] = UseTrialSetFcn(hObject,~)
        trialsetstate = get(hObject,'String');
        currenttrialset = find(strcmp(trialsetstate,trialsetnames));
        % assigns which trials are used
        trialstate = trialsetstates(currenttrialset,:);
        for v=1:numtrials
            set(UseTrial(v),'Value',trialstate(v))
        end            
        % assigns trialsetstateplot buttons
        for v=1:length(trialsetnames)
            set(UseTrialSetPlot(v),'Value',v == currenttrialset);
        end        
        updateimage
    end

% Label for trialsetplot radio button
TrialSetPlotTxt = uicontrol(DisplayFigure,'Style','Text');
set(TrialSetPlotTxt,'Units','centimeters','String','Plot','Position',[6.5 16 1 .5])

% Create trialsetplot
for v = 1:numtrialsets
    UseTrialSetPlot(v) = uicontrol(DisplayFigure,'Style','checkbox');
    set(UseTrialSetPlot(v),'Units','centimeters','Position',[6.75 15-.5*rem(v-1,10) .5 .5]);
    set(UseTrialSetPlot(v),'Callback',@UseTrialSetPlotCallback);
end
    % Callback for usetrialsetplot checkbox
    function [] = UseTrialSetPlotCallback(~,~)
        updateimage;
    end


%% Interface for modifying trialsets
TrialSetName = uicontrol(DisplayFigure,'Style','Edit');
set(TrialSetName,'Units','centimeters','Position',[4.5 9.5 3 .5]);
TrialSetNameAdd = uicontrol(DisplayFigure,'Style','PushButton');
set(TrialSetNameAdd,'String','add')
set(TrialSetNameAdd,'Units','centimeters','Position',[4.5 10 1.5 .5]);
set(TrialSetNameAdd,'Callback',@AddTrialSet)
TrialSetNameSub = uicontrol(DisplayFigure,'Style','PushButton');
set(TrialSetNameSub,'String','subtract')
set(TrialSetNameSub,'Units','centimeters','Position',[6.0 10 1.5 .5]);
set(TrialSetNameSub,'Callback',@SubTrialSet)
    function [] = AddTrialSet(~,~)
        % get name of trialset
        t = get(TrialSetName,'String');
        if sum(strcmp(trialsetnames,t))==0 && ~isempty(t)
            % add trials
            numtrialsets = numtrialsets+1;
            v = numtrialsets;            
            trialsetnames{v} = t;
            trialsetstates(v,:) = trialstate;
            % add radiobutton and checkbox
            UseTrialSet(v) = uicontrol(DisplayFigure,'Style','radio','String',trialsetnames{v},'Parent',TrialSetGroup);
            set(UseTrialSet(v),'Units','centimeters','Position',[floor((v-1)/10) 4.5-.5*rem(v-1,10) 2.5 .5]);
            set(UseTrialSet(v),'Callback',@UseTrialSetFcn);
            UseTrialSetPlot(v) = uicontrol(DisplayFigure,'Style','checkbox');
            set(UseTrialSetPlot(v),'Units','centimeters','Position',[6.75 15-.5*rem(v-1,10) .5 .5]);
            set(UseTrialSetPlot(v),'Callback',@UseTrialSetPlotCallback);
            % set current trialset
            set(TrialSetGroup,'SelectedObject',UseTrialSet(v))
            SaveSetNames
        end
    end
    function [] = SubTrialSet(~,~)
        % get name of trialset
        t = get(TrialSetName,'String');
        if sum(strcmp(trialsetnames,t))==1 && ~isempty(t)            
            % delete trials 
            v = find(strcmp(trialsetnames,t));
            trialsetnames(v) = [];
            trialsetstates(v,:) = [];
            delete(UseTrialSet(v))
            UseTrialSet(v) = [];
            delete(UseTrialSetPlot(v))
            UseTrialSetPlot(v) = [];
            % reposition buttons after deletion
            numtrialsets = length(trialsetnames);
            for v=1:numtrialsets
                set(UseTrialSet(v),'Units','centimeters','Position',[floor((v-1)/10) 4.5-.5*rem(v-1,10) 2.5 .5]);
                set(UseTrialSetPlot(v),'Units','centimeters','Position',[6.75 15-.5*rem(v-1,10) .5 .5]);
            end
            SaveSetNames
        end            
    end
    function [] = SaveSetNames()
        savefile = fullfile(folder,'trialtypes.mat');
        save(savefile,'trialsetnames','trialsetstates')
    end
%% function that updates the image display
    cellresponsemap = [];
    celllist = []; %initialize as global relative to updateimage
    % so that its contents can be accessed by other functions
    function [] = updateimage()
        updatestatus
        % gather criteria for selecting valid transients
        sncutoff = get(SigmaCutoff,'Value');
        pctcutoff = get(PercentCutoff,'Value');
        startframe = str2num(get(StartingFrame,'String'));
        endframe = str2num(get(EndingFrame,'String'));   

        % gather trialsets to display
        markedsets = [];
        % scan for stuff
        for v=1:length(trialsetnames)
            markedsets(v) = get(UseTrialSetPlot(v),'Value');
        end
            % if no marked sets are selected, assume all trials to be shown    
            if sum(markedsets)==0
                markedsets(1) = 1;
            end
        trialsetstouse = find(markedsets);
        %%% for each trialstate....    
        for p = 1:length(trialsetstouse)
            currenttrialstate =trialsetstates(trialsetstouse(p),:);            
            % select by valid trials and valid frames
            if useonsetonly
                newmap = onsetmap(:,logical(currenttrialstate),startframe:endframe);
            else            
                newmap = snmap(:,logical(currenttrialstate),startframe:endframe);
            end
            
            % only select cells greater than or equal to cutoff
            newmap = newmap>=sncutoff;
            % see if there is a transient in that region / convert to binary
            % signal
            cellresponsemap = (sum(newmap,3))>0;
            % sum over all trials
            numresponses = sum(cellresponsemap,2);
            validcells(p,:) = numresponses>=((pctcutoff/100)*sum(trialstate)); 
        end
        % create color value based on accepted response type
        for q=1:size(validcells,2)
            colorvalue(q) = bin2dec(num2str(validcells(:,q)'));
        end
        
        
        % create overlayimg containing colored cells
        % create cellidmap containing the id for each colored cell
        overlayimg = zeros(size(imgdata));  
        validcells = sum(validcells,1)>0;
        if isempty(validcells)
            return
        end
        celllist = find(validcells);
        cellidmap = zeros(datadim(1),datadim(2));
        for q=1:length(celllist)
            cellnumber = celllist(q);
            xrad = cellxrad(cellnumber);  % the x and y labels need to be sorted out
            yrad = cellyrad(cellnumber);
            xmin = max(cellx(cellnumber)-xrad,1);
            xmax = min(cellx(cellnumber)+xrad,size(imgdata,1));
            ymin = max(celly(cellnumber)-yrad,1);
            ymax = min(celly(cellnumber)+yrad,size(imgdata,2));
            % color overlay image by number of cellresponses
            overlayimg(ymin:ymax,xmin:xmax) = max(overlayimg(ymin:ymax,xmin:xmax),colorvalue(cellnumber)*(cellspatialfilters{cellnumber}>0));
            cellidmap(ymin:ymax,xmin:xmax) = cellnumber*(cellspatialfilters{cellnumber}>0);
        end
        % display raw image, then overlay image of cells
        hold on
        DataImg = imshow(imgdata,[],'Parent',ImageAxes);        
        customcmap = [[0 0 0];jet(2^length(trialsetstouse)-1)];
        overlayimg = ind2rgb(1+max(overlayimg,[],3),customcmap);
        MaskImg = imshow(overlayimg,'Parent',ImageAxes);        
        set(MaskImg,'AlphaData',.5)

        % create a legend
        % http://stackoverflow.com/questions/12902709/how-to-add-legend-in-i
        % magesc-plot-in-matlab
        N = size(customcmap,1);
        L = line(ones(N-1),ones(N-1));        
        set(L,{'color'},mat2cell(customcmap(2:end,:),ones(1,N-1),3));
        % make legend titles
        usedsets = trialsetnames(trialsetstouse);
        for q = 1:N-1
            a = dec2bin(q,length(usedsets)) - '0';
            legendnames{q} = [usedsets{logical(a)}];
        end
        colorlegend = legend(L,legendnames);
        hold off
        updatestatus
    end

%% function that updates the status display
    function [] = updatestatus()
        if isempty(get(StatusText,'String'))
            set(StatusText,'String','processing...');
        else
            set(StatusText,'String','');
        end
        pause(.1)
    end

%% function for mouse click based cell examination
CellNumTxt = uicontrol(DisplayFigure,'Style','Text');
set(CellNumTxt,'String','Cell Number','Units','centimeters')
set(CellNumTxt,'Position',[3 5.5 4 .5])
set(DisplayFigure,'WindowButtonDownFcn',@MouseClick)

    function MouseClick(~,~)
        x = get(ImageAxes,'CurrentPoint');
        pts = round(x(1,1:2));
        if pts(1)<1 || pts(2)<1 || pts(1)>datadim(1) || pts(2)>datadim(2)
            return
        end
        cellselected = cellidmap(pts(2),pts(1));
        if cellselected>0
            clicktype = get(gcbf,'SelectionType');
            if strcmp(clicktype, 'open') % double click  
                markedsets = [];
                % scan for trials that need to be plotted
                for v=1:length(trialsetnames)
                    markedsets(v) = get(UseTrialSetPlot(v),'Value');
                end         
                % change plot type depending on number of checked
                % trialstates
                hottrials = cellresponsemap(cellselected,:);
                if sum(markedsets) == 0
                    Plotcellsummary(folder,cellselected,[],[],1);
                elseif sum(markedsets) == 1
                    Plotcellsummary(folder,cellselected,trialstate,[startframe endframe],1,[],markerframe);                    
                elseif sum(markedsets) > 1
                    Plotcellsets(folder,cellselected,trialsetstates(find(markedsets),:),trialsetnames(find(markedsets)),[startframe endframe],1,[],markerframe);
                end
            else % single click
                set(CellNumTxt,'String',['Cell Number: ' num2str(cellselected)]);
            end
        end
    end

%% function for slider based cell examination
CellScrollButton = uicontrol(DisplayFigure,'Style','PushButton');
set(CellScrollButton,'String','Cell Viewer','Units','centimeters')
set(CellScrollButton,'Position',[.5 5.5 2 .5])
set(CellScrollButton,'Callback',@CellScrollButtonFcn);
    function CellScrollButtonFcn(~,~)
        markedsets = [];
        % scan for stuff
        for v=1:length(trialsetnames)
            markedsets(v) = get(UseTrialSetPlot(v),'Value');
        end
        
        % create an info string
        % format string for trials used
        if sum(markedsets) == 0
            tstring = '';
            for i=1:numtrials
                if trialstate(i)
                    tstring = [tstring num2str(i) ','];
                end
            end
            tstring(end) = '';
            infodata{1} = ['trials used:' tstring ];
        elseif sum(markedsets) >= 1
            for v = 1:numtrialsets
                if get(UseTrialSet(v),'Value')
                    trialsetstr = get(UseTrialSet(v),'String');
                    break
                end
            end            
            infodata{1} = ['trialset used: ' trialsetstr];            
        end    
        % format string for cutoff threshold
        sigmacutoff = num2str(get(SigmaCutoff,'Value'));
        infodata{2} = ['amplitude threshold (sigma):' sigmacutoff];
        % format string for % response threshold
        totaltrials = sum(trialstate);
        pctcutoff = get(PercentCutoff,'Value');
        mintrials = ceil(totaltrials*pctcutoff/100);
        infodata{3} = ['minimum number of responsive trials: ' num2str(mintrials) ' of ' num2str(totaltrials)];
        % format string for window of frames examined
        startframe = str2num(get(StartingFrame,'String'));
        endframe = str2num(get(EndingFrame,'String'));
        infodata{4} = ['frames used: ' num2str(startframe) ' to ' num2str(endframe)];       
        
        % spawn the right cell player parameters
        if sum(markedsets) == 0
            modularcellplayer(folder,celllist,[],[],[],[],infodata,markerframe)
        elseif sum(markedsets) == 1
            modularcellplayer(folder,celllist,markedsets,trialstate,[],[startframe endframe],infodata,markerframe)
        elseif sum(markedsets) > 1
            modularcellplayer(folder,celllist,markedsets,trialsetstates(find(markedsets),:),trialsetnames(find(markedsets)),[startframe endframe],infodata,markerframe)
        end
    end

%% function for saving plots
savename = folder( (max(strfind(folder,filesep))+1):end );
SaveName = uicontrol(DisplayFigure,'Style','Edit');
set(SaveName,'String',savename,'Units','centimeters');
set(SaveName,'Callback',@SaveNameFcn);
set(SaveName,'Position',[3 .5 4 .5])
set(SaveName,'HorizontalAlignment','left')
    function [] = SaveNameFcn(~,~)
        savename = get(SaveName,'String');
    end
SavePlotButton = uicontrol(DisplayFigure,'Style','PushButton');
set(SavePlotButton,'String','Save plot','Units','centimeters');
set(SavePlotButton,'Callback',@SavePlotFcn);
set(SavePlotButton,'Position',[.5 .5 2 .5]);
    % function for saving button
    function [] = SavePlotFcn(~,~)
        % setup print display
        PrintFigure = figure;
        set(PrintFigure,'units','centimeters')
        set(PrintFigure,'position',[5 5 24.5 17])
        set(PrintFigure,'PaperOrientation','landscape')
        set(PrintFigure,'PaperUnits','centimeters','PaperPositionMode','manual','PaperPosition',[1 1 24.5 17])
        set(PrintFigure,'color','w')

        % image display
        PrintAxes = axes;
        set(PrintAxes,'units','centimeters','position',[8 .5 16 16])
        set(PrintAxes,'XLim',[1 datadim(1)],'YLim',[1 datadim(2)])
        % original way -- overlay two images
        % hold on
        % imshow(imgdata,[],'parent',PrintAxes);
        % MaskImg = imshow(overlayimg,'parent',PrintAxes);
        % set(MaskImg,'AlphaData',.25)        
        % hold off
        imglowdata = double(imgdata)/(double(max(imgdata(:))));
        rgbdata = repmat(imglowdata,[1 1 3]);
        rgbdata(overlayimg>0) = overlayimg(overlayimg>0);
        imshow(rgbdata,[],'parent',PrintAxes)
        
        % text display
        InfoTxt = uicontrol(PrintFigure,'Style','Text');
        set(InfoTxt,'String','Cell Number','Units','centimeters')
        set(InfoTxt,'HorizontalAlignment','left')
        set(InfoTxt,'Position',[.5 .5 7 15])
        set(InfoTxt,'BackgroundColor','w')

        % format string for folder analyzed
        setname = folder( (max(strfind(folder,filesep))+1):end );
        infodata{1} = setname;
        % format string for trials used
        tstring = '';
        for i=1:numtrials
            if trialstate(i)
                tstring = [tstring num2str(i) ','];
            end
        end
        tstring(end) = '';
        infodata{3} = ['trials used:' tstring ];
        % format string for cutoff threshold
        sigmacutoff = num2str(get(SigmaCutoff,'Value'));
        infodata{5} = ['amplitude threshold (sigma):' sigmacutoff];
        % format string for % response threshold
        totaltrials = sum(trialstate);
        pctcutoff = get(PercentCutoff,'Value');
        mintrials = ceil(totaltrials*pctcutoff/100);
        infodata{7} = ['minimum number of responsive trials: ' num2str(mintrials) ' of ' num2str(totaltrials)];
        % format string for window of frames examined
        startframe = str2num(get(StartingFrame,'String'));
        endframe = str2num(get(EndingFrame,'String'));
        infodata{9} = ['frames used: ' num2str(startframe) ' to ' num2str(endframe)];
        
        set(InfoTxt,'FontSize',12)
        set(InfoTxt,'String',infodata)

        % print the graph
        printstuff = 1 ;
        if printstuff == 1
            printfolder = fullfile(folder,'plots');
            if ~exist(printfolder)
                mkdir(printfolder)
            end
            savefile = fullfile(folder,'plots',[savename '.pdf']);
            gsprinter(PrintFigure,savefile);
        end

    end

%% Pause to update display
pause(.5)
% enable controls
set(findall(DisplayFigure,'Style','Edit'),'Enable','off')
set(findall(DisplayFigure,'Style','RadioButton'),'Enable','off')
set(findall(DisplayFigure,'Style','PushButton'),'Enable','off')
set(findall(DisplayFigure,'Style','Slider'),'Enable','off')
set(findall(DisplayFigure,'Style','Checkbox'),'Enable','off')
set(findall(DisplayFigure,'Style','Text'),'Enable','off')

%% calculate transients only once
% create a nomrmalized cell signal for each set
pause(.1)
tic
if ~isempty(celltransientsfile)
    set(StatusText,'String','Calculating Transients')
    temp = load(celltransientsfile);
    normcellsigx = temp.normcellsigx;
    normcellsigy = temp.normcellsigy;
    normcellsigtransients = temp.normcellsigtransients;
    normcellsigvar = temp.normcellsigvar;
else
    set(StatusText,'String','Loading Saved Transients')
    % initialize an array of cells and trials
    normcellsigx = cell(numcells,numtrials);
    normcellsigy = normcellsigx;
    normcellsigtransients = normcellsigx;
    normcellsigvar = normcellsigx;
    % analyze each cell
    for currentcell=1:numcells
        set(StatusText,'String',['Calculating Transients for cell ' num2str(currentcell) ' of ' num2str(numcells)]);
        pause(.01)
        for currenttrial=1:numtrials
            celltrace = cellsig{currentcell,currenttrial};
            % calculate time series
            xtrace = frameperiod*(1:length(celltrace));
            % median filter
            [celltrace,~] = CalciumTraceFilter(celltrace,'median',4);
            % calculate baseline and normalize
            % note that this is slow because we run sequential skewness
            % functions twice....

            % calculate baseline
            fo = FastFindfo(celltrace);
            % normalize 
            celltrace = (celltrace-fo)./(fo); 
            % calculate transients
            [transients,v,~] = FastFindTransients(celltrace,xtrace); 

            % store in cell array
            normcellsigy{currentcell,currenttrial} = celltrace; 
            normcellsigx{currentcell,currenttrial} = xtrace;
            normcellsigtransients{currentcell,currenttrial} = transients; 
            normcellsigvar{currentcell,currenttrial} = v;
        end
    end 
    save(fullfile(folder,'celltransients.mat'),'normcellsigy','normcellsigx','normcellsigtransients','normcellsigvar');
end
% populate a master SN map and onset map
snmap = zeros(numcells,numtrials,numframes);
onsetmap = snmap;
for currentcell=1:numcells
    for currenttrial=1:numtrials
         transients = normcellsigtransients{currentcell,currenttrial};
         if ~isempty(transients.x)
             for i=1:length(transients.x);
                 snmap(currentcell,currenttrial,round(transients.x{i}/frameperiod)) = transients.sn{i};
                 onsetmap(currentcell,currenttrial,transients.start(i)) = transients.sn{i};             
             end
         end
    end
end
toc

% calculate display of cells
updateimage

% enable controls
set(findall(DisplayFigure,'Style','Edit'),'Enable','on')
set(findall(DisplayFigure,'Style','RadioButton'),'Enable','on')
set(findall(DisplayFigure,'Style','PushButton'),'Enable','on')
set(findall(DisplayFigure,'Style','Slider'),'Enable','on')
set(findall(DisplayFigure,'Style','Checkbox'),'Enable','on')
set(findall(DisplayFigure,'Style','Text'),'Enable','on')
set(StatusText,'String','')
pause(.1)
    
end % end main function