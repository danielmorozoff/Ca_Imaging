function modularmovieplayer(MovieAxis,DataImg,originaldata,triallist,frameperiod)
%% initialize filtered data variables
data = originaldata;
filtereddata = [];
dffdata = [];
dfffiltereddata = [];
datadim = size(data);
dataintensityrange = [min(data(:)) .7*max(data(:))];
set(MovieAxis,'CLim',dataintensityrange);

%% initialize variables
% frames to play
minframe = 1;
maxframe = datadim(3);
startframe = minframe;
endframe = maxframe;
currentframe = 1;
playbackperiod = .005;
    
% break warning
stopflag = 1;

%% setup figure
DisplayFigure = figure;
set(DisplayFigure,'units','centimeters','position',[5 2 8 5])
set(DisplayFigure,'Name','Movie Controls','NumberTitle','off')
%% video controls
PlayButton = uicontrol(DisplayFigure,'Style','PushButton');
set(PlayButton,'String','play','Units','centimeters');
set(PlayButton,'Callback',@PlayFcn);
set(PlayButton,'Position',[.5 .5 2 .5]);
    function [] = PlayFcn(~,~)
        stopflag = ~stopflag;
        if ~stopflag
            set(PlayButton,'String','stop')
        else
            set(PlayButton,'String','play')
        end
        set(MovieAxis,'CLim',dataintensityrange);
        % loop the playback
        while ~stopflag
            for currentframe=startframe:1:endframe
                set(CurrentFrameText,'String',currentframe);
                set(DataImg,'CData',data(:,:,currentframe));
                pause(playbackperiod)
                if stopflag == 1
                    return
                end
            end
        end
    end
% StopButton = uicontrol(DisplayFigure,'Style','PushButton');
% set(StopButton,'String','stop','Units','centimeters');
% set(StopButton,'Callback',@StopFcn);
% set(StopButton,'Position',[.5 1 2 .5]);
%     function [] = StopFcn(~,~)
%         stopflag = 1;
%     end

BackButton = uicontrol(DisplayFigure,'Style','PushButton');
set(BackButton,'String','back','Units','centimeters');
set(BackButton,'Callback',@BackFcn);
set(BackButton,'Position',[.5 1.5 2 .5]);
    function [] = BackFcn(~,~)
        currentframe = max(minframe,currentframe-1);
        set(CurrentFrameText,'String',currentframe);
        set(DataImg,'CData',data(:,:,currentframe));
    end

ForwardButton = uicontrol(DisplayFigure,'Style','PushButton');
set(ForwardButton,'String','forward','Units','centimeters');
set(ForwardButton,'Callback',@ForwardFcn);
set(ForwardButton,'Position',[.5 2 2 .5]);
    function [] = ForwardFcn(~,~)
        currentframe = min(maxframe,currentframe+1);
        set(CurrentFrameText,'String',currentframe);
        set(DataImg,'CData',data(:,:,currentframe));
    end

FrameRate = uicontrol(DisplayFigure,'Style','Edit');
set(FrameRate,'String',num2str(1/playbackperiod),'Units','centimeters')
set(FrameRate,'Callback',@FrameRateCallback)
set(FrameRate,'Position',[5.5 2 1 .5])
    function [] = FrameRateCallback(~,~)
        playbackperiod = 1/str2num(get(FrameRate,'String'));
    end
FrameRateTxt = uicontrol(DisplayFigure,'Style','Text');
set(FrameRateTxt,'String','framerate','Units','centimeters')
set(FrameRateTxt,'Position',[3 2 2 .5])

StartingFrame = uicontrol(DisplayFigure,'Style','Edit');
set(StartingFrame,'String','1','Units','centimeters')
set(StartingFrame,'Callback',@StartingFrameCallback)
set(StartingFrame,'Position',[5.5 1 1 .5])
    function [] = StartingFrameCallback(~,~)
        markerin = get(StartingFrame,'String');
        sloc = strfind(markerin,'s');
        if isempty(sloc) % frames were input
            markerframe = round(str2num(markerin));
        else
            markerin(sloc) = [];
            markerframe = round(str2num(markerin)/frameperiod);
        end        
        startframe = max(1,markerframe);
    end
StartingFrameTxt = uicontrol(DisplayFigure,'Style','Text');
set(StartingFrameTxt,'String','start','Units','centimeters')
set(StartingFrameTxt,'Position',[3 1 2 .5])

EndingFrame = uicontrol(DisplayFigure,'Style','Edit');
set(EndingFrame,'String',num2str(datadim(3)),'Units','centimeters')
set(EndingFrame,'Callback',@EndingFrameCallback)
set(EndingFrame,'Position',[5.5 1.5 1 .5])
    function [] = EndingFrameCallback(~,~)
        markerin = get(EndingFrame,'String');
        sloc = strfind(markerin,'s');
        if isempty(sloc) % frames were input
            markerframe = round(str2num(markerin));
        else
            markerin(sloc) = [];
            markerframe = round(str2num(markerin)/frameperiod);
        end                            
        endframe = min(datadim(3),markerframe);
    end
EndingFrameTxt = uicontrol(DisplayFigure,'Style','Text');
set(EndingFrameTxt,'String','end','Units','centimeters')
set(EndingFrameTxt,'Position',[3 1.5 2 .5])

CurrentFrameText = uicontrol(DisplayFigure,'Style','text','HorizontalAlignment','left');
set(CurrentFrameText,'String','1','Units','centimeters');
set(CurrentFrameText,'Position',[5.5 .5 1 0.5]);

CurrentFrameLabel = uicontrol(DisplayFigure,'Style','text');
set(CurrentFrameLabel,'String','frame','Units','centimeters');
set(CurrentFrameLabel,'Position',[3 .5 2 .5]);

%% button for keyboard
KeyboardButton = uicontrol(DisplayFigure,'Style','PushButton');
set(KeyboardButton,'String','keyboard','Units','centimeters');
set(KeyboardButton,'Callback',@KeyboardFcn);
set(KeyboardButton,'Position',[5.5 3 2 .5]);
    % function for kuwahara cells
    function [] = KeyboardFcn(~,~)
        keyboard
    end
%% processed video to display
usefilterindex = 0;
usedffindex = 0;
usedfffilterindex = 0;

UseFilterButton = uicontrol(DisplayFigure,'Style','checkbox','Value',0);
set(UseFilterButton,'String','show filter','Units','centimeters');
set(UseFilterButton,'Callback',@UseFilter)
set(UseFilterButton,'Position',[3 3 2 .5])
set(UseFilterButton,'Enable','off')
    function [] = UseFilter(~,~)
        usefilterindex = get(UseFilterButton,'Value');
        if usefilterindex
            usedffindex = 0; usefilterindex = 1; usedfffilterindex = 0;
            set(UseDFFButton,'Value',0);
            set(UseDFFFilterButton,'Value',0);
        elseif ~usefilterindex
            usefilterindex = 0;
        end
        updatevideo;
    end
UseDFFButton = uicontrol(DisplayFigure,'Style','checkbox','Value',0);
set(UseDFFButton,'String','show dff','Units','centimeters');
set(UseDFFButton,'Callback',@UseDFF)
set(UseDFFButton,'Position',[3 3.5 2 .5])
set(UseDFFButton,'Enable','off')
    function [] = UseDFF(~,~)
        usedffindex = get(UseDFFButton,'Value');
        if usedffindex
            usedffindex = 1; usefilterindex = 0; usedfffilterindex = 0;
            set(UseFilterButton,'Value',0);
            set(UseDFFFilterButton,'Value',0);            
        elseif ~usedffindex
            usedffindex = 0;
        end
        updatevideo;
    end
UseDFFFilterButton = uicontrol(DisplayFigure,'Style','checkbox','Value',0);
set(UseDFFFilterButton,'String','show both','Units','centimeters');
set(UseDFFFilterButton,'Callback',@UseDFFFilter)
set(UseDFFFilterButton,'Position',[3 4 2 .5])
set(UseDFFFilterButton,'Enable','off')
    function [] = UseDFFFilter(~,~)
        usedfffilterindex = get(UseDFFFilterButton,'Value');
        if usedfffilterindex
            usedffindex = 0; usefilterindex = 0; usedfffilterindex = 1;
            set(UseFilterButton,'Value',0);
            set(UseDFFButton,'Value',0);            
        elseif ~usedfffilterindex
            usedfffilterindex = 0;
        end
        updatevideo;
    end

    % subfunction that checks the dff / filter indices and updates the data
    % range
    function [] = updatevideo()
        if usedfffilterindex
            data = dfffiltereddata;
            dataintensityrange = [0 1.0];
        elseif usedffindex
            data = dffdata;
            dataintensityrange = [0 1.75];        
        elseif usefilterindex
            data = filtereddata;
        %    dataintensityrange = [min(data(:)) max(data(:))];        
            dataintensityrange = [min(data(:)) .4*max(data(:))];        
        elseif ~usedfffilterindex && ~usedffindex && ~usefilterindex
            data = originaldata;
            dataintensityrange = [min(data(:)) .4*max(data(:))];        
        end
        set(MovieAxis,'CLim',dataintensityrange)
    end  

FilterButton = uicontrol(DisplayFigure,'Style','PushButton');
set(FilterButton,'String','calc filter','Units','centimeters');
set(FilterButton,'Callback',@FilterFcn);
set(FilterButton,'Position',[.5 3 2 .5]);
    function [] = FilterFcn(~,~)
        [filtereddata] = KalmanPadStackFilter(originaldata);
        set(UseFilterButton,'Enable','on')
        set(DFFFilterButton,'Enable','on')
        set(FilterButton,'Enable','off')
    end
DFFButton = uicontrol(DisplayFigure,'Style','PushButton');
set(DFFButton,'String','calc dff','Units','centimeters');
set(DFFButton,'Callback',@DFFFcn);
set(DFFButton,'Position',[.5 3.5 2 .5]);
    function [] = DFFFcn(~,~)
        mediandata = double( repmat(median(data,3),[1 1 maxframe]) );
        dffdata = (double(originaldata) - mediandata) ./ (mediandata);
        set(UseDFFButton,'Enable','on')
        set(DFFButton,'Enable','off')
    end

DFFFilterButton = uicontrol(DisplayFigure,'Style','PushButton');
set(DFFFilterButton,'String','calc dfffilt','Units','centimeters');
set(DFFFilterButton,'Callback',@DFFFilterFcn);
set(DFFFilterButton,'Position',[.5 4 2 .5]);
set(DFFFilterButton,'Enable','off');
    function [] = DFFFilterFcn(~,~)
        mediandata = double( repmat(median(filtereddata,3),[1 1 maxframe]) );
        dfffiltereddata = (double(filtereddata) - mediandata) ./ (mediandata);        
        set(UseDFFFilterButton,'Enable','on')
        set(DFFFilterButton,'Enable','off')
    end

%% trial indicator



TrialUsed = uicontrol(DisplayFigure,'Style','Edit');
set(TrialUsed,'String','1','Units','centimeters')
set(TrialUsed,'Callback',@TrialUsedCallback)
set(TrialUsed,'Position',[5.5 3.5 1 .5])
    function [] = TrialUsedCallback(~,~)
        % reset all buttons
        set(FilterButton,'Enable','on')
        set(DFFButton,'Enable','on')
        set(DFFFilterButton,'Enable','off')
        set(UseFilterButton,'Value',0,'Enable','off');
        set(UseDFFButton,'Value',0,'Enable','off');            
        set(UseDFFFilterButton,'Value',0,'Enable','off');
        usefilterindex = 0;
        usedffindex = 0;
        usedfffilterindex = 0;
        % validate input trial number
        currenttrial = round(str2num(get(TrialUsed,'String')));
        currenttrial = max(1,currenttrial);
        currenttrial = min(length(triallist),currenttrial);
        set(TrialUsed,'String',num2str(currenttrial));
        % load dataset        
        filetrial = triallist{currenttrial};       
        originaldata = imreadtiffstack(filetrial);
        data = originaldata;
            if isempty(strfind(filetrial,'Image_Registration'))
                isRobert = 1;
            else
                isRobert = 0;
            end
            % correct for edge line issues
            if isRobert
                border = 0;
            else 
                border = 20;
            end
            data = NullImageBorders(border,border,0,0,data);        
        % reset main variables
        filtereddata = [];
        dffdata = [];
        dfffiltereddata = [];
        datadim = size(data);
        minframe = 1;
        maxframe = datadim(3);
%         startframe = minframe;
%         endframe = maxframe;
%         currentframe = 1;           
        updatevideo;
    end
TrialUsedTxt = uicontrol(DisplayFigure,'Style','Text');
set(TrialUsedTxt,'String','trial','Units','centimeters')
set(TrialUsedTxt,'Position',[5.5 4 2 .5])






%% setup indicators
 

% %% setup controls to pick cells
% cellcontours = [];
% 
% AddButton = uicontrol(DisplayFigure,'Style','PushButton');   
% set(AddButton,'String','pick cells','Units','centimeters');
% set(AddButton,'Callback',@AddButtonCallback);
% set(AddButton,'Position',[4 5 2 .5]);
%     function [] = AddButtonCallback(~,~)
%             % highlighter tool to identify blob
%             doneflag = 0;
%             BW = {};
%             while doneflag == 0
%                 roihandles = imfreehand(gca);
%                 cellmask = createMask(roihandles,ContourImg);
%                 if sum(cellmask(:)) == 0
%                     doneflag = 1;
%                 else
%                     BW{end+1} = cellmask;
%                 end
%             end
%             % identify the selected regions
%             if isempty(cellcontours)
%                 numcells = 0;
%             else
%                 numcells = size(cellcontours,3);
%             end                
%             for q = 1:length(BW)
%                 cellcontours(:,:,numcells+q) = (numcells+q)*(BW{q});
%             end            
%             customcmap = [[0 0 0];jet(size(cellcontours,3))];
%             contourmask = ind2rgb(1+max(cellcontours,[],3),customcmap);
%             hold on
%             DataImg = imshow(data(:,:,currentframe),[],'Parent',MovieAxis);
%             ContourImg = imshow(contourmask,'Parent',MovieAxis);
%             set(ContourImg,'AlphaData',  .2*(max(contourmask,[],3) > 0));
%             hold off
%     end





end