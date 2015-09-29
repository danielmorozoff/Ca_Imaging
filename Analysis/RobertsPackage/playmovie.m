function playmovie(DisplayFigure)
file = '';
%% load the video
if isempty(file)
    file = Filemacpc('C:\Users\rpjb\Dropbox\project - geniculate newdata\GCaMP6s\02182013-01\mid-sevensucrosekcl\TSeries-02182013-1544-017\aligned.tif');
end
if ~exist(file)
    disp('File not found.')
    return
end

originaldata = imreadtiffstack(file);
data = originaldata;
filtereddata = [];
dffdata = [];
dfffiltereddata = [];
datadim = size(data);
dataintensityrange = [min(data(:)) max(data(:))];





%% initialize variables
% frames to play
minframe = 1;
maxframe = datadim(3);
startframe = minframe;
endframe = maxframe;
currentframe = 1;
frameperiod = .005;
% colormap
    
% break warning
stopflag = 0;

%% setup figure
% DisplayFigure = figure;
set(DisplayFigure,'units','centimeters','position',[5 5 24.5 17])

MovieAxis = axes;
set(MovieAxis,'units','centimeters','position', [8 .5 16 16])
DataImg = imshow(zeros(datadim(1),datadim(2)),[],'Parent',MovieAxis);
set(MovieAxis,'CLim',dataintensityrange)
%% create a contour mask
hold on
contourmask = zeros(datadim(1),datadim(2),3);
ContourImg = imshow(contourmask,'Parent',MovieAxis);
set(ContourImg,'AlphaData',  max(contourmask,[],3) > 0)
hold off

%% video controls
PlayButton = uicontrol(DisplayFigure,'Style','PushButton');
set(PlayButton,'String','play','Units','centimeters');
set(PlayButton,'Callback',@PlayFcn);
set(PlayButton,'Position',[.5 .5 2 .5]);
    function [] = PlayFcn(~,~)
        for currentframe=startframe:1:endframe
            set(CurrentFrameText,'String',currentframe);
            set(DataImg,'CData',data(:,:,currentframe));
            pause(frameperiod)
            if stopflag == 1
                stopflag = 0;
                return
            end
        end
    end
StopButton = uicontrol(DisplayFigure,'Style','PushButton');
set(StopButton,'String','stop','Units','centimeters');
set(StopButton,'Callback',@StopFcn);
set(StopButton,'Position',[.5 1 2 .5]);
    function [] = StopFcn(~,~)
        stopflag = 1;
    end

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
set(FrameRate,'String',num2str(1/frameperiod),'Units','centimeters')
set(FrameRate,'Callback',@FrameRateCallback)
set(FrameRate,'Position',[5.5 2 1 .5])
    function [] = FrameRateCallback(~,~)
        frameperiod = 1/str2num(get(FrameRate,'String'));
    end
FrameRateTxt = uicontrol(DisplayFigure,'Style','Text');
set(FrameRateTxt,'String','framerate','Units','centimeters')
set(FrameRateTxt,'Position',[3 2 2 .5])

StartingFrame = uicontrol(DisplayFigure,'Style','Edit');
set(StartingFrame,'String','1','Units','centimeters')
set(StartingFrame,'Callback',@StartingFrameCallback)
set(StartingFrame,'Position',[5.5 1 1 .5])
    function [] = StartingFrameCallback(~,~)
        startframe = max(1,str2num(get(StartingFrame,'String')));
    end
StartingFrameTxt = uicontrol(DisplayFigure,'Style','Text');
set(StartingFrameTxt,'String','start','Units','centimeters')
set(StartingFrameTxt,'Position',[3 1 2 .5])

EndingFrame = uicontrol(DisplayFigure,'Style','Edit');
set(EndingFrame,'String',num2str(datadim(3)),'Units','centimeters')
set(EndingFrame,'Callback',@EndingFrameCallback)
set(EndingFrame,'Position',[5.5 1.5 1 .5])
    function [] = EndingFrameCallback(~,~)
        endframe = min(datadim(3),str2num(get(EndingFrame,'String')));
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
set(KeyboardButton,'Position',[3 3 2 .5]);
    % function for kuwahara cells
    function [] = KeyboardFcn(~,~)
        keyboard
    end
%% processed video to display
usefilterindex = 0;
usedffindex = 0;

UseFilterButton = uicontrol(DisplayFigure,'Style','checkbox','Value',0);
set(UseFilterButton,'String','show filter','Units','centimeters');
set(UseFilterButton,'Callback',@UseFilter)
set(UseFilterButton,'Position',[3 4.5 2 .5])
set(UseFilterButton,'Enable','off')
    function [] = UseFilter(~,~)
        usefilterindex = get(UseFilterButton,'Value');
        if usefilterindex
            usefilterindex = 1;
        elseif ~usefilterindex
            usefilterindex = 0;
        end
        updatevideo;
    end
UseDFFButton = uicontrol(DisplayFigure,'Style','checkbox','Value',0);
set(UseDFFButton,'String','show dff','Units','centimeters');
set(UseDFFButton,'Callback',@UseDFF)
set(UseDFFButton,'Position',[3 4 2 .5])
set(UseDFFButton,'Enable','off')
    function [] = UseDFF(~,~)
        usedffindex = get(UseDFFButton,'Value');
        if usedffindex
            usedffindex = 1;
        elseif ~usedffindex
            usedffindex = 0;
        end
        updatevideo;
    end

    % subfunction that checks the dff / filter indices and updates the data
    % range
    function [] = updatevideo()
        if usedffindex && usefilterindex
            data = dfffiltereddata;
            data(data>5) = 5;
        elseif usedffindex && ~usefilterindex
            data = dffdata;
        elseif ~usedffindex && usefilterindex
            data = filtereddata;
        elseif ~usedffindex && ~usefilterindex
            data = originaldata;
        end
        dataintensityrange = [min(data(:)) max(data(:))];
        set(MovieAxis,'CLim',dataintensityrange)
    end  

FilterButton = uicontrol(DisplayFigure,'Style','PushButton');
set(FilterButton,'String','calc filter','Units','centimeters');
set(FilterButton,'Callback',@FilterFcn);
set(FilterButton,'Position',[.5 3 2 .5]);
    function [] = FilterFcn(~,~)
        [filtereddata] = KalmanPadStackFilter(originaldata);
        set(UseFilterButton,'Enable','on')
        set(FilterButton,'Enable','off')
    end
DFFButton = uicontrol(DisplayFigure,'Style','PushButton');
set(DFFButton,'String','calc dff','Units','centimeters');
set(DFFButton,'Callback',@DFFFcn);
set(DFFButton,'Position',[.5 3.5 2 .5]);
    function [] = DFFFcn(~,~)
        % the first ten frames look bad, so repeat the first ten frames
        d = size(originaldata);
        mediandata = double( repmat(median(data,3),[1 1 d(3)]) );
        dffdata = (double(originaldata) - mediandata) ./ (mediandata);
        set(UseDFFButton,'Enable','on')
        set(DFFButton,'Enable','off')
        set(DFFFilterButton,'Enable','on');

    end

DFFFilterButton = uicontrol(DisplayFigure,'Style','PushButton');
set(DFFFilterButton,'String','calc dfffilt','Units','centimeters');
set(DFFFilterButton,'Callback',@DFFFilterFcn);
set(DFFFilterButton,'Position',[.5 4 2 .5]);
set(DFFFilterButton,'Enable','off');
    function [] = DFFFilterFcn(~,~)
        % the first ten frames look bad, so repeat the first ten frames
        d = size(dffdata);
        r = zeros(d(1),d(2),d(3)+10);
        r(:,:,1:10) = double(dffdata(:,:,1:10));
        r(:,:,11:end) = double(dffdata);
        % run kalman stack and truncate first ten frames
        dfffiltereddata = KalmanStackFilter(r,.8,.05);
        dfffiltereddata = dffdata(:,:,11:end);
        set(UseFilterButton,'Enable','on')
        set(DFFFilterButton,'Enable','off')
    end




%% setup indicators
 

%% setup controls to pick cells
cellcontours = [];

AddButton = uicontrol(DisplayFigure,'Style','PushButton');   
set(AddButton,'String','pick cells','Units','centimeters');
set(AddButton,'Callback',@AddButtonCallback);
set(AddButton,'Position',[4 5 2 .5]);
    function [] = AddButtonCallback(~,~)
            % highlighter tool to identify blob
            doneflag = 0;
            BW = {};
            while doneflag == 0
                roihandles = imfreehand(gca);
                cellmask = createMask(roihandles,ContourImg);
                if sum(cellmask(:)) == 0
                    doneflag = 1;
                else
                    BW{end+1} = cellmask;
                end
            end
            % identify the selected regions
            if isempty(cellcontours)
                numcells = 0;
            else
                numcells = size(cellcontours,3);
            end                
            for q = 1:length(BW)
                cellcontours(:,:,numcells+q) = (numcells+q)*(BW{q});
            end            
            customcmap = [[0 0 0];jet(size(cellcontours,3))];
            contourmask = ind2rgb(1+max(cellcontours,[],3),customcmap);
            hold on
            DataImg = imshow(data(:,:,currentframe),[],'Parent',MovieAxis);
            ContourImg = imshow(contourmask,'Parent',MovieAxis);
            set(ContourImg,'AlphaData',  .2*(max(contourmask,[],3) > 0));
            hold off
    end





end