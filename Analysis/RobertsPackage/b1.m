% B1 Visual interface to identify cells using spatial or temporal ICA.
%
% b1 is a standalone program that is used to load imaging data and provide
% a visual interface for ICA analysis. a user first loads a directory 
% containing imaging data and an average image. a user then selects a box
% size that can loosely enclose a single cell body.
%
% a user can then run an automated "cell picker" which uses the expected 
% circular nature of cells to identify approximate cell locations. spatial
% ICA is then performed on the box of pixels neighboring the cell for the
% first available trial dataset. identified cell are then overlayed in 
% green, and regions without an identified cell are shown in red.
%
% a user can also manually select cells by clicking on various locations on
% the displayed image.
%
% lastly, spatial ICA only works well on images with bright cells, which
% correspond to depolarized neurons. this might only happen in a subset of
% trial datasets. so, a user can search other trials, which looks
% sequentially through datasets to find remaining neurons.
%
% cell data can then be saved on "save cells" for future loading. or when
% completed, the user can select "save traces" which then extracts traces
% based on the calculated spatial filters. 
%
% type: function
%
% inputs:
%   debug - can be either a number or a filepath, which indicates what cell
%     data to load
%
% outputs: none
% 
% dependencies on custom functions:
%   Filemacpc
%   FindImagingFiles
%   NullImageBorders
%   imreadtiffstack
%   cellmap
%   temporalica


function [] = b1(debug)
%% declare variables that have full scope

% switch for debug
if nargin == 0
    debug = 0;
end

% variables for cell parameters to be saved
celltraces = {};
cellxrad = [];
cellyrad = [];
cellspatialfilters = {};
cellx = [];
celly = [];
cellcentroidx = [];
cellcentroidy = [];
cellarea = [];
celltrial = [];

% variables for manually selected cells
currentx = [];
currenty = [];
currentxrad = [];
currentyrad = [];
badx = [];
bady = [];
baddone = [];

% set the size of cell boxes
% choose the box size for analysis
xrad = 10; yrad = 10;

% variables for images
avgimg = zeros(5);
[imgheight imgwidth] = size(avgimg);
data = [];
frameperiod = 0;

% variables for loading files
avgfile = [];
filepath = ''; 
filelist = {};
filetrial = [];
currenttrial = 1;
tracefile = '';
isRobert = '';

% horizontal border to trim Jayaram's datasets
border = 0;

%% create image display

DisplayFigure = figure;
set(DisplayFigure,'units','centimeters','position',[5 8.5 24.5 17])
set(DisplayFigure,'Name','Cell Picker','NumberTitle','off')
set(DisplayFigure,'Toolbar','figure')

% image display for large imagekey
ImageAxis = axes;
set(ImageAxis,'units','centimeters','position', [8 .5 16 16])
MainImg = imshow(zeros(500),[],'Parent',ImageAxis);
set(ImageAxis,'NextPlot','add');
OverlayImg = imshow(zeros([imgheight imgwidth]),[],'Parent',ImageAxis);
set(OverlayImg,'AlphaData',.25);
set(ImageAxis,'NextPlot','replacechildren');

% image display for single cell
CellImgAxis = axes;
set(CellImgAxis,'units','centimeters','position', [.5 9.5 7 7])
imshow(zeros(5),[],'Parent',CellImgAxis);

%% create status display
% status update
StatusText = uicontrol(DisplayFigure,'Style','text','HorizontalAlignment','left');
set(StatusText,'String','','Units','centimeters');
set(StatusText,'Position',[.5 7 7 0.5]);

% file loaded
FileLoadText = uicontrol(DisplayFigure,'Style','text','HorizontalAlignment','left');
set(FileLoadText,'String','','Units','centimeters');
set(FileLoadText,'Position',[.5 8 7 1])

%% load data
FileButton = uicontrol(DisplayFigure,'Style','PushButton');
set(FileButton,'String','load avg file','Units','centimeters');
set(FileButton,'Callback',@FileFcn);
set(FileButton,'Position',[.5 6 2 .5]);
    % function for loading files
    function [] = FileFcn(~,~)
        ButtonManager('off');              
        set(StatusText,'String','Looking for files');
        pause(.1)
        % interface for finding file
        if debug>0
            if isstr(debug) && exist(debug,'dir')
                filepath = debug;
            elseif debug == 1
                filepath = Filemacpc('C:\Users\rpjb\Dropbox\Share (from jayaram)\02132013\');
            elseif debug == 2
                filepath = Filemacpc('C:\Users\rpjb\Dropbox\project - geniculate newdata\GCaMP6s\02182013-01');
            elseif debug == 3
                filepath = Filemacpc('C:\Users\rpjb\Dropbox\project - geniculate newdata\GCaMP6s\02232013-01');                
            elseif debug == 4
                filepath = Filemacpc('C:\Users\rpjb\Dropbox\project - geniculate newdata\GCaMP6s\03042013-01');
            elseif debug == 5
                filepath = Filemacpc('C:\Users\rpjb\Dropbox\Share (from jayaram)\doughnut cells\2013_05_22');
            end
        else            
            filepath = uigetdir(filepath,'Select the directory containing images');    
            if filepath==0
                set(StatusText,'String','No directory chosen')
                set(FileButton,'Enable','on')
                filepath = '';
                return
            end
        end
        [avgfile,filelist,~,tracefile,~,isRobert] = FindImagingFiles(filepath);
        if isempty(avgfile) || isempty(filelist)
            set(FileButton,'Enable','on')
            filepath = '';
            return
        end
        % display image in figure
        avgimg = mat2gray(imread(avgfile)); % rescales to [0 1]
        [imgheight imgwidth] = size(avgimg);
        mag = 200;
        set(ImageAxis,'XLim',[0 imgwidth],'YLim',[0 imgheight])        
        set(ImageAxis,'Clim',[0 1])
        set(MainImg,'CData',avgimg)
        % load trial
        LoadDataFile;

        set(StatusText,'String','Files loaded')
        ButtonManager
        set(FileButton,'Visible','off')
        
        % framerate detection
        if ~isRobert % for scanimage files
            x =imfinfo(filelist{1});
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
            [samplefolder,~,~] = fileparts(filelist{1});
            frameperiod = getframeperiod(samplefolder);
        end
        
    end

%% sliders for image contrast
LowImgButton = uicontrol(DisplayFigure,'Style','slider');
set(LowImgButton,'Units','centimeters','Position',[5.5 5 2 .5])
set(LowImgButton,'Min',0,'Max',100,'Value',0,'Callback',@UpdateImgFcn)
set(LowImgButton,'Enable','off')
HighImgButton = uicontrol(DisplayFigure,'Style','slider');
set(HighImgButton,'Units','centimeters','Position',[5.5 5.5 2 .5])
set(HighImgButton,'Min',0,'Max',100,'Value',100,'Callback',@UpdateImgFcn)
set(HighImgButton,'Enable','off')
    function [] = UpdateImgFcn(~,~)
        % get and validate bounds
        n = round(get(LowImgButton,'Value'));
        m = round(get(HighImgButton,'Value'));        
        n = min(max(0,n),m); 
        m = max(min(100,m),n);
        set(LowImgButton,'Value',n);
        set(HighImgButton,'Value',m);
        % adjust image lookup table
        set(ImageAxis,'CLim',[n,m]/100);
    end

%% button for box size
BoxButton = uicontrol(DisplayFigure,'Style','PushButton');
set(BoxButton,'String','box size','Units','centimeters');
set(BoxButton,'Callback',@BoxFcn);
set(BoxButton,'Position',[3 5.5 2 .5]);
    % function for choosing button size
    function [] = BoxFcn(~,~)
        ButtonManager
        set(StatusText,'String','Click and drag rectangle around a cell')
        rect = getrect(ImageAxis);
        xrad = round(rect(3)/2);
        yrad = round(rect(4)/2);
        set(StatusText,'String',[num2str(xrad) 'by ' num2str(yrad) ' box chosen'])
        ButtonManager
    end

%% button for using kuwahara to select cells
KuwaharaButton = uicontrol(DisplayFigure,'Style','PushButton');
set(KuwaharaButton,'String','autopick cells','Units','centimeters');
set(KuwaharaButton,'Callback',@KuwaharaFcn);
set(KuwaharaButton,'Position',[.5 5.5 2 .5]);
    % function for kuwahara cells
    function [] = KuwaharaFcn(hObject,eventdata)
        ButtonManager;      
        edgetype = 'canny';
        set(StatusText,'String','Analyzing avgimg for cells')
        kimg = kuwahara(double(avgimg),9);
    
        edgeimg = edge(kimg,edgetype);
        % zero out borders
        edgeimg = NullImageBorders(border,border,0,0,edgeimg);
        % fill holes and remove non cells
        fillimg = imfill(edgeimg,'holes');
        cellimg = bwmorph(fillimg,'open');
        edgeimg = edge(cellimg,'log');
        % calculate an overlay image and display
        overlayimg = max( double(avgimg), double(max(avgimg(:)))*double(edgeimg));
        set(OverlayImg,'CData',overlayimg);
        set(StatusText,'String','Cells autopicked')        

        % label each putative cell
        labelimg = bwlabel(cellimg);
        numcells = max(labelimg(:));

        % identify the x and y coordinates for each cell
        stats = regionprops(cellimg,'Centroid');
        statsloc = [stats.Centroid];
        statsloc = reshape(statsloc,[2 length(statsloc)/2])';
        currentx = [currentx; round(statsloc(:,1))];
        currenty = [currenty; round(statsloc(:,2))];
        UpdateList
        ButtonManager;      
        
    end

%% button for manually selecting cells
ManuallyPickButton = uicontrol(DisplayFigure,'Style','PushButton');
set(ManuallyPickButton,'String','pick cells','Units','centimeters');
set(ManuallyPickButton,'Callback',@ManuallyPickFcn);
set(ManuallyPickButton,'Position',[.5 5 2 .5]);

    % function for manually selecting cells
    function ManuallyPickFcn(~,~)
        ButtonManager;
        set(StatusText,'String','Manually selecting cells')
        currenttrial = 1;
        [newcurrentx,newcurrenty] = getpts(ImageAxis);
        currentx = [currentx; round(newcurrentx)]; currenty = [currenty; round(newcurrenty)];
        % remove points selected out of the image window
        badindices = (currentx<1) | (currentx>imgwidth) | (currenty<1) | (currenty>imgheight);
        currentx(badindices) = [];
        currenty(badindices) = [];
        set(StatusText,'String','Finished manually selecting cells')

        savefile = fullfile(filepath,'autosavepts.mat');
        save(savefile,'currentx','currenty');
        
        % uncomment this line if you like to perform ICA immediately after
        % selecting cells
%         ICAFcn;
        UpdateList;
        ButtonManager;
    end
%% button for performing slower ICA with graphical output
useslowica = 0;
UseSlowICA = uicontrol(DisplayFigure,'Style','Checkbox','Value',useslowica);
set(UseSlowICA,'Units','centimeters','String','show ica')
set(UseSlowICA,'Position',[5.5 3 2 .5]);
set(UseSlowICA,'Callback',@UseSlowICAfcn)
    function [] = UseSlowICAfcn(~,~)
        if useslowica
            useslowica = 0;
        else
            useslowica = 1;
        end
    end

%% button for spawning movie controls
MovieButton = uicontrol(DisplayFigure,'Style','PushButton');
set(MovieButton,'String','movie controls','Units','centimeters');
set(MovieButton,'Callback',@MovieFcn);
set(MovieButton,'Position',[3 5 2 .5]);

% function for spawning movie
    function MovieFcn(~,~)
        modularmovieplayer(ImageAxis,MainImg,data,filelist,frameperiod)
        UpdateList
    end

%% button for removing cells
RemoveButton = uicontrol(DisplayFigure,'Style','PushButton');
set(RemoveButton,'String','remove cells','Units','centimeters');
set(RemoveButton,'Callback',@RemoveFcn);
set(RemoveButton,'Position',[.5 4.5 2 .5]);

% function for manually selecting cells to be removed
    function RemoveFcn(~,~)
        ButtonManager;
        set(StatusText,'String','Manually selecting cells')
        [rx,ry] = getpts(ImageAxis);
        rx = round(rx); ry = round(ry);
        % screen for points selected out of the image window
        badindices = (rx<1) | (rx>imgwidth) | (ry<1) | (ry>imgheight);
        rx(badindices) = [];
        ry(badindices) = [];

        
        removecelllist = [];
        % for each selected cell
        
        for i=1:length(rx)
            [~, removecelllist(i)] = min(sum([(cellcentroidx-rx(i)).^2;(cellcentroidy-ry(i)).^2],1));
        end
        % clear main variables for each cell
        cellxrad(removecelllist) = [];
        cellyrad(removecelllist) = [];
        celltraces(removecelllist) = [];
        cellspatialfilters(removecelllist) = [];
        cellx(removecelllist) = [];
        celly(removecelllist) = [];
        cellcentroidx(removecelllist) = [];
        cellcentroidy(removecelllist) = [];
        cellarea(removecelllist) = [];
        celltrial(removecelllist) = [];
        UpdateList;
                
        set(StatusText,'String','Finished manually selecting cells')
        ButtonManager;
    end
    
%% button for running ICA
ICAButton = uicontrol(DisplayFigure,'Style','PushButton');
set(ICAButton,'String','run ica','Units','centimeters');
set(ICAButton,'Callback',@ICAFcn);
set(ICAButton,'Position',[.5 3.5 2 .5]);

    % function for ica button
    function [] = ICAFcn(~,~)
        ButtonManager;      
        set(StatusText,'String','Running ICA')
        FindCells;
        set(StatusText,'String','ICA complete')
        UpdateList;        
        ButtonManager;      
    end

%% button for removing bad points
RemoveBadButton = uicontrol(DisplayFigure,'Style','PushButton');
set(RemoveBadButton,'String','remove bad pts','Units','centimeters');
set(RemoveBadButton,'Callback',@RemoveBadFcn);
set(RemoveBadButton,'Position',[3 3.5 2 .5]);

    % function for removing bad points
    function [] = RemoveBadFcn(~,~)
        ButtonManager;      
        set(StatusText,'String','Removing Bad points')
        badx = []; bady = []; baddone = [];
        set(StatusText,'String','')
        UpdateList;        
        ButtonManager;      
    end

%% button for iterative ICA
IterICAButton = uicontrol(DisplayFigure,'Style','PushButton');
set(IterICAButton,'String','search trials','Units','centimeters');
set(IterICAButton,'Callback',@IterICAFcn);
set(IterICAButton,'Position',[.5 3 2 .5]);
    % function for iterativeICA
    function [] = IterICAFcn(~,~)
        currenttrial = 1;
        % loop through
        while (currenttrial < length(filelist)) && (~isempty(badx))
            if sum(baddone)==length(baddone)
                return
            end
            % load points that didn't work prior, but not already scanned
            % for all trials
            currentx = badx(baddone==0);
            currenty = bady(baddone==0);
            % autosave points
            savefile = fullfile(filepath,'autosavepts.mat');
            save(savefile,'currentx','currenty');
            % load next trial
            currenttrial = currenttrial+1;
            filetrial = filelist{currenttrial};

            LoadDataFile;            
            set(StatusText,'String',['Loaded Trial ' num2str(currenttrial) ' of ' num2str(length(filelist))]);
            pause(.1)
            % run ICA
            FindCells
            UpdateList;            

        end
        
        if isempty(badx)
            baddone = [];
            set(StatusText,'String','All cells found')
        else
            % mark all bad cells so that they won't be analyzed again
            baddone = ones(length(badx),1);
            set(StatusText,'String','All available trials used')   
        end
    end

%% button for manual keyboard
ManualButton = uicontrol(DisplayFigure,'Style','PushButton');
set(ManualButton,'String','keyboard','Units','centimeters');
set(ManualButton,'Callback',@ManualFcn);
set(ManualButton,'Position',[3 3 2 .5]);
    % function for popping up keyboard
    function [] = ManualFcn(~,~)
        keyboard
    end

%% button for loading data
LoadButton = uicontrol(DisplayFigure,'Style','PushButton');
set(LoadButton,'String','load cells','Units','centimeters');
set(LoadButton,'Callback',@LoadFcn);
set(LoadButton,'Position',[.5 2 2 .5]);
    function LoadFcn(~,~)
        ButtonManager;
        set(StatusText,'String','Loading data');
        pause(.1)
        savefile = fullfile(filepath,'celldata.mat');
        if exist(savefile) > 0        
            S = load(savefile);
            cellxrad = S.cellxrad;
            cellyrad = S.cellyrad;
            xrad = S.xrad; yrad = S.yrad;
            celltraces = S.celltraces;
            cellspatialfilters = S.cellspatialfilters;
            cellx = S.cellx;
            celly = S.celly;
            cellcentroidx = S.cellcentroidx;
            cellcentroidy = S.cellcentroidy;
            cellarea = S.cellarea;
            celltrial = S.celltrial;                                
            if size(avgimg,1) == 5
                set(StatusText,'String','Load avgimg first')
            else
                UpdateList()
                set(StatusText,'String','Data loaded')
            end
            % load autosave points
            loadfile = fullfile(filepath,'autosavepts.mat');
            if exist(loadfile)>0
                pts = load(loadfile);
                currentx = pts.currentx;
                currenty = pts.currenty;
            end
            if ~isempty(currentx)
                set(StatusText,'String','Please run ica as there are points to be analyzed.');
            end            
        else % case when celldata.mat doesn't exist
            set(StatusText,'String','celldata file does not exist');
        end
        ButtonManager;        
    end

%% button for saving cell contour data
SaveButton = uicontrol(DisplayFigure,'Style','PushButton');
set(SaveButton,'String','save cells','Units','centimeters');
set(SaveButton,'Callback',@SaveFcn);
set(SaveButton,'Position',[.5 1.5 2 .5]);
    function SaveFcn(~,~)
        ButtonManager;
        savefile = fullfile(filepath,'celldata.mat');
        save(savefile,'celltraces','cellspatialfilters','cellx','celly','cellcentroidx','cellcentroidy','cellarea','celltrial','xrad','yrad','cellxrad','cellyrad')
        set(StatusText,'String','Data saved')
        ButtonManager;
    end

%% button for saving trace data
SaveTraceButton = uicontrol(DisplayFigure,'Style','PushButton');
set(SaveTraceButton,'String','save traces','Units','centimeters');
set(SaveTraceButton,'Callback',@SaveTraceFcn);
set(SaveTraceButton,'Position',[3 1.5 2 .5]);
    function SaveTraceFcn(~,~)
        ButtonManager;
        savefile = fullfile(filepath,'celltraces.mat');          
        % open each trial
        for t=1:length(filelist)
            set(StatusText,'String',['Extracting traces from trial ' num2str(t) ' of ' num2str(length(filelist))])
            pause(.1)
            data = imreadtiffstack(filelist{t});            
            for i=1:length(cellx)
                % identify the imagestack based on provided coordinates
                volumestack = cellmap(cellx(i),celly(i),cellxrad(i),cellyrad(i),data);                
                maskstack = repmat(cellspatialfilters{i},[1 1 size(volumestack,3)]);
                filtstack = double(volumestack) .* double(maskstack);
                timecourse = squeeze(sum(sum(filtstack,1),2));
                cellsig{i,t} = timecourse;
            end
        end
        save(savefile,'cellsig','cellspatialfilters','cellx','celly','cellcentroidx','cellcentroidy','cellarea','xrad','yrad','cellxrad','cellyrad')
        set(StatusText,'String','Data saved')
        ButtonManager;
    end

%% Button for making trial rasters

PlotTrialButton = uicontrol(DisplayFigure,'Style','PushButton');
set(PlotTrialButton,'String','plot trials','Units','centimeters');
set(PlotTrialButton,'Callback',@PlotTrialFcn);
set(PlotTrialButton,'Position',[3 .5 2 .5]);
    function PlotTrialFcn(~,~)
        ButtonManager;
        set(StatusText,'String','Generating plots....')
        pause(.1)
        for t = 1:length(filelist)
            [TrialFigure,RasterFigure] = Plottrial(filepath,t);
            close(TrialFigure)
            close(RasterFigure)
        end
        set(StatusText,'String','Data saved')
        ButtonManager;
    end

%% Button for making individual cell trace summaries
PlotCellButton = uicontrol(DisplayFigure,'Style','PushButton');
set(PlotCellButton,'String','plot cells','Units','centimeters');
set(PlotCellButton,'Callback',@PlotCellFcn);
set(PlotCellButton,'Position',[.5 .5 2 .5]);
    function PlotCellFcn(~,~)
        ButtonManager;
        set(StatusText,'String','Generating plots....')
        pause(.1)
        if exist(tracefile)>0
            localdata = load(tracefile);
            numcells = size(localdata.cellsig,1);
            for i=1:numcells
                set(StatusText,'String',['Generating plot of cell: ' num2str(i) ' of ' num2str(numcells)])
                PlotFigure = Plotcellsummary(filepath,i);
                close(PlotFigure)
            end
        end        
        set(StatusText,'String','Data saved')
        ButtonManager;
    end

%% maintenance function for buttons
    % loads datafile
    function LoadDataFile()
        set(FileLoadText,'String',avgfile)
        filetrial = filelist{currenttrial};       
        data = imreadtiffstack(filetrial);
        % correct for edge line issues
        if isRobert
            border = 0;
        else 
            border = 20;
        end
        data = NullImageBorders(border,border,0,0,data);
    end
 
    % blocks buttons while functions execute
    function ButtonManager(varargin)
        buttonlist = {PlotCellButton,PlotTrialButton,BoxButton,IterICAButton,FileButton,KuwaharaButton,ManuallyPickButton,RemoveButton,ICAButton,ManualButton,LoadButton,SaveButton,SaveTraceButton,MovieButton,RemoveBadButton,LowImgButton,HighImgButton};
        if isempty(varargin)
            for i=1:length(buttonlist)            
                if strcmp(get(buttonlist{i},'Enable'),'off')
                    set(buttonlist{i},'Enable','on')
                else
                    set(buttonlist{i},'Enable','off')
                end
            end        
        elseif strcmp(varargin,'off')            
            for i=1:length(buttonlist)            
                set(buttonlist{i},'Enable','off')
            end            
        elseif strcmp(varargin,'on')
            for i=1:length(buttonlist)            
                set(buttonlist{i},'Enable','on')
            end            
        end
    end



    % updates list of points 
    % assigns red crosses for unsuccessful points
    function [] = UpdateList()
       % calculate overlay image -- value of zero means pass image through
       contourmask = zeros(imgwidth,imgheight);
       % for each good cell set to value 2
       for i=1:size(celltraces,2)
            xpt = cellx(i); ypt = celly(i);
            % convert from local spatial coordinates to global
            boxx = max(xpt-cellxrad(i),1):min(xpt+cellxrad(i),imgwidth);
            boxy = max(ypt-cellyrad(i),1):min(ypt+cellyrad(i),imgheight);
            contourmask(boxy,boxx) = contourmask(boxy,boxx) + 2*(cellspatialfilters{i} > 0);  
       end              
       % for each bad cell set cross to value 1
       for i=1:length(badx)
           contourmask(bady(i)+[-round(yrad/4):round(yrad/4)],badx(i)) = 1;
           contourmask(bady(i), badx(i)+[-round(xrad/4):round(xrad/4)]) = 1;
       end
       % for each unprocessed cell set to value 3
       for i=1:size(currentx)
           contourmask(currenty(i)+[-round(yrad/4):round(yrad/4)],currentx(i)) = 3;
           contourmask(currenty(i), currentx(i)+[-round(xrad/4):round(xrad/4)]) = 3;
       end

       % create lookup table
       customcmap = [[0 0 0];[1 0 0];[0 1 0];[0 0 1];[1 1 0]];
       contourmask = ind2rgb(1+contourmask,customcmap);
       % display image    
       set(ImageAxis,'CLim',[0 1])
       set(MainImg,'CData',avgimg);
       set(OverlayImg,'CData',contourmask)
       set(OverlayImg,'AlphaData',.25);
    end

% runs the ica
    % create a containing figure for slowica components
    if ~exist('plotfigure','var')
        disp('ICA figure did not exist')
        plotfigure = figure;                        
        figure(DisplayFigure); %bring display figure to front
    end
    function [] = FindCells()
        % keep a tab on points that don't have detected cells
        badpoints = zeros([1 length(currentx)]);   
        % loop through all unselected points
        for i=1:length(currentx)
            set(StatusText,'String',['Current Trial ' num2str(currenttrial) ' of ' num2str(length(filelist)) ', ICA cell ' num2str(i) ' of ' num2str(length(currentx))])
            pause(.1)
            % identify the imagestack based on provided coordinates
            if isRobert
                volumestack = cellmap(currentx(i),currenty(i),xrad,yrad,data);                
                % perform kalman smoothing
                volumestack = KalmanStackFilter(double(volumestack),.8,.05);
%                 volumestack = cellmap(currentx(i),currenty(i),xrad,yrad,data);                                
            else % perform spatial median filter for jayaram's line corrected dataset
                medianwindow = 4; borderpad = round(medianwindow/2);
                [volumestack,xrange,yrange] = cellmap(currentx(i),currenty(i),xrad+borderpad,yrad+borderpad,data);                
                volumestack = SpatialMedianFilter(volumestack,medianwindow);
                % remap the padded volume stack back to the correct region
                xmin = max(currentx(i)-xrad,1);
                xmax = min(currentx(i)+xrad,imgheight);
                ymin = max(currenty(i)-yrad,1);
                ymax = min(currenty(i)+yrad,imgwidth);    
                volumestack = volumestack(yrange<=ymax & yrange>=ymin,xrange<=xmax & xrange>=xmin,:);                
                volumestack = KalmanStackFilter(double(volumestack),.8,.05);
            end           
            % calculate trace and image for a given image stack
%            [outputtrace,outputimage] = spatialica(volumestack);
                if ~useslowica
                    [outputtrace,outputimage] = temporalica(volumestack);
                else
                    [outputtrace,outputimage] = Slowtemporalica(volumestack,plotfigure);                    
                end
                
            if ~isempty(outputimage) %check if outputimage is empty
                imshow(sum(outputimage,3),[],'Parent',CellImgAxis);
                % keep a counter for number of cells in list
                counter = size(celltraces,2) + 1;          
                startcounter = counter;
                % mark cells that border the ICA window edges or have too small
                % area
                areathreshold = .25*xrad*yrad;                                
                % this was code for handling multiple cells detected in a FOV                
                for j=1:size(outputtrace,2)
                    space = outputimage(:,:,j);                    
                    areaofspace = sum(sum(space>0));
                    if (areaofspace<areathreshold) || (sum(space(1,:))+sum(space(end,:))+sum(space(:,1))+sum(space(:,end))>0)
                        % don't count cells that are below area size or on
                        % borders
                        continue
                    else
                        % allocate to cell structure to save parameters
                        cellxrad(counter) = xrad;
                        cellyrad(counter) = yrad;
                        celltraces{counter} = outputtrace(:,j);
                        cellspatialfilters{counter} = outputimage(:,:,j);
                        cellx(counter) = currentx(i);
                        celly(counter) = currenty(i);                
                        stats = regionprops(outputimage(:,:,j)>0,'Area','Centroid');
                        if length(stats)>1 % get the largest particle
                            stats = stats(find(max([stats.Area])));
                        end
                        cellcentroidx(counter) = max(currentx(i)-xrad,1) + round(stats(1).Centroid(1)) - 1 ;
                        cellcentroidy(counter) = max(currenty(i)-yrad,1) + round(stats(1).Centroid(2)) - 1 ;
                        cellarea(counter) = stats.Area(1);
                        celltrial(counter) = currenttrial;
                        counter = counter+1;
                    end
                end
                % mark as a bad point if no cells were added
                if startcounter == counter
                    badpoints(i) = 1;
                end
            else
                imshow(zeros(5),[],'Parent',CellImgAxis);
                badpoints(i) = 1;
            end
        end
        % add recently found bad points to list, and prune duplicates
            badx = [badx; currentx(logical(badpoints))]; bady = [bady; currenty(logical(badpoints))];
            baddone = [baddone; zeros(sum(badpoints),1)];                
            [uniquebadpts, uniquebadind, ~] = unique([badx bady],'rows','first');
            if isempty(uniquebadpts)
                badx = []; bady = [];
                baddone = [];                
            else
                badx = uniquebadpts(:,1); bady = uniquebadpts(:,2);                                
                baddone = baddone(uniquebadind);
            end

            % find bad points that were not addressed
            goodx = currentx(~logical(badpoints)); goody = currenty(~logical(badpoints));           
            if ~isempty(goodx) && ~isempty(badx)
                [uniquebadpts, uniquebadind] = setdiff([badx bady],[goodx goody],'rows');
            end
            % update bad points 
            if isempty(uniquebadpts) 
                badx = []; bady = []; baddone = [];
            else
                badx = uniquebadpts(:,1); bady = uniquebadpts(:,2);                
                baddone = baddone(uniquebadind);
            end
            currentx = []; currenty = [];        
    end
% initial settings of buttons, only allow loading of file
ButtonManager;
set(FileButton,'Enable','on')

% initialize without requiring a click
if ~isempty(debug)
    FileFcn;
end

end %end main function