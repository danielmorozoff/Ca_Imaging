% modularcellplayer(folder,celllist, , , ,[startframe endframe])

function [] = modularcellplayer(folder,celllist,markedsets,trialstate,trialsetnames,whichframes,infodata,markerframe)
    %% validate inputs
    if isempty(celllist)
        disp('no cells selected')
        return
    end
    if ~exist('infodata','var')
        infodata = '';
    end
    numcells = length(celllist);
    cellselected = 1;
    
    %% setup figure
    DisplayFigure = figure;
    set(DisplayFigure,'units','centimeters','position',[5 2 8 5])
    set(DisplayFigure,'Name','Cell Controls','NumberTitle','off')
    
    % print button
    PrintButton = uicontrol(DisplayFigure,'Style','PushButton');
    set(PrintButton,'String','Save Cell','Units','centimeters');
    set(PrintButton,'Position',[.5 3 2 .5])
    set(PrintButton,'Callback',@printbuttonfcn)
    function printbuttonfcn(~,~)
        filename = fullfile(folder,'plots',[printprefix '_cell' num2str(cellselected) '.pdf']);
        gsprinter(CellFigure,filename);
    end
    printprefix = '';
    PrintPrefix = uicontrol(DisplayFigure,'Style','Edit');
    set(PrintPrefix,'String','','Units','centimeters');
    set(PrintPrefix,'Position',[.5 3.5 2 .5])
    set(PrintPrefix,'Callback',@printprefixfcn)
    function printprefixfcn(~,~)
        printprefix = get(PrintPrefix,'String');
    end
    
    % print all button
    PrintAllButton = uicontrol(DisplayFigure,'Style','PushButton');
    set(PrintAllButton,'String','Save All','Units','centimeters');
    set(PrintAllButton,'Position',[.5 2.5 2 .5])
    set(PrintAllButton,'Callback',@printallbuttonfcn)
    function printallbuttonfcn(~,~)
        for cellindex = 1:numcells
            % display cell
            set(Cellindex,'Value',cellindex);
            CellindexFcn;
            % print cell
            printfolder = fullfile(folder,'plots');
            if ~exist(printfolder)
                mkdir(printfolder)
            end            
            filename = fullfile(folder,'plots',[printprefix '_cell' num2str(cellselected) '.pdf']);
            gsprinter(CellFigure,filename);
        end
    end
    
    % text display
    InfoTxt = uicontrol(DisplayFigure,'Style','Text');
    set(InfoTxt,'String','Cell Number','Units','centimeters')
    set(InfoTxt,'HorizontalAlignment','left')
    set(InfoTxt,'Position',[.5 .5 7 1.5])
    set(InfoTxt,'FontSize',8)
    set(InfoTxt,'String',infodata)

    % create figure for plot display
    CellFigure = figure;
    %% cell controls    
    cellindex = 1;
    Cellindex = uicontrol(DisplayFigure,'Style','slider');
    set(Cellindex,'Units','centimeters','Position',[3 3 4 .5])
    set(Cellindex,'Min',1,'Max',numcells,'Value',1,'Callback',@CellindexFcn)
    set(Cellindex,'SliderStep',[1,5]/(numcells-1))
        function [] = CellindexFcn(~,~)
            cellindex = round(get(Cellindex,'Value'));
            cellindex = min(max(1,cellindex),numcells);
            set(CellindexEdit,'String',num2str(cellindex));
            cellselected = celllist(cellindex);
            updateoutput(cellselected)
        end
    CellindexTxt = uicontrol(DisplayFigure,'Style','Text');
    set(CellindexTxt,'String','cell order','Units','centimeters')
    set(CellindexTxt,'Position',[3 3.5 2 .5])
    CellindexEdit = uicontrol(DisplayFigure,'Style','Edit');
    set(CellindexEdit,'String',1,'Units','centimeters');
    set(CellindexEdit,'Position',[5 3.5 2 .5],'Callback',@CellindexEditFcn)
        function [] = CellindexEditFcn(~,~)
            cellindex = round(str2num(get(CellindexEdit,'String')));
            cellindex = min(max(1,cellindex),numcells);
            set(CellindexEdit,'String',num2str(cellindex));
            set(Cellindex,'Value',cellindex);
            cellselected = celllist(cellindex);
            updateoutput(cellselected)
        end

    %% output function    
    function updateoutput(cellselected)
        if sum(markedsets) == 0
            Plotcellsummary(folder,cellselected,[],[],0,CellFigure,markerframe);
        elseif sum(markedsets) == 1
            Plotcellsummary(folder,cellselected,trialstate,whichframes,0,CellFigure,markerframe);                    
        elseif sum(markedsets) > 1
            Plotcellsets(folder,cellselected,trialstate,trialsetnames,whichframes,0,CellFigure,markerframe);
        end
    end
end

