function [output] = gsprinter(FigureHandle,savefile)
    % enable execution block
    set(FigureHandle,'CloseRequestFcn','')
    set(FigureHandle,'Name','Printing, please wait')
    pause(.01)
    % determine filenames
    pngsavefile = savefile;
    pngsavefile((end-2:end)) = 'png';
    pdfsavefile = savefile;
    % print figures
    try
        export_fig(FigureHandle,pdfsavefile,'-q101','-painters','-nocrop')
        export_fig(FigureHandle,pngsavefile,'-q101','-painters','-nocrop')
        output = 'ghostscript';
    catch
        print(FigureHandle, '-dpdf', '-painters', pdfsavefile)
        print(FigureHandle, '-dpng', '-painters', pngsavefile)
        output = 'native';
    end
    % remove execution block
    set(FigureHandle,'CloseRequestFcn',@delete)
    set(FigureHandle,'Name','')
    pause(.01)    
end

