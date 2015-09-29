function pdf2pngbatch(folder)
    pdffiles = dir([folder filesep '*.pdf']);
    for i=1:length(pdffiles)
        filename = fullfile(folder,pdffiles(i).name);
        pdf2png(filename);
    end
end