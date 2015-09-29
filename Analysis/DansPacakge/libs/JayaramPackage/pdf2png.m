function [output] = pdf2png(savefilepdf)
    savefilepdf = 'C:\Users\rpjb\Dropbox\Share (from jayaram)\02132013\plots\cell113.pdf';
    % validation of input file
    ext = savefilepdf((end-2):end);
    if ~strcmp(ext,'pdf') || ~exist(savefilepdf)
        disp('filename not a pdf')
    end
    % format output png filename
    savefilepng = savefilepdf;
    savefilepng((end-2):end) = 'png';
    savefilepdf = ['"' savefilepdf '"'];
    savefilepng = ['"' savefilepng '"'];
    % create output png
    evalstr = ['-dNOPAUSE -dBATCH -sDEVICE=pngalpha -r300 -sOutputFile=' savefilepng ' ' savefilepdf];
    ghostscript(evalstr);
end