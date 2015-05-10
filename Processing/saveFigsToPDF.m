function saveFigsToPDF (dirToLoad,fName)
disp('looking up dir')
    files = dir(dirToLoad);
   fprintf('Dir contains %d files\n',length(files)); 
   nameList = {} ;
   cnt = 1;
   
    for i=1:length(files)
       fileName = files(i,1).name;
       if strfind(fileName,'.fig')
%            if cnt > 1 break;end dev purposes
           h = hgload(strcat(dirToLoad,'/',fileName));
           title(fileName)
           pdfName = files(i,1).name(1:end-4);
           pdfName = strcat(pdfName,'.pdf');
           pdfName = strcat(dirToLoad,'/',pdfName,' ');
           %requires export_fig lib and ghostscript bin
           export_fig(h,pdfName,'-pdf') ;
           nameList{cnt} =   [pdfName ' '];
           close(h);
           cnt = cnt+1;
       end
    end
    close all;
    finaleName = strcat(fName,'.pdf');
    finaleName = [finaleName ' ']
    nameList = cell2mat(nameList);
    matrixName = mat2str(nameList);
    matrixName(regexp(matrixName,'['',;]'))=' ';
    matrixName = strcat(' ',matrixName(2:end-1));
    
    command = [ 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=', finaleName,  matrixName ]
    [status,cmdout] = system(command)

end