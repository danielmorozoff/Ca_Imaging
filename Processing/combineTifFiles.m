%Combine tiff
function [tif_files,str]=combineTifFiles(folderDir,outputFileName)
    tif_files = dir(fullfile(strcat(folderDir,'*.tif')));
    fprintf('Found %d tiff files\n',size(tif_files));
     %Create new dir
     mkdir(outputFileName);
     i=1;
    while length(tif_files)>=1
%        disp(tif_files(i,1));
       file = tif_files(i,1);
       singleFileName = file.name(1:end-8);
        %find all files with this name
%         indx = findIndicesOfObject(singleFileName,tif_files);
        indx = [1 2 3 4];
        %Combine all files at the specified indx
        outImage= [];
        tifFile=[];
        for ind=1:length(indx)
            fname = strcat(folderDir,tif_files(indx(ind)).name);
            info = imfinfo(fname);
            num_images = numel(info);
            if(i==1)
               %add header info
               tifFile = Tiff(strcat(outputFileName,'/',singleFileName,'.tif'),'w');
%                propToSet = info.ImageDescription;
               
               setTag(tifFile,headerInfo);
%                structureTag.ImageDescription = info.ImageDescription;
             
            end
            for k = 1:num_images
                [appendImg,map] = imread(fname,k);
                tifFile.write(uint16(appendImg), strcat(outputFileName,'/',singleFileName,'.tif'), ...
                'writemode', 'append','Compression','none','RowsPerStrip',8);
            end
        end
        tifFile.close();
%         fname2 = strcat(outputFileName,'/',singleFileName,'.tif');
%          info = imfinfo(fname2);
%          info.ImageDescription
         %Delete the files
         tif_files(indx) = [];
        
    end
end

