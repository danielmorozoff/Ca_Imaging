% lindseyprep

% get input folder
folder = Filemacpc('/Users/rpjb/Desktop/foliate_2/');

% expect a single aligned file
file = FindFiles(folder,'\S+.tif',2);
oldfile = fullfile(folder,file.name);
imginfo = imfinfo(oldfile);
numframes = length(imginfo);

% expect a 30,15,30 timing in seconds and 2Hz frame rate
framerate = 2; % frames per second
interval = 75;
startframe = [1:framerate*interval:numframes]';
endframe = [framerate*interval:framerate*interval:numframes]';
if length(startframe)>length(endframe)
    endframe(end+1) = numframes-1; %leave out last frame because it is bad
end

% chunk and save individual trials
for i=1:length(startframe)
    [~,basename,extname] = fileparts(file.name);    
    newfile = fullfile(folder,['aligned','_',num2str(i),extname]);
    
    % load data and save to separate chunk file
    for j=startframe(i):endframe(i)
        data = imread(oldfile,j);
        imwrite(data,newfile,'Compression','none','WriteMode','append');        
        stackdata(:,:,j-startframe(i)+1) = double(data);
    end    
       
    % create AVG image
    avgdata = mean(stackdata,3);
    newfile = fullfile(folder,['AVG_', basename,extname]);
    imwrite(data,newfile,'Compression','none','WriteMode','append');        
    
    % create STD image
    stddata = std(stackdata,0,3);
    newfile = fullfile(folder,['STD_', basename,extname]);
    imwrite(data,newfile,'Compression','none','WriteMode','append');             
    
    % clear stackdata
    stackdata = [];
end

