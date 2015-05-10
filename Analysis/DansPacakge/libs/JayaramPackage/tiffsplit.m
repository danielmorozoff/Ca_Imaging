% TIF file splitter

input_file = '/Users/rpjb/Dropbox/project - calcium sample/dave/Stack.tif';

% calculate sets to use
imginfo = imfinfo(input_file);
numframes = length(imginfo);
setsize = 40;
numsets = ceil(numframes/setsize);
for set=1:numsets
    % calculate frames to use
    minframe = 1+setsize*(set-1);
    maxframe = min(numframes,setsize*set);

    % create new file to save
    [basedir,basename,basetype] = fileparts(input_file);
    output_file = fullfile(basedir,[basename,'_',num2str(set),basetype]);
        
    % read/write frames
    for frame = minframe:maxframe
        x{set}(1+frame-minframe) = frame;
        data = imread(input_file,frame);
        imwrite(data,output_file,'Compression','none','WriteMode','append');
    end   
    
end
    
    
