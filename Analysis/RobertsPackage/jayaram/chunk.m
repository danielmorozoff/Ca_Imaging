% this was a script used to truncate the cellsig to the first 70 frames as
% there was a problem with the data acquisition.
folder = 'C:\Users\rpjb\Dropbox\Share (from jayaram)\check\to_truncate\WS25\04202013\fov_01004';
load(fullfile(folder,'celltraces.mat'))

oldcellsig = cellsig;
for i=1:size(cellsig,1)
    for j=1:size(cellsig,2)
        cellsig{i,j} = oldcellsig{i,j}(1:70);
    end
end

savefile = fullfile(folder,'newcelltraces.mat');
save(savefile,'cellsig','cellspatialfilters','cellx','celly','cellcentroidx','cellcentroidy','cellarea','xrad','yrad','cellxrad','cellyrad')