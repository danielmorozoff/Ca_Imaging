%% this function takes tuning input as data and calculates cell type

function [output] = analysisgcamp6(input)
% expects 9 columns as input
output = classifycells(input,1);

% total number of cells
totalcells = sum(sum(input,2)>0);

% number of each tastant irrespective of tuning breadth
for i=1:5
    totalnumbers(i) = length(find(output(:,i)==1));
end
% umami cells
umamicells = find(output(:,4)==1);
umamiscores = output(find(output(:,4)==1),:);


% sweet cells
sweetcells = find(output(:,1)==1);
sweetscores = output(find(output(:,1)==1),:);

% salt cells
saltcells = find(output(:,3)==1);
saltscores = output(find(output(:,3)==1),:);

% sour cells
sourcells = find(output(:,5)==1);
sourscores = output(find(output(:,5)==1),:);

% bitter cells
bittercells = find(output(:,2)==1);
bitterscores = output(find(output(:,2)==1),:);

bitteronlycells = bittercells(ismember(bitterscores,[0 1 0 0 0],'rows'));
bittersouronlycells = bittercells(ismember(bitterscores,[0 1 0 0 1],'rows'));

bitteronlynum = length(bitteronlycells);
bittersouronlynum = length(bittersouronlycells);
keyboard
% aitcsensitivecells = find(input(:,8)==1);
end