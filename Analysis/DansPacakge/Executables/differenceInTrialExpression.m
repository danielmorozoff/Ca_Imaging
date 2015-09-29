function [t1_out,t2_out,inter_out]=differenceInTrialExpression(finalOut,mouseExp, trial1, trial2)
table=[];
if mouseExp==11 table = readtable('s11_experiment_10W_10S_10W_10pause_20150326103013_summary.txt','Delimiter','tab');
elseif mouseExp==10 table = readtable('s10_experiment_10W_10S_10W_10pause_20150326205548_summary.txt','Delimiter','tab');
end
stim_array= generateStimArray(table);
planeToview = 4;
 %Get Unique stimulants
[repeated,unique_stims] =  hist(stim_array,unique(stim_array));
[smallest,small_indx] = min(repeated);
[largest,large_indx] = max(repeated);    


degreeSize = length(stim_array);
plot_w_max = ceil(degreeSize/2);
cnt = 0;

planeToView = 4;

[trial1_color,trial1_circle]=plotCellActivityHeatMap('session_maxproj_chan_01',finalOut, trial1, planeToView,false,true, 10,'normal-noravg'); 
[trial2_color,trial2_circle]=plotCellActivityHeatMap('session_maxproj_chan_01',finalOut, trial2, planeToView,false,true, 10,'normal-noravg'); 

t1 = [trial1_circle,trial1_color];
t2 = [trial2_circle,trial2_color];
%Remove zero cells
t1(find(t1(:,3)==0),:) = [];
t2(find(t2(:,3)==0),:) = [];


[inter_grp,ia,~] =  intersect(t1(:,1:2),t2(:,1:2),'rows');


[uniqueTrial_1,i1] =  setdiff(t1(:,1:2),inter_grp,'rows');
[uniqueTrial_2,i2] = setdiff(t2(:,1:2),inter_grp,'rows');

t1_out = t1(i1,:); % unique t1
t2_out= t2(i2,:); %unique t2 
inter_out= t1(ia,:); %overlapping
% size(inter_grp)
% size(uniqueTrial_1)
% size(uniqueTrial_2)

end