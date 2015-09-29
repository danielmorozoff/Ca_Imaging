function  activityHeatMapForAllTrials(finalOut, planeToView,plotZero, mouseExp)
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

for i = 1:4
    cnt = cnt+1;
    subplot(2,2 ,cnt);
    plotCellActivityHeatMap('session_maxproj_chan_01',finalOut, i, planeToView,true,plotZero, 10,'normal-noravg') ;
end

end