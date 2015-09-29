function [differenceInCells,ordering,index_normal,index_scores]= cellDifferenceMatrix(finalOut,mouseExp,sortBy) 
    %Extract stimulus
    table=[];
    if mouseExp==11 table = readtable('s11_experiment_10W_10S_10W_10pause_20150326103013_summary.txt','Delimiter','tab');
    elseif mouseExp==10 table = readtable('s10_experiment_10W_10S_10W_10pause_20150326205548_summary.txt','Delimiter','tab');
    end
    stim_array= generateStimArray(table); 
    
    differenceInCells = [];
    planeToview = 4;
    colormap('cool')
%     for(j=1:length(stim_array))
    j = 1; % doesnt matter just to calculate cells...
        cell_signals = finalOut{2,2}{1,planeToview}{j,1}{4,2}(:,:,:);
        [frames cells data] = size(cell_signals);
        for(i= 1:cells)
            
            summedDiff_cubic  = singleCellAcrossTrials(i,finalOut,false,mouseExp);
            differenceInCells = [summedDiff_cubic; differenceInCells];
        end
%     end

%Whole set normalized
   normalizedDiff = [];
   index_normal = [];
%     normalizedDiff = double(differenceInCells)./double(max(abs(differenceInCells(:))));
%     [dummy index_normal] = sort(sum(normalizedDiff(:,105:140),2),'descend');
    
    
%     [dummy index_normal] = sort(sum(differenceInCells(:,105:140),2),'descend');
    

ordering = [];
if strcmp(sortBy,'2area')
    water_weight_coef = 1;
    signal_weight_coef = 1;
    
    %Absolute value taken in order to account for 2 way signal propagation.
    [lowInFirstWater lowInFirstWater_ind] = sort(sum(water_weight_coef*abs(differenceInCells(:,1:70)),2),'ascend');%Lowest in water window;
    [highInStim highInStim_ind] = sort(sum(signal_weight_coef*abs(differenceInCells(:,105:140)),2),'descend'); %Highest signal window
    
    
    ordering = [lowInFirstWater_ind,highInStim_ind];
    
    %Add index across rows for each unit
    cells_indxAr = [1:cells];
    index_scores = [];
    
    
    for(i=1:cells)
        cell = cells_indxAr(1,i);
        [~, index_water] = ismember(cell, ordering(:,1));
        [~, index_stim] = ismember(cell, ordering(:,2));
        
        index_scores(i,1) = cell;
        index_scores(i,2) = index_water+index_stim;
    end
    
    [dummy,sorted_ind] = sort(index_scores(:,2),'ascend');
    index_scores = index_scores(sorted_ind,:);
    imagesc(differenceInCells(sorted_ind,:));
    title('Cell difference matrix');
    set(gca,'ytick',1:cells);
    set(gca,'YTickLabel',index_scores(:,1));
    colorbar;
elseif strcmp(sortBy,'stimarea')
    [highInStim highInStim_ind] = sort(sum(differenceInCells(:,105:140),2),'descend'); %Highest signal window
    imagesc(differenceInCells(highInStim_ind,:));
    title('Cell difference matrix');
end
end