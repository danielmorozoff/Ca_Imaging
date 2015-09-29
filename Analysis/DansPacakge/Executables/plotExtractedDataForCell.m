function []=plotExtractedDataForCell(cellIndex,finalOut,plot,mouseExp)
  
    %Extract stimulus
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
   
   
    cellSignalAr_smallest = [];
    cellSignalAr_largest = [];
    for(j=1:length(stim_array))
        cell_signals = finalOut{2,2}{1,planeToview}{j,1}{4,2}(:,:,:);
        [frames cells data] = size(cell_signals);

        z_data_full= [];
        
            data= [cell_signals{[1:frames],cellIndex,3}];
            
            data_avg= mean([cell_signals{[1:frames],cellIndex,3}]);
            data_sv = std([cell_signals{[1:frames],cellIndex,3}]);
            z_data = (data-data_avg)./data_sv;
            
            if(mod(stim_array(1,j),unique_stims(1,large_indx))==0)
                cellSignalAr_largest = [z_data;cellSignalAr_largest];
            else
                cellSignalAr_smallest = [z_data;cellSignalAr_smallest];
            end
    end
    if plot
       subplot(2,1,1)
      
       large = sum(cellSignalAr_largest(:,105:140),2);
       [dummy2, index2] = sort(large,'descend');
       imagesc(cellSignalAr_largest(index2,:));
       t1= strcat(num2str(unique_stims(large_indx)), ' valve');
       title(t1);
       
       subplot(2,1,2)
       small = sum(cellSignalAr_smallest(:,105:140),2);
       [dummy2, index2] = sort(small,'descend');
       imagesc(cellSignalAr_smallest(index2,:));
       t2= strcat(num2str(unique_stims(small_indx)), ' valve');
       title(t2)
    end
    
end