function [allCells_stack]=cellSumAcrosSpecificTrial(source,cellIndex, allCells, finalOut,plot,mouseExp,sigType)
  
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
    if ~allCells
        
        for(j=1:length(stim_array))
            cell_signals = finalOut{2,2}{1,planeToview}{j,1}{4,2}(:,:,:);
            [frames cells data] = size(cell_signals);

            

                data= [cell_signals{[1:frames],cellIndex,3}];
                grad_data= del2(data);
                data_avg= mean([cell_signals{[1:frames],cellIndex,3}]);
                data_sv = std([cell_signals{[1:frames],cellIndex,3}]);
                
                if strcmp(sigType,'df/f');
                        z_data_perCell= grad_data./data_avg;
%                         df_f_avg = mean(df_f);
%                         df_f_sv = std(df_f);
%                         z_data_perCell= (df_f-df_f_avg).^2./df_f_sv;
                     elseif strcmp(sigType,'df')
                         z_data_perCell= grad_data;
                     elseif strcmp(sigType,'df-dfav')
                        z_data_perCell = grad_data - mean(grad_data);
                    elseif strcmp(sigType,'normal')
                        min_d = min(data);
                        max_d = max(data);
                        z_data_perCell= (data-min_d)./(max_d-min_d);
                    elseif strcmp(sigType,'normal-noravg')
                        min_d = min(data);
                        max_d = max(data);
                        z_data_perCell= (data-min_d)./(max_d-min_d) ;
                        z_data_perCell = z_data_perCell - mean(z_data_perCell);
                    elseif strcmp(sigType,'standard');
                        z_data_perCell = (data-data_avg)./data_sv;
                    elseif strcmp(sigType,'standard-stdavg');
                        z_data_perCell = (data-data_avg)./data_sv;
                        z_data_perCell = z_data_perCell- mean(z_data_perCell);
                    elseif strcmp(sigType,'raw-rawavg');
                        z_data_perCell = (data-data_avg);
                    elseif strcmp(sigType,'raw');
                        z_data_perCell = (data);
                     end
                    
                
                if(mod(stim_array(1,j),unique_stims(1,large_indx))==0)
                    cellSignalAr_largest = [z_data_perCell;cellSignalAr_largest];
                else
                    cellSignalAr_smallest = [z_data_perCell;cellSignalAr_smallest];
                end
        end
    elseif allCells 
        %Consider all cells
        allCells_signal = [];
              cell_signals = finalOut{2,2}{1,planeToview}{1,1}{4,2}(:,:,:);
              [frames_stack cells_stack ~] = size(cell_signals);

        allCells_stack = zeros( size(stim_array,2),cells_stack, frames_stack);
        
        for(j=1:length(stim_array))
            cell_signals = finalOut{2,2}{1,planeToview}{j,1}{4,2}(:,:,:);
            [frames cells data] = size(cell_signals);
               allCells_signal = [];
               for(i=1:cells)
                    data= [cell_signals{[1:frames],i,3}];
                    grad_data= del2(data);

                    data_avg= mean([cell_signals{[1:frames],i,3}]);
                    data_sv = std([cell_signals{[1:frames],i,3}]);
                     if strcmp(sigType,'df/f');
                        z_data_perCell= grad_data./data_avg;
%                         df_f_avg = mean(df_f);
%                         df_f_sv = std(df_f);
%                         z_data_perCell= (df_f-df_f_avg).^2./df_f_sv;
                     elseif strcmp(sigType,'df')
                         z_data_perCell= grad_data;
                     elseif strcmp(sigType,'df-dfav')
                        z_data_perCell = grad_data - mean(grad_data);
                    elseif strcmp(sigType,'normal')
                        min_d = min(data);
                        max_d = max(data);
                        z_data_perCell= (data-min_d)./(max_d-min_d);
                    elseif strcmp(sigType,'normal-noravg')
                        min_d = min(data);
                        max_d = max(data);
                        z_data_perCell= (data-min_d)./(max_d-min_d) ;
                        z_data_perCell = z_data_perCell - mean(z_data_perCell);
                    elseif strcmp(sigType,'standard');
                        z_data_perCell = (data-data_avg)./data_sv;
                    elseif strcmp(sigType,'standard-stdavg');
                        z_data_perCell = (data-data_avg)./data_sv;
                        z_data_perCell = z_data_perCell- mean(z_data_perCell);
                    elseif strcmp(sigType,'raw-rawavg');
                        z_data_perCell = (data-data_avg);
                    elseif strcmp(sigType,'raw');
                        z_data_perCell = (data);
                     end
                    
                    allCells_signal = [z_data_perCell;allCells_signal];
               end
               %Sort Signal
               if  strcmp(sigType,'df') || strcmp(sigType,'df/f') || strcmp(sigType,'df-dfav')
                    [~,sortInd] = sort(sum(abs(allCells_signal(:,105:140)),2),'descend');
               else
                   [~,sortInd] = sort(sum(allCells_signal(:,105:140),2),'descend');
               end
               
               allCells_signal = allCells_signal(sortInd,:);
                if(mod(stim_array(1,j),unique_stims(1,large_indx))==0)
                    allCells_stack(j,:,:) = allCells_signal;
                else
                    allCells_stack(j,:,:) = allCells_signal;
                end
        end

        %allCells_stack - trialNum x cells x frames
        %reundant work but whatever...
        large_stim=[];
        small_stim=[];
       for(i = 1:length(stim_array))
           if(unique_stims(1,large_indx)== stim_array(i))
            large_stim = [allCells_stack(i,:,:) ; large_stim];
           else
             small_stim = [allCells_stack(i,:,:) ; small_stim];
           end
       end
       if strcmp(source,'small')
           objToLoad = small_stim;
       else 
           objToLoad = large_stim;
       end
       play_stim = permute(objToLoad,[2 3 1]);
       degreeSize = size(play_stim,3);
         if(mod(degreeSize,2)==0)
               plot_w_max = degreeSize/2;
               plot_h_max = degreeSize/2;
           else
               plot_w_max = ceil(degreeSize/2);
               plot_h_max = ceil(degreeSize/2)-1;
         end
       colormap('jet')
       cnt = 0;   
     
       for(i=1:plot_h_max)
        for(j=1:plot_w_max)
          cnt = cnt+1;
          if(cnt>(plot_w_max+plot_h_max)) break; end
          subplot(2,plot_w_max ,cnt);
          imagesc(play_stim(:,:,cnt))
          if(strcmp(sigType,'normal-noravg') || strcmp(sigType,'normal' )) caxis([0 1]);
          elseif strcmp(sigType,'raw-rawavg') ||  strcmp(sigType,'df') ||  strcmp(sigType,'df/f') || strcmp(sigType,'df-dfav')  caxis([0 max(max(reshape(allCells_stack(1,:,:),[size(allCells_stack,2) size(allCells_stack,3)])))])
          
          end
          colorbar
          tl = strcat('Trial  ',num2str(cnt))
          title(tl);
%           pause(1);
        end
         if(cnt>(plot_w_max+plot_h_max)) break; end
       end
       
    end
    if plot && ~allCells
     
        freezeColors
       subplot(3,1,1)
       colormap('jet')
       large = sum(cellSignalAr_largest(:,105:140),2);
       [dummy2, index2] = sort(large,'descend');
       
       sorted_large= cellSignalAr_largest(index2,:);
       sum_col_large = sum(sorted_large);
       imagesc(sum_col_large);
       
       t1= strcat(num2str(unique_stims(large_indx)), ' valve summed');
       title(t1);
       freezeColors
       subplot(3,1,2)
       colormap('jet')
       small = sum(cellSignalAr_smallest(:,105:140),2);
       [dummy2, index2] = sort(small,'descend');
       
       sorted_small= cellSignalAr_smallest(index2,:);
       sum_col_small = sum(sorted_small);
       imagesc(sum_col_small);
       t2= strcat(num2str(unique_stims(small_indx)), ' valve summed');
       title(t2)
       
       grp = [cellSignalAr_smallest;cellSignalAr_largest];
       grp_avg =  mean(grp);
       grp_std = std(grp);
       
       large_norm = sum_col_large-grp_avg./grp_std;
       small_norm = sum_col_small-grp_avg./grp_std;
       
       
       diff = large_norm-small_norm;
       freezeColors;
       subplot(3,1,3);
       colormap('cool')
       imagesc(diff);
       t2= 'Standardized valve difference';
       title(t2)
       colorbar;
    end
    
end