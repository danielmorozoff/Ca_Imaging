function [summedDiff_cubic,cellSignalAr_quinine,cellSignalAr_sucrose,comb]=singleCellAcrossTrials(cellIndex,finalOut,plot,mouseExp,sigType)
  
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
            data_to_load= [] ;
                    if strcmp(sigType,'df/f');
                        data_to_load= grad_data./data_avg;
%                         df_f_avg = mean(df_f);
%                         df_f_sv = std(df_f);
%                         z_data_perCell= (df_f-df_f_avg).^2./df_f_sv;
                    elseif strcmp(sigType,'df-dfav')
                        data_to_load = grad_data - mean(grad_data);
                    elseif strcmp(sigType,'normal')
                        min_d = min(data);
                        max_d = max(data);
                        data_to_load= (data-min_d)./(max_d-min_d);
                    elseif strcmp(sigType,'normal-noravg')
                        min_d = min(data);
                        max_d = max(data);
                        data_to_load= (data-min_d)./(max_d-min_d) ;
                        data_to_load = data_to_load - mean(data_to_load);
                    elseif strcmp(sigType,'standard');
                        data_to_load = (data-data_avg)./data_sv;
                    elseif strcmp(sigType,'standard-stdavg');
                        data_to_load = (data-data_avg)./data_sv;
                        data_to_load = data_to_load- mean(data_to_load);
                    elseif strcmp(sigType,'raw-rawavg');
                        data_to_load = (data-data_avg);
                    elseif strcmp(sigType,'raw');
                        data_to_load = (data);
                     end
                    
            
            if(mod(stim_array(1,j),unique_stims(1,large_indx))==0)
                cellSignalAr_largest = [data_to_load;cellSignalAr_largest];
            else
                cellSignalAr_smallest = [data_to_load;cellSignalAr_smallest];
            end
    end
    
   
%     figure;
%     figTitle = strcat('Cell: ',num2str(cellIndex));
%     title(figTitle);
    
   %reduce largest data group to match smallest
    cellSignalAr_largest_reduced = cellSignalAr_largest(1:smallest,:);
 
   %normalize pairwise trials against one another since scales are different
%    diff = [];
%    for(i=1:7)
%        comb = [cellSignalAr_quinine(i), cellSignalAr_sucrose(i,:)];
%        normal_qui = double(cellSignalAr_quinine)./double(max(comb));
%        normal_suc = double(cellSignalAr_sucrose)./double(max(comb));
%        diff(i,:) = .5*(normal_qui(i,:)-normal_suc(i,:)).^3;
%    end

   %Normalized against entire combined object
   diff = [];
   normal_largest = [];
   normal_smallest = [];
   
   comb = [cellSignalAr_largest_reduced; cellSignalAr_smallest];
   col_maxes = max(abs(comb)); %abs because you want to consider the negative deviations
   for(i=1:size(cellSignalAr_largest_reduced,1))
       normal_largest = [cellSignalAr_largest_reduced(i,:)./col_maxes;normal_largest];
       normal_smallest = [cellSignalAr_smallest(i,:)./col_maxes;normal_smallest];
   end

   %DIfferent normalization  min-max 
%    normal_largest =  cellSignalAr_largest_reduced - min(cellSignalAr_largest_reduced)/(max(cellSignalAr_largest_reduced) - min(cellSignalAr_largest_reduced));
%    normal_smallest = cellSignalAr_smallest - min(cellSignalAr_smallest)/(max(cellSignalAr_smallest) - min(cellSignalAr_smallest));
%    
   square_diff = .5*(normal_largest-normal_smallest).^2;
   cubic_diff = .5*(normal_largest-normal_smallest).^3;
   summedDiff_square = sum(square_diff);
   summedDiff_cubic = sum(cubic_diff);
   if plot
    subplot(4,1,1);
    
    colormap('jet')
    larg = sum(normal_largest(:,70:140),2);
    [dummy, index] = sort(larg,'descend');
    imagesc(normal_largest(index,:));
    t1= strcat(num2str(unique_stims(large_indx)), 'valve');
    title(t1)
    cb1= colorbar;
    if(strcmp(sigType,'normal-noravg') || strcmp(sigType,'normal' )) caxis([0 1]);
    elseif strcmp(sigType,'raw-rawavg') ||  strcmp(sigType,'df') ||  strcmp(sigType,'df/f') || strcmp(sigType,'df-dfav')  caxis([0 max(max(reshape(allCells_stack(1,:,:),[size(allCells_stack,2) size(allCells_stack,3)])))])
    end
    
    freezeColors;
%     
    subplot(4,1,2)
    smal = sum(normal_smallest(:,105:140),2);
    [dummy2, index2] = sort(smal,'descend');
    imagesc(normal_smallest(index2,:));
    t2= strcat(num2str(unique_stims(small_indx)), 'valve');
    title(t2)
    cb2= colorbar;
     if(strcmp(sigType,'normal-noravg') || strcmp(sigType,'normal' )) caxis([0 1]);
    elseif strcmp(sigType,'raw-rawavg') ||  strcmp(sigType,'df') ||  strcmp(sigType,'df/f') || strcmp(sigType,'df-dfav')  caxis([0 max(max(reshape(allCells_stack(1,:,:),[size(allCells_stack,2) size(allCells_stack,3)])))]) 
     end
    
    freezeColors;
    
    subplot(4,1,3)
    colormap('bone');
    imagesc(summedDiff_square);
    title('Normalized Square Difference')
    cb3= colorbar;
    
    freezeColors;
    subplot(4,1,4)
    colormap('cool');
    imagesc(summedDiff_cubic);
    title('Normalized Cubic Difference')
    cb4 = colorbar;
   end
   
   
%     x0=10;y0=10;width=550;height=20;
%     set(gcf,'units','points','position',[x0,y0,width,height])
%     d = log10(diff);
%     mn = min(d(:));
%     rng = max(d(:))-mn;
%     d = 1+63*(d-mn)/rng; % Self scale data

%     hC = colorbar;
%     L = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000];
%     % Choose appropriate
%     % or somehow auto generate colorbar labels
%     l = 1+63*(log10(L)-mn)/rng; % Tick mark positions
%     set(hC,'Ytick',l,'YTicklabel',L);
    
    
    
end