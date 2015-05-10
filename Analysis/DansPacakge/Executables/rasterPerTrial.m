function [z_data_full,index2]= rasterPerTrial(trialNum,plane,finalOut,sigType)

 cell_signals = finalOut{2,2}{1,plane}{trialNum,1}{4,2}(:,:,:);
 [frames cells data] = size(cell_signals);

 colormap('jet');
    z_data_full= [];
    for(i= 1:cells)
        data= [cell_signals{[1:frames],i,3}];

        data_avg= mean([cell_signals{[1:frames],i,3}]);
        data_sv = std([cell_signals{[1:frames],i,3}]);
        data_to_load =[];
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
         z_data_full = [data_to_load;z_data_full];
    %     plot(x,z_data,'Color',jet_map(i,:));
    end
    
    %  order by total signal
      subplot(1,1,1);
      z_data_full_max2 = sum(z_data_full(:,70:140),2);
      [dummy2, index2] = sort(z_data_full_max2,'descend');
      imagesc(z_data_full(index2, :));
      title('Raster by summed signal 5 seconds post stim onset')
      set(gca,'ytick',1:cells);
      set(gca,'YTickLabel',index2);
      if(strcmp(sigType,'normal-noravg') || strcmp(sigType,'normal' )) caxis([0 1]);
      elseif strcmp(sigType,'raw-rawavg') ||  strcmp(sigType,'df') ||  strcmp(sigType,'df/f') || strcmp(sigType,'df-dfav')  caxis([0 max(max(reshape(allCells_stack(1,:,:),[size(allCells_stack,2) size(allCells_stack,3)])))])
      end
      colorbar;
      
%       %Order by highest signal point
%      colormap('jet');
%      subplot(2,1,2);
%      z_data_full_max = max(z_data_full, [], 2);
%      [dummy, index] = sort(z_data_full_max,'descend');
%      imagesc(z_data_full(index, :));
%      title('Raster by top signal')
%      set(gca,'ytick',1:cells);
%      set(gca,'YTickLabel',index);
%      colorbar;
end