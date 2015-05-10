function [cellSignalAr_acek,cellSignalAr_quinine]=stateSpaceGen (finalOut,sigType)
 table = readtable('s10_experiment_10W_10S_10W_10pause_20150326205548_summary.txt','Delimiter','tab');
 stim_array= generateStimArray(table);
 
  cellSignalAr_acek = [];
    cellSignalAr_quinine = [];
    planeToview = 4;
    length(stim_array)
    for(j=1:length(stim_array))
        cell_signals = finalOut{2,2}{1,planeToview}{j,1}{4,2}(:,:,:);
        [frames cells data] = size(cell_signals);
        
    col_qui=[];
    col_acek=[];
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
         
         data_to_load = data_to_load(1,70:140);
         
        if(mod(stim_array(1,j),6)==0)
            col_qui = [data_to_load;col_qui];
        else
            col_acek = [data_to_load;col_acek];
        end
        
    end
      cellSignalAr_quinine = [cellSignalAr_quinine col_qui];
      cellSignalAr_acek = [cellSignalAr_acek col_acek];  
    end
    

end