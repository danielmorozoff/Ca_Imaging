function getUniqueCellIndicesForTrialType(finalOut,mouseExp)


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

cell_signals = finalOut{2,2}{1,planeToview}{trialNum,1}{4,2}(:,:,:);

    for(i=1:length(stim_array))
        type = stim_array(i);
        if type == 6
            
        elseif type == 4 

        end
        
        for(cellIndex =1:cells)
            data= [cell_signals{[1:frames],cellIndex,3}];
            data_to_load = [];

                if strcmp(sigType,'df/f');
                        data_to_load= gradient(data)./data_avg;
    %                         df_f_avg = mean(df_f);
    %                         df_f_sv = std(df_f);
    %                         z_data_perCell= (df_f-df_f_avg).^2./df_f_sv;
                    elseif strcmp(sigType,'df-dfav')
                        grad_data = gradient(data);
                        data_to_load = grad_data - mean(grad_data);
                    elseif strcmp(sigType,'avg')
                        data_to_load = data_avg;
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
        summedData =  sum(data_to_load(1,105:140));
        if(summedData<0) summedData=0;end


       
        end
    end
end