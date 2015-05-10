function [cellSignalAr_quinine,cellSignalAr_acek,summedDiff]=collapsedMultiCellAcrossTrials(finalOut)
 table = readtable('s11_experiment_10W_10S_10W_10pause_20150326103013_summary.txt','Delimiter','tab');
 stim_array= generateStimArray(table);
 
  cellSignalAr_acek = [];
    cellSignalAr_quinine = [];
    planeToview = 4;
    for(j=1:length(stim_array))
        cell_signals = finalOut{2,2}{1,planeToview}{j,1}{4,2}(:,:,:);
        [frames cells data] = size(cell_signals);
        
%         jet_map = colormap(hsv(cells));
        dataAvg_col= [];
        
        for(i= 1:cells)
        data= [cell_signals{[1:frames],i,3}];

        data_avg= mean([cell_signals{[1:frames],i,3}]);
        data_sv = std([cell_signals{[1:frames],i,3}]);
        z_data = (data-data_avg)./data_sv;
        z_avg = mean(z_data);
        dataAvg_col = [z_avg;dataAvg_col];

        end
        if(mod(stim_array(1,j),6)==0)
            cellSignalAr_quinine = [dataAvg_col,cellSignalAr_quinine];
        else
            cellSignalAr_acek = [dataAvg_col,cellSignalAr_acek];
        end
    end
    colormap('jet')
    subplot(1,3,1)
    z_data_full_max = sum(cellSignalAr_quinine,2);
    [dummy2, index] = sort(z_data_full_max,'descend');
    imagesc(cellSignalAr_quinine(index,:));
    set(gca,'ytick',1:cells);
    set(gca,'YTickLabel',index);
    title('Quinine');
    freezeColors;
    subplot(1,3,2)
    imagesc(cellSignalAr_acek(index,:));
    set(gca,'ytick',1:cells);
    set(gca,'YTickLabel',index);
    title('AceK');
    freezeColors;
%Find maximal difference... in 11C7 combination...
%     subplot(1,3,3)
%     comb_vec = mat2str([1:11]);
%     comb_vec = comb_vec(find(~isspace(comb_vec)));
%     comb_vec = comb_vec(2:end-1);
% 
%     maxDiffInd = 0;
%     maxDiff = 0;
%       for ii = idx'
%         temp_mat = cellSignalAr_quinine(:,ii);
%         curDiff = .5*(temp_mat-cellSignalAr_acek).^2;
%         curDiffSum = sum(sum(curDiff));
%         fprintf('CurDiff %d  vs maxDiff %d \n',curDiffSum,maxDiff);
%         if(ii~=1)
%             if(curDiffSum > maxDiff) maxDiff = curDiffSum;maxDiffInd = ii; end
%         end
%       end
    colormap('cool');
    subplot(1,3,3)
    difference_mat = .5*(cellSignalAr_quinine(index,1:7)-cellSignalAr_acek(index,:)).^3;  
    summedDiff = sum(difference_mat,2);
    [sorted_sum,sum_indx] = sort(summedDiff,'descend');
    imagesc(sorted_sum);
    
    set(gca,'ytick',1:cells);
    set(gca,'YTickLabel',sum_indx);
    title('Difference');
    
    
end