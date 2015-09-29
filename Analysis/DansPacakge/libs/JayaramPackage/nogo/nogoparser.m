% nogo parser

function out = nogoparser(fn,fdir)

%% select filename
if isempty(fn)
    fn = 'WS27_QUI_GO_Suc_NOGO_12_20131120190811_summary.txt';
    fdir = '/Volumes/usb/JC_GO_NOGO/11202013';
    fn = fullfile(fdir,fn);
end
%% gather data
fid = fopen(fn);
% grab valve labels
tastelabels = textscan(fid, 'tastant%d is %s','BufSize',4096*16);
% skip definitions
numlinestoskip = 1;
% grab scores
data = textscan(fid, '%d, %d, %d, %d','BufSize',4096*16,'HeaderLines',numlinestoskip);
fclose(fid);


out.num = data{1};
out.valve = data{2};
out.fluid = tastelabels{2};
out.type = data{3};
out.result = data{4};


    

