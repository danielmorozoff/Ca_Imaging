function stim_array = generateStimArray(table)
if(nargin==0) table = readtable('s11_experiment_10W_10S_10W_10pause_20150326103013_summary.txt','Delimiter','tab');end
%Extract stimulus
stim_array= [];
    for(i=8:1:height(table))
        valInCell = char(table(i,1).tastant1Is);
        values = strsplit(char(valInCell),',');
        num = str2num(char(values(1,2)));
        if(num==0) break; end
        stim_array = [num,stim_array]; 
    end
end