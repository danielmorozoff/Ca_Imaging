function modifyName(prefix,direct)
    files= dir(direct);
    for i=1:length(files)
        filename = files(i,1).name; 
        if length(filename)>2
            movefile([direct '/' filename], [direct '/' prefix filename]);
        end
    end
end