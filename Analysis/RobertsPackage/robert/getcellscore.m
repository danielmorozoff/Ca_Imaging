function [celltype] = getcellscore(score,stims)

    % organize tastants to the following order
    % 1: sweet
    % 2: bitter
    % 3: salt
    % 4: sour
    % 5: umami
    baseorder = {'sweet','bitter','low salt','sour','umami'};
    order = arrayfun(@(taste) min(find(ismember(stims,taste))),baseorder);    
    combo = score(order);

    % populate a combination table
    combinations = [];
    for k = 1:5
        pos = nchoosek([1:5],k);
        for j = 1:size(pos,1)
            row = [0 0 0 0 0];
            row(pos(j,:)) = 1;
            combinations = [combinations;row];
        end
    end

    celltype = find(ismember(combinations,combo,'rows'));
    
end

