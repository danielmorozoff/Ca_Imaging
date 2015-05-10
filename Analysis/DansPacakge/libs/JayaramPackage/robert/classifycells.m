function [output] = classifycells(input,skipsucrose)
    % expects
    % 1: sweet
    % 2: bitter
    % 3: salt
    % 4: umami
    % 5: sour

    %% provides a taste score for each input cell response profile
    if ~exist('skipsucrose')
        skipsucrose = 0;
    end
    for i=1:size(input,1)
        row = input(i,:);
        score(i,1) = issweet(row);
        score(i,2) = isbitter(row);
        score(i,3) = issalt(row);
        score(i,4) = issour(row);
        score(i,5) = isumami(row);    
    end
    output = score;


    %% subfunction logic for determining taste responsiveness
    function [tf] = issweet(row)
        if skipsucrose
            % responds to acek
            if (row(1) == 1) 
                tf = 1;
            else
                tf = 0;
            end            
        else
            % responds to acek and sucrose
            if (row(1) == 1) && (row(7) ==1 ) 
                tf = 1;
            else
                tf = 0;
            end
        end
    end            
            
    function [tf] = isbitter(row)
        % responds to bitter
        if (row(2) == 1)
            tf = 1;
        else
            tf = 0;
        end
    end
    
    function [tf] = issalt(row)
        if (row(3) == 1)
            tf = 1;
        else
            tf = 0;
        end
    end

    function [tf] = isumami(row)
        if (row(4) == 1)
            tf = 1;
        else
            tf = 0;
        end
    end

    function [tf] = issour(row)
        if (row(5) == 1)
            tf = 1;
            if length(row)>=12
                if (row(2)==1) && (row(12) == 0)
                    tf = 0;
                end
            end
        else
            tf = 0;
        end
    end
end

