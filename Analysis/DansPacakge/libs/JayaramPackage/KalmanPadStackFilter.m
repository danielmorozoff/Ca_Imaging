function [ filtereddata ] = KalmanPadStackFilter( originaldata )
        d = size(originaldata);
        r = zeros(d(1),d(2),d(3)+10);
        r(:,:,1:10) = double(originaldata(:,:,1:10));
        r(:,:,11:end) = double(originaldata);
        % run kalman stack and truncate first ten frames
        filtereddata = KalmanStackFilter(r,.8,.05);
        filtereddata = filtereddata(:,:,11:end);
end

