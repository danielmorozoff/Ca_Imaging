% tophat calculation
% goal is define a 2D filter that when multiplied with an random image,
% generates 0.  and when multiplied with a center heavy image generates a
% large value.

in = Simg(:,:,1);
datadim = size(in);
xrad = (datadim(1)-1)/2;
yrad = (datadim(2)-1)/2;
for i=1:size(Simg,3)
    in = Simg(:,:,i);
    xweights = repmat([1:datadim(1)]',[1 datadim(2)]);
    yweights = repmat([1:datadim(2)],[datadim(1) 1]);
    xcenter = mean(mean(xweights));
    ycenter = mean(mean(in.*yweights));
    out = sqrt(xcenter^2+ycenter^2);
    disp(out)
end