function imageStack=SpatialMedianFilter(imageStack,window)
%Set up variables
width = size(imageStack,1);
height = size(imageStack,2);
stacksize = size(imageStack,3);

for i=1:stacksize
    imageStack(:,:,i) = medfilt2(imageStack(:,:,i),[window window]);
end

end
