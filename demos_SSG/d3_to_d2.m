function [img rows cols bands] = d3_to_d2(dat, rot)
if rot == 0, [rows cols bands] = size(dat);
else [cols rows bands] = size(dat);
end
img = zeros(bands, rows*cols);
for i = 1 : bands
    if rot == 0, x = dat(:,:,i);
    else x = dat(:,:,i)';
    end
    x = reshape(x,1,rows*cols);
    img(i,:) = x;
end
end