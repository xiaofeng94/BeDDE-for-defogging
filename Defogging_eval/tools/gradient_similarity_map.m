function gradSimMap = gradient_similarity_map(fogImg, clearImg)

if ndims(clearImg) == 3 %images are colorful
    Y1 = 0.299 * double(clearImg(:,:,1)) + 0.587 * double(clearImg(:,:,2)) + 0.114 * double(clearImg(:,:,3));
    Y2 = 0.299 * double(fogImg(:,:,1)) + 0.587 * double(fogImg(:,:,2)) + 0.114 * double(fogImg(:,:,3));
else %images are grayscale
    Y1 = clearImg;
    Y2 = fogImg;
end
Y1 = double(Y1);
Y2 = double(Y2);

% Downsample the image
F = 2;
[rows, cols] = size(clearImg(:,:,1));
aveKernel = fspecial('average',F);

aveY1 = conv2(Y1, aveKernel,'same');
aveY2 = conv2(Y2, aveKernel,'same');
Y1 = aveY1(1:F:rows,1:F:cols);
Y2 = aveY2(1:F:rows,1:F:cols);

% calculate similarity map
dx = [3 0 -3; 10 0 -10;  3  0 -3]/16;
dy = [3 10 3; 0  0   0; -3 -10 -3]/16;
IxY1 = conv2(Y1, dx, 'same');     
IyY1 = conv2(Y1, dy, 'same');    
gradientMap1 = sqrt(IxY1.^2 + IyY1.^2);

IxY2 = conv2(Y2, dx, 'same');     
IyY2 = conv2(Y2, dy, 'same');    
gradientMap2 = sqrt(IxY2.^2 + IyY2.^2);

T = 160; %fixed
gradSimMap = (2*gradientMap1.*gradientMap2 + T) ./(gradientMap1.^2 + gradientMap2.^2 + T);

end