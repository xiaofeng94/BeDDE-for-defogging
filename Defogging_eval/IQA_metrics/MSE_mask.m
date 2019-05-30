function diffVal = MSE_mask(img1,img2,mask)

imgRef = double(img1);
imgDist = double(img2);
tempMask = mask;

[rows, cols, chl] = size(imgRef);

if any(size(imgRef) ~= size(imgDist))
    imgDist = imresize(imgDist, [rows, cols], 'bilinear');
    imgDist = max(imgDist,0);
end

if any(size(tempMask) ~= [rows, cols])
    tempMask = imresize(tempMask, [rows, cols], 'bilinear');
end

if chl == 1
    diffVal = mse(imgRef(tempMask),imgDist(tempMask));
else
    tempMask_rgb = repmat(tempMask, [1,1,3]);
    diffVal = mse(imgRef(tempMask_rgb),imgDist(tempMask_rgb));
%     diffVal = sum((imgRef(tempMask_rgb)-imgDist(tempMask_rgb)).^2)/(rows*cols*chl);
end
    
end