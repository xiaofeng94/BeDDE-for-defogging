function score = psnr_mask(defogImg, clearImg, mask)

chl = size(defogImg,3);

if chl == 3
    mask_rgb = repmat(mask, [1,1,3]);
    score = psnr(defogImg(mask_rgb), clearImg(mask_rgb));
else
    score = psnr(defogImg(mask), clearImg(mask));
end

end