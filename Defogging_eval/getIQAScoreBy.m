function score = getIQAScoreBy(clearImg, defogImg, mask, eval_method)

    if strcmp(eval_method, 'GMSD')
        score = GMSD_mask(double(rgb2gray(clearImg)), double(rgb2gray(defogImg)), mask);
    elseif strcmp(eval_method, 'MDSI')
        score = MDSI_mask(clearImg, defogImg, mask); % ignore this metrics
    elseif strcmp(eval_method, 'SSIM')
        score = ssim_mask(rgb2gray(defogImg), rgb2gray(clearImg), mask); % keep uint8 input
    elseif strcmp(eval_method, 'PSNR')
        curMask_rgb = repmat(mask,[1,1,3]);
        score = psnr(defogImg(curMask_rgb), clearImg(curMask_rgb)); % keep uint8 input
    elseif strcmp(eval_method, 'VSI')
        score = VSI_mask(defogImg, clearImg, mask); % keep uint8 input
    elseif strcmp(eval_method, 'FSIM')
        score = FeatureSIM_mask(defogImg, clearImg, mask, 'FSIM'); % keep uint8 input
    elseif strcmp(eval_method, 'FSIMc')
        score = FeatureSIM_mask(defogImg, clearImg, mask, 'FSIMc'); % keep uint8 input
    elseif strcmp(eval_method, 'MSE')
        score = MSE_mask(defogImg, clearImg, mask); % keep uint8 input
    elseif any(strcmp(eval_method, {'Visual', 'VI'}))
        score = visual_index(defogImg, clearImg, mask);
    elseif any(strcmp(eval_method, {'Realness', 'RI'}))
        score = realness_index(defogImg, clearImg, mask);
    elseif strcmp(eval_method, 'FADE')
        score = FADE(defogImg);
    else
        error(sprintf('No implementation for %s', eval_method));
    end

end