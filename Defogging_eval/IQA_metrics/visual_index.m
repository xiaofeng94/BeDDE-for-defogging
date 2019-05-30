function score = visual_index(fogImg, clearImg, mask)

fogImg_dou = im2double(fogImg);
clearImg_dou = im2double(clearImg);

% compute transmission
omega = 1;
win_size = 15;
dark_chl_clear = get_dark_channel(clearImg_dou, win_size);
air_clear = get_atmosphere(clearImg_dou, dark_chl_clear);
trans_clear = get_transmission_estimate(clearImg_dou, air_clear, omega, win_size);

dark_chl_fog = get_dark_channel(fogImg_dou, win_size);
air_fog = get_atmosphere(fogImg_dou, dark_chl_fog);
trans_fog = get_transmission_estimate(fogImg_dou, air_fog, omega, win_size);

% downsample featrue maps
F = 2;
[rows, cols] = size(clearImg_dou(:,:,1));
aveKernel = fspecial('average',F);

ave_trans_clear = conv2(trans_clear, aveKernel,'same');
ave_trans_fog = conv2(trans_fog, aveKernel,'same');
trans_clear_dn = ave_trans_clear(1:F:rows,1:F:cols);
trans_fog_dn = ave_trans_fog(1:F:rows,1:F:cols);

% apply mask
mask_dn  = imresize(mask, size(trans_clear_dn));
trans_clear_dn_mask = trans_clear_dn(mask_dn);

% calculate similarity score
T1 = mean(mean(trans_clear_dn_mask));  %fixed
transMap = (2 * trans_clear_dn .* trans_fog_dn + T1) ./ (trans_clear_dn.^2 + trans_fog_dn.^2 + T1);


gradSimMap = gradient_similarity_map(fogImg, clearImg);

alpha_1 = 0.4; alpha_2 = 1;

maxTransMap = max(1-trans_clear_dn, 1-trans_fog_dn);
featureVector = real((real(transMap(mask_dn)).^alpha_1) .* (real(gradSimMap(mask_dn)).^alpha_2) .* maxTransMap(mask_dn));

score = sum(featureVector)./sum(maxTransMap(mask_dn));

end

