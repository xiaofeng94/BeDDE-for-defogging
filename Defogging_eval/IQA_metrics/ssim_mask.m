function [score, ssimmap] = ssim_mask(A, REF, mask)
% SSIM with mask
% refer to ssim for SSIM details

[~, ssimmap] = ssim(A, REF);

chl = size(A,3);

if chl == 3
    mask_rgb = repmat(mask, [1,1,3]);
    score = mean(ssimmap(mask_rgb));
else
    score = mean(ssimmap(mask));
end

end