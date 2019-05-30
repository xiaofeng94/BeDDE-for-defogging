function ChromSimMap = chromine_similarity_map(fogImg, clearImg, type)

if strcmp(type, 'YIQ')
%     Y1 = 0.299 * double(clearImg(:,:,1)) + 0.587 * double(clearImg(:,:,2)) + 0.114 * double(clearImg(:,:,3));
%     Y2 = 0.299 * double(fogImg(:,:,1)) + 0.587 * double(fogImg(:,:,2)) + 0.114 * double(fogImg(:,:,3));
    I1 = 0.596 * double(clearImg(:,:,1)) - 0.274 * double(clearImg(:,:,2)) - 0.322 * double(clearImg(:,:,3));
    I2 = 0.596 * double(fogImg(:,:,1)) - 0.274 * double(fogImg(:,:,2)) - 0.322 * double(fogImg(:,:,3));
    Q1 = 0.211 * double(clearImg(:,:,1)) - 0.523 * double(clearImg(:,:,2)) + 0.312 * double(clearImg(:,:,3));
    Q2 = 0.211 * double(fogImg(:,:,1)) - 0.523 * double(fogImg(:,:,2)) + 0.312 * double(fogImg(:,:,3));

    ChromA1 = I1; ChromA2 = I2;
    ChromB1 = Q1; ChromB2 = Q2;

    consA = 200;
    consB = 200;
else % LMN is default
%     L1 = 0.06 * double(clearImg(:,:,1)) + 0.63 * double(clearImg(:,:,2)) + 0.27 * double(clearImg(:,:,3));
%     L2 = 0.06 * double(fogImg(:,:,1)) + 0.63 * double(fogImg(:,:,2)) + 0.27 * double(fogImg(:,:,3));
    M1 = 0.30 * double(clearImg(:,:,1)) + 0.04 * double(clearImg(:,:,2)) - 0.35 * double(clearImg(:,:,3));
    M2 = 0.30 * double(fogImg(:,:,1)) + 0.04 * double(fogImg(:,:,2)) - 0.35 * double(fogImg(:,:,3));
    N1 = 0.34 * double(clearImg(:,:,1)) - 0.60 * double(clearImg(:,:,2)) + 0.17 * double(clearImg(:,:,3));
    N2 = 0.34 * double(fogImg(:,:,1)) - 0.60 * double(fogImg(:,:,2)) + 0.17 * double(fogImg(:,:,3));

    ChromA1 = M1; ChromA2 = M2;
    ChromB1 = N1; ChromB2 = N2;
    
    consA = 130;
    consB = 130;
end

% Downsample the image
F = 2;
[rows, cols] = size(clearImg(:,:,1));
aveKernel = fspecial('average',F);

aveChromA1 = conv2(ChromA1, aveKernel,'same');
aveChromA2 = conv2(ChromA2, aveKernel,'same');
ChromA1_dn = aveChromA1(1:F:rows,1:F:cols);
ChromA2_dn = aveChromA2(1:F:rows,1:F:cols);

aveChromB1 = conv2(ChromB1, aveKernel,'same');
aveChromB2 = conv2(ChromB2, aveKernel,'same');
ChromB1_dn = aveChromB1(1:F:rows,1:F:cols);
ChromB2_dn = aveChromB2(1:F:rows,1:F:cols);

% calculate feature map
ChromAMap = (2 * ChromA1_dn .* ChromA2_dn + consA) ./ (ChromA1_dn.^2 + ChromA2_dn.^2 + consA);
ChromBMap = (2 * ChromB1_dn .* ChromB2_dn + consB) ./ (ChromB1_dn.^2 + ChromB2_dn.^2 + consB);

ChromSimMap = ChromAMap.*ChromBMap;

end


