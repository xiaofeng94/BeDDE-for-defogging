function sim = VSI_mask(image1, image2, mask)
% ========================================================================
% Visual Saliency based Index
% Copyright(c) 2013 Lin ZHANG
% All Rights Reserved.
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereQ
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
% For more information, please refer to
% Lin Zhang et al., "VSI: A Visual Saliency Induced Index for Perceptual
% Image Quality Assessment", submitted to TIP
%----------------------------------------------------------------------
%
%Input : (1) image1: the first image being compared, which is a RGB image
%        (2) image2: the second image being compared, which is a RGB image
%
%Output: sim: the similarity score between two images, a real number
%        
%-----------------------------------------------------------------------

constForVS = 1.27;%fixed
constForGM = 386;%fixed
constForChrom = 130;%fixed
alpha = 0.40;%fixed
lambda = 0.020;%fixed
sigmaF = 1.34;%fixed donot change
omega0 = 0.0210;%fixed
sigmaD = 145;%fixed
sigmaC = 0.001;%fixed

%compute the visual saliency map using SDSP
saliencyMap1 = SDSP(image1,sigmaF,omega0,sigmaD,sigmaC);
saliencyMap2 = SDSP(image2,sigmaF,omega0,sigmaD,sigmaC);

[rows, cols, junk] = size(image1);
%transform into an opponent color space
L1 = 0.06 * double(image1(:,:,1)) + 0.63 * double(image1(:,:,2)) + 0.27 * double(image1(:,:,3));
L2 = 0.06 * double(image2(:,:,1)) + 0.63 * double(image2(:,:,2)) + 0.27 * double(image2(:,:,3));
M1 = 0.30 * double(image1(:,:,1)) + 0.04 * double(image1(:,:,2)) - 0.35 * double(image1(:,:,3));
M2 = 0.30 * double(image2(:,:,1)) + 0.04 * double(image2(:,:,2)) - 0.35 * double(image2(:,:,3));
N1 = 0.34 * double(image1(:,:,1)) - 0.60 * double(image1(:,:,2)) + 0.17 * double(image1(:,:,3));
N2 = 0.34 * double(image2(:,:,1)) - 0.60 * double(image2(:,:,2)) + 0.17 * double(image2(:,:,3));
%%%%%%%%%%%%%%%%%%%%%%%%%
% Downsample the image
%%%%%%%%%%%%%%%%%%%%%%%%%
% minDimension = min(rows,cols);
% F = max(1,round(minDimension / 256));
F = 2;
aveKernel = fspecial('average',F);

aveM1 = conv2(M1, aveKernel,'same');
aveM2 = conv2(M2, aveKernel,'same');
M1 = aveM1(1:F:rows,1:F:cols);
M2 = aveM2(1:F:rows,1:F:cols);

aveN1 = conv2(N1, aveKernel,'same');
aveN2 = conv2(N2, aveKernel,'same');
N1 = aveN1(1:F:rows,1:F:cols);
N2 = aveN2(1:F:rows,1:F:cols);

aveL1 = conv2(L1, aveKernel,'same');
aveL2 = conv2(L2, aveKernel,'same');
L1 = aveL1(1:F:rows,1:F:cols);
L2 = aveL2(1:F:rows,1:F:cols);

aveSM1 = conv2(saliencyMap1, aveKernel,'same');
aveSM2 = conv2(saliencyMap2, aveKernel,'same');
saliencyMap1 = aveSM1(1:F:rows,1:F:cols);
saliencyMap2 = aveSM2(1:F:rows,1:F:cols);
%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the gradient map
%%%%%%%%%%%%%%%%%%%%%%%%%
dx = [3 0 -3; 10 0 -10;  3  0 -3]/16;
dy = [3 10 3; 0  0   0; -3 -10 -3]/16;

IxL1 = conv2(L1, dx, 'same');     
IyL1 = conv2(L1, dy, 'same');    
gradientMap1 = sqrt(IxL1.^2 + IyL1.^2);

IxL2 = conv2(L2, dx, 'same');     
IyL2 = conv2(L2, dy, 'same');    
gradientMap2 = sqrt(IxL2.^2 + IyL2.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the VSI
%%%%%%%%%%%%%%%%%%%%%%%%%
% apply mask
mask_dwon = imresize(mask, size(saliencyMap1));
saliencyMap1 = saliencyMap1(mask_dwon);
saliencyMap2 = saliencyMap2(mask_dwon);
gradientMap1 = gradientMap1(mask_dwon);
gradientMap2 = gradientMap2(mask_dwon);
M1 = M1(mask_dwon);
M2 = M2(mask_dwon);
N1 = N1(mask_dwon);
N2 = N2(mask_dwon);

VSSimMatrix = (2 * saliencyMap1 .* saliencyMap2 + constForVS) ./ (saliencyMap1.^2 + saliencyMap2.^2 + constForVS);
gradientSimMatrix = (2*gradientMap1.*gradientMap2 + constForGM) ./(gradientMap1.^2 + gradientMap2.^2 + constForGM);

weight = max(saliencyMap1, saliencyMap2);

ISimMatrix = (2 * M1 .* M2 + constForChrom) ./ (M1.^2 + M2.^2 + constForChrom);
QSimMatrix = (2 * N1 .* N2 + constForChrom) ./ (N1.^2 + N2.^2 + constForChrom);

SimMatrixC = (gradientSimMatrix .^ alpha) .* VSSimMatrix .* real((ISimMatrix .* QSimMatrix) .^ lambda) .* weight;
sim = sum(sum(SimMatrixC)) / sum(weight(:));

return;
% %===================================
% function VSMap = SDSP(image,sigmaF,omega0,sigmaD,sigmaC)
% % ========================================================================
% % SDSP algorithm for salient region detection from a given image.
% % Copyright(c) 2013 Lin ZHANG, School of Software Engineering, Tongji
% % University
% % All Rights Reserved.
% % ----------------------------------------------------------------------
% % Permission to use, copy, or modify this software and its documentation
% % for educational and research purposes only and without fee is here
% % granted, provided that this copyright notice and the original authors'
% % names appear on all copies and supporting documentation. This program
% % shall not be used, rewritten, or adapted as the basis of a commercial
% % software or hardware product without first obtaining permission of the
% % authors. The authors make no representations about the suitability of
% % this software for any purpose. It is provided "as is" without express
% % or implied warranty.
% %----------------------------------------------------------------------
% %
% % This is an implementation of the algorithm for calculating the
% % SDSP (Saliency Detection by combining Simple Priors).
% %
% % Please refer to the following paper
% %
% % Lin Zhang, Zhongyi Gu, and Hongyu Li,"SDSP: a novel saliency detection 
% % method by combining simple priors", ICIP, 2013.
% % 
% %----------------------------------------------------------------------
% %
% %Input : image: an uint8 RGB image with dynamic range [0, 255] for each
% %channel
% %        
% %Output: VSMap: the visual saliency map extracted by the SDSP algorithm.
% %Data range for VSMap is [0, 255]. So, it can be regarded as a common
% %gray-scale image.
% %        
% %-----------------------------------------------------------------------
% %convert the image into LAB color space
% [oriRows, oriCols, junk] = size(image);
% image = double(image);
% dsImage(:,:,1) = imresize(image(:,:,1), [256, 256],'bilinear');
% dsImage(:,:,2) = imresize(image(:,:,2), [256, 256],'bilinear');
% dsImage(:,:,3) = imresize(image(:,:,3), [256, 256],'bilinear');
% lab = RGB2Lab(dsImage); 
% 
% LChannel = lab(:,:,1);
% AChannel = lab(:,:,2);
% BChannel = lab(:,:,3);
% 
% LFFT = fft2(double(LChannel));
% AFFT = fft2(double(AChannel));
% BFFT = fft2(double(BChannel));
% 
% [rows, cols, junk] = size(dsImage);
% LG = logGabor(rows,cols,omega0,sigmaF);
% FinalLResult = real(ifft2(LFFT.*LG));
% FinalAResult = real(ifft2(AFFT.*LG));
% FinalBResult = real(ifft2(BFFT.*LG));
% 
% SFMap = sqrt(FinalLResult.^2 + FinalAResult.^2 + FinalBResult.^2);
% 
% %the central areas will have a bias towards attention
% coordinateMtx = zeros(rows, cols, 2);
% coordinateMtx(:,:,1) = repmat((1:1:rows)', 1, cols);
% coordinateMtx(:,:,2) = repmat(1:1:cols, rows, 1);
% 
% centerY = rows / 2;
% centerX = cols / 2;
% centerMtx(:,:,1) = ones(rows, cols) * centerY;
% centerMtx(:,:,2) = ones(rows, cols) * centerX;
% SDMap = exp(-sum((coordinateMtx - centerMtx).^2,3) / sigmaD^2);
% 
% %warm colors have a bias towards attention
% maxA = max(AChannel(:));
% minA = min(AChannel(:));
% normalizedA = (AChannel - minA) / (maxA - minA);
% 
% maxB = max(BChannel(:));
% minB = min(BChannel(:));
% normalizedB = (BChannel - minB) / (maxB - minB);
% 
% labDistSquare = normalizedA.^2 + normalizedB.^2;
% SCMap = 1 - exp(-labDistSquare / (sigmaC^2));
% % VSMap = SFMap .* SDMap;
% VSMap = SFMap .* SDMap .* SCMap;
% 
% VSMap =  imresize(VSMap, [oriRows, oriCols],'bilinear');
% VSMap = mat2gray(VSMap);
% return;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function labImage = RGB2Lab(image)
% 
% image = double(image);
% normalizedR = image(:,:,1) / 255;
% normalizedG = image(:,:,2) / 255;
% normalizedB = image(:,:,3) / 255;
% 
% RSmallerOrEqualto4045 = normalizedR <= 0.04045;
% RGreaterThan4045 = 1 - RSmallerOrEqualto4045;
% tmpR = (normalizedR / 12.92) .* RSmallerOrEqualto4045;
% tmpR = tmpR + power((normalizedR + 0.055)/1.055,2.4) .* RGreaterThan4045;
% 
% GSmallerOrEqualto4045 = normalizedG <= 0.04045;
% GGreaterThan4045 = 1 - GSmallerOrEqualto4045;
% tmpG = (normalizedG / 12.92) .* GSmallerOrEqualto4045;
% tmpG = tmpG + power((normalizedG + 0.055)/1.055,2.4) .* GGreaterThan4045;
% 
% BSmallerOrEqualto4045 = normalizedB <= 0.04045;
% BGreaterThan4045 = 1 - BSmallerOrEqualto4045;
% tmpB = (normalizedB / 12.92) .* BSmallerOrEqualto4045;
% tmpB = tmpB + power((normalizedB + 0.055)/1.055,2.4) .* BGreaterThan4045;
% 
% X = tmpR*0.4124564 + tmpG*0.3575761 + tmpB*0.1804375;
% Y = tmpR*0.2126729 + tmpG*0.7151522 + tmpB*0.0721750;
% Z = tmpR*0.0193339 + tmpG*0.1191920 + tmpB*0.9503041;
% 
% epsilon = 0.008856;	%actual CIE standard
% kappa   = 903.3;		%actual CIE standard
%  
% Xr = 0.9642;	%reference white D50
% Yr = 1.0;		%reference white
% Zr = 0.8251;	%reference white
% 
% xr = X/Xr;
% yr = Y/Yr;
% zr = Z/Zr;
% 
% xrGreaterThanEpsilon = xr > epsilon;
% xrSmallerOrEqualtoEpsilon = 1 - xrGreaterThanEpsilon;
% fx = power(xr, 1.0/3.0) .* xrGreaterThanEpsilon;
% fx = fx + (kappa*xr + 16.0)/116.0 .* xrSmallerOrEqualtoEpsilon;
% 
% yrGreaterThanEpsilon = yr > epsilon;
% yrSmallerOrEqualtoEpsilon = 1 - yrGreaterThanEpsilon;
% fy = power(yr, 1.0/3.0) .* yrGreaterThanEpsilon;
% fy = fy + (kappa*yr + 16.0)/116.0 .* yrSmallerOrEqualtoEpsilon;
% 
% zrGreaterThanEpsilon = zr > epsilon;
% zrSmallerOrEqualtoEpsilon = 1 - zrGreaterThanEpsilon;
% fz = power(zr, 1.0/3.0) .* zrGreaterThanEpsilon;
% fz = fz + (kappa*zr + 16.0)/116.0 .* zrSmallerOrEqualtoEpsilon;
% 
% [rows,cols,junk] = size(image);
% labImage = zeros(rows,cols,3);
% labImage(:,:,1) = 116.0 * fy - 16.0;
% labImage(:,:,2) = 500.0 * (fx - fy);
% labImage(:,:,3) = 200.0 * (fy - fz);
% return;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function LG = logGabor(rows,cols,omega0,sigmaF)
%      [u1, u2] = meshgrid(([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
% 			            ([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2)));
%      mask = ones(rows, cols);
%      for rowIndex = 1:rows
%          for colIndex = 1:cols
%              if u1(rowIndex, colIndex)^2 + u2(rowIndex, colIndex)^2 > 0.25
%                  mask(rowIndex, colIndex) = 0;
%              end
%          end
%      end
%      u1 = u1 .* mask;
%      u2 = u2 .* mask;
%      
%      u1 = ifftshift(u1);  
%      u2 = ifftshift(u2);
%      
%      radius = sqrt(u1.^2 + u2.^2);    
%      radius(1,1) = 1;
%             
%      LG = exp((-(log(radius/omega0)).^2) / (2 * (sigmaF^2)));  
%      LG(1,1) = 0; 
% return;