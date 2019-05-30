function score = realness_index(fogImg, clearImg, mask)

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

% calculate phase congruency
PC1 = phasecong2(Y1);
PC2 = phasecong2(Y2);

maxPCMap = max(PC1,PC2);

T1 = 0.85;  %fixed
PCSimMap = (2 * PC1 .* PC2 + T1) ./ (PC1.^2 + PC2.^2 + T1);

ChromSimMap = chromine_similarity_map(fogImg, clearImg, 'LMN');

% apply mask
mask_dn  = imresize(mask, size(maxPCMap));

alpha_3 = 1; alpha_4 = 0.02;
featureVector = real( PCSimMap(mask_dn).^alpha_3 .* real(ChromSimMap(mask_dn)).^alpha_4 .* maxPCMap(mask_dn) );

score = sum(featureVector)./sum(maxPCMap(mask_dn));

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ResultPC]=phasecong2(im)
% ========================================================================
% Copyright (c) 1996-2009 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby  granted, free of charge, to any  person obtaining a copy
% of this software and associated  documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% The software is provided "as is", without warranty of any kind.
% References:
%
%     Peter Kovesi, "Image Features From Phase Congruency". Videre: A
%     Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
%     Summer 1999 http://mitpress.mit.edu/e-journals/Videre/001/v13.html

nscale          = 4;     % Number of wavelet scales.    
norient         = 4;     % Number of filter orientations.
minWaveLength   = 6;     % Wavelength of smallest scale filter.    
mult            = 2;   % Scaling factor between successive filters.    
sigmaOnf        = 0.55;  % Ratio of the standard deviation of the
                             % Gaussian describing the log Gabor filter's
                             % transfer function in the frequency domain
                             % to the filter center frequency.    
dThetaOnSigma   = 1.2;   % Ratio of angular interval between filter orientations    
                             % and the standard deviation of the angular Gaussian
                             % function used to construct filters in the
                             % freq. plane.
k               = 2.0;   % No of standard deviations of the noise
                             % energy beyond the mean at which we set the
                             % noise threshold point. 
                             % below which phase congruency values get
                             % penalized.
epsilon         = .0001;                % Used to prevent division by zero.

thetaSigma = pi/norient/dThetaOnSigma;  % Calculate the standard deviation of the
                                        % angular Gaussian function used to
                                        % construct filters in the freq. plane.

[rows,cols] = size(im);
imagefft = fft2(im);              % Fourier transform of image

zero = zeros(rows,cols);
EO = cell(nscale, norient);       % Array of convolution results.                                 

estMeanE2n = [];
ifftFilterArray = cell(1,nscale); % Array of inverse FFTs of filters

% Pre-compute some stuff to speed up filter construction

% Set up X and Y matrices with ranges normalised to +/- 0.5
% The following code adjusts things appropriately for odd and even values
% of rows and columns.
if mod(cols,2)
    xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
else
    xrange = [-cols/2:(cols/2-1)]/cols;	
end

if mod(rows,2)
    yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
else
    yrange = [-rows/2:(rows/2-1)]/rows;	
end

[x,y] = meshgrid(xrange, yrange);

radius = sqrt(x.^2 + y.^2);       % Matrix values contain *normalised* radius from centre.
theta = atan2(-y,x);              % Matrix values contain polar angle.
                                  % (note -ve y is used to give +ve
                                  % anti-clockwise angles)
				  
radius = ifftshift(radius);       % Quadrant shift radius and theta so that filters
theta  = ifftshift(theta);        % are constructed with 0 frequency at the corners.
radius(1,1) = 1;                  % Get rid of the 0 radius value at the 0
                                  % frequency point (now at top-left corner)
                                  % so that taking the log of the radius will 
                                  % not cause trouble.

sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;    % save a little memory

% Filters are constructed in terms of two components.
% 1) The radial component, which controls the frequency band that the filter
%    responds to
% 2) The angular component, which controls the orientation that the filter
%    responds to.
% The two components are multiplied together to construct the overall filter.

% Construct the radial filter components...

% First construct a low-pass filter that is as large as possible, yet falls
% away to zero at the boundaries.  All log Gabor filters are multiplied by
% this to ensure no extra frequencies at the 'corners' of the FFT are
% incorporated as this seems to upset the normalisation process when
% calculating phase congrunecy.
lp = lowpassfilter([rows,cols],.45,15);   % Radius .45, 'sharpness' 15

logGabor = cell(1,nscale);

for s = 1:nscale
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor{s} = logGabor{s}.*lp;        % Apply low-pass filter
    logGabor{s}(1,1) = 0;                 % Set the value at the 0 frequency point of the filter
                                          % back to zero (undo the radius fudge).
end

% Then construct the angular filter components...

spread = cell(1,norient);

for o = 1:norient
  angl = (o-1)*pi/norient;           % Filter angle.

  % For each point in the filter matrix calculate the angular distance from
  % the specified filter orientation.  To overcome the angular wrap-around
  % problem sine difference and cosine difference values are first computed
  % and then the atan2 function is used to determine angular distance.

  ds = sintheta * cos(angl) - costheta * sin(angl);    % Difference in sine.
  dc = costheta * cos(angl) + sintheta * sin(angl);    % Difference in cosine.
  dtheta = abs(atan2(ds,dc));                          % Absolute angular distance.
  spread{o} = exp((-dtheta.^2) / (2 * thetaSigma^2));  % Calculate the
                                                       % angular filter component.
end

% The main loop...
EnergyAll(rows,cols) = 0;
AnAll(rows,cols) = 0;

for o = 1:norient                    % For each orientation.
  sumE_ThisOrient   = zero;          % Initialize accumulator matrices.
  sumO_ThisOrient   = zero;       
  sumAn_ThisOrient  = zero;      
  Energy            = zero;      
  for s = 1:nscale,                  % For each scale.
    filter = logGabor{s} .* spread{o};   % Multiply radial and angular
                                         % components to get the filter. 
    ifftFilt = real(ifft2(filter))*sqrt(rows*cols);  % Note rescaling to match power
    ifftFilterArray{s} = ifftFilt;                   % record ifft2 of filter
    % Convolve image with even and odd filters returning the result in EO
    EO{s,o} = ifft2(imagefft .* filter);      

    An = abs(EO{s,o});                         % Amplitude of even & odd filter response.
    sumAn_ThisOrient = sumAn_ThisOrient + An;  % Sum of amplitude responses.
    sumE_ThisOrient = sumE_ThisOrient + real(EO{s,o}); % Sum of even filter convolution results.
    sumO_ThisOrient = sumO_ThisOrient + imag(EO{s,o}); % Sum of odd filter convolution results.
    if s==1                                 % Record mean squared filter value at smallest
      EM_n = sum(sum(filter.^2));           % scale. This is used for noise estimation.
      maxAn = An;                           % Record the maximum An over all scales.
    else
      maxAn = max(maxAn, An);
    end
  end                                       % ... and process the next scale

  % Get weighted mean filter response vector, this gives the weighted mean
  % phase angle.

  XEnergy = sqrt(sumE_ThisOrient.^2 + sumO_ThisOrient.^2) + epsilon;   
  MeanE = sumE_ThisOrient ./ XEnergy; 
  MeanO = sumO_ThisOrient ./ XEnergy; 

  % Now calculate An(cos(phase_deviation) - | sin(phase_deviation)) | by
  % using dot and cross products between the weighted mean filter response
  % vector and the individual filter response vectors at each scale.  This
  % quantity is phase congruency multiplied by An, which we call energy.

  for s = 1:nscale,       
      E = real(EO{s,o}); O = imag(EO{s,o});    % Extract even and odd
                                               % convolution results.
      Energy = Energy + E.*MeanE + O.*MeanO - abs(E.*MeanO - O.*MeanE);
  end

  % Compensate for noise
  % We estimate the noise power from the energy squared response at the
  % smallest scale.  If the noise is Gaussian the energy squared will have a
  % Chi-squared 2DOF pdf.  We calculate the median energy squared response
  % as this is a robust statistic.  From this we estimate the mean.
  % The estimate of noise power is obtained by dividing the mean squared
  % energy value by the mean squared filter value

  medianE2n = median(reshape(abs(EO{1,o}).^2,1,rows*cols));
  meanE2n = -medianE2n/log(0.5);
  estMeanE2n(o) = meanE2n;

  noisePower = meanE2n/EM_n;                       % Estimate of noise power.

  % Now estimate the total energy^2 due to noise
  % Estimate for sum(An^2) + sum(Ai.*Aj.*(cphi.*cphj + sphi.*sphj))

  EstSumAn2 = zero;
  for s = 1:nscale
    EstSumAn2 = EstSumAn2 + ifftFilterArray{s}.^2;
  end

  EstSumAiAj = zero;
  for si = 1:(nscale-1)
    for sj = (si+1):nscale
      EstSumAiAj = EstSumAiAj + ifftFilterArray{si}.*ifftFilterArray{sj};
    end
  end
  sumEstSumAn2 = sum(sum(EstSumAn2));
  sumEstSumAiAj = sum(sum(EstSumAiAj));

  EstNoiseEnergy2 = 2*noisePower*sumEstSumAn2 + 4*noisePower*sumEstSumAiAj;

  tau = sqrt(EstNoiseEnergy2/2);                     % Rayleigh parameter
  EstNoiseEnergy = tau*sqrt(pi/2);                   % Expected value of noise energy
  EstNoiseEnergySigma = sqrt( (2-pi/2)*tau^2 );

  T =  EstNoiseEnergy + k*EstNoiseEnergySigma;       % Noise threshold

  % The estimated noise effect calculated above is only valid for the PC_1 measure. 
  % The PC_2 measure does not lend itself readily to the same analysis.  However
  % empirically it seems that the noise effect is overestimated roughly by a factor 
  % of 1.7 for the filter parameters used here.

  T = T/1.7;        % Empirical rescaling of the estimated noise effect to 
                    % suit the PC_2 phase congruency measure
  Energy = max(Energy - T, zero);          % Apply noise threshold

  EnergyAll = EnergyAll + Energy;
  AnAll = AnAll + sumAn_ThisOrient;
end  % For each orientation
ResultPC = EnergyAll ./ AnAll;
return;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOWPASSFILTER - Constructs a low-pass butterworth filter.
%
% usage: f = lowpassfilter(sze, cutoff, n)
% 
% where: sze    is a two element vector specifying the size of filter 
%               to construct [rows cols].
%        cutoff is the cutoff frequency of the filter 0 - 0.5
%        n      is the order of the filter, the higher n is the sharper
%               the transition is. (n must be an integer >= 1).
%               Note that n is doubled so that it is always an even integer.
%
%                      1
%      f =    --------------------
%                              2n
%              1.0 + (w/cutoff)
%
% The frequency origin of the returned filter is at the corners.
%
% See also: HIGHPASSFILTER, HIGHBOOSTFILTER, BANDPASSFILTER
%

% Copyright (c) 1999 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% October 1999
% August  2005 - Fixed up frequency ranges for odd and even sized filters
%                (previous code was a bit approximate)

function f = lowpassfilter(sze, cutoff, n)
    
    if cutoff < 0 || cutoff > 0.5
	error('cutoff frequency must be between 0 and 0.5');
    end
    
    if rem(n,1) ~= 0 || n < 1
	error('n must be an integer >= 1');
    end

    if length(sze) == 1
	rows = sze; cols = sze;
    else
	rows = sze(1); cols = sze(2);
    end

    % Set up X and Y matrices with ranges normalised to +/- 0.5
    % The following code adjusts things appropriately for odd and even values
    % of rows and columns.
    if mod(cols,2)
	xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
    else
	xrange = [-cols/2:(cols/2-1)]/cols;	
    end

    if mod(rows,2)
	yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
    else
	yrange = [-rows/2:(rows/2-1)]/rows;	
    end
    
    [x,y] = meshgrid(xrange, yrange);
    radius = sqrt(x.^2 + y.^2);        % A matrix with every pixel = radius relative to centre.
    f = ifftshift( 1 ./ (1.0 + (radius ./ cutoff).^(2*n)) );   % The filter
    return;
end



