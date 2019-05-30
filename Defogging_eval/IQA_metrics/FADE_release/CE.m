function [CE_gray CE_by CE_rg] = CE(I)

% Input 
    % I: input image    
% Output
    % Perceived Contrast Energy for gray, blue-yellow, red-green color channels

% parameters    
    % In this code, I followed Groen's experiment [28]: fixed sigma, weight, t1, t2,and t3    
    % sigma: 0.16 degree =~ 3.25pixel, filter windows =~ 20 x 20 = (-3.25*3:1:3.25*3)
    % semisaturation: 0.1 
    % t1, t2, t3 follows the EmpiricalNoiseThreshold used in [28]
    

    % Basic parameters
    sigma               = 3.25;
    semisaturation      = 0.1; 
    t1                  = 9.225496406318721e-004 *255; %0.2353; %
    t2                  = 8.969246659629488e-004 *255; %0.2287; %
    t3                  = 2.069284034165411e-004 *255; %0.0528; %
    border_s            = 20;

    % Gaussian & LoG & Retification, Normalization(?)
    break_off_sigma     = 3;
    filtersize          = break_off_sigma*sigma;
    x                   = -filtersize:1:filtersize;
    Gauss               = 1/(sqrt(2 * pi) * sigma)* exp((x.^2)/(-2 * sigma * sigma) );
    Gauss               = Gauss./(sum(Gauss));
    Gx                  = (x.^2/sigma^4-1/sigma^2).*Gauss;  %LoG
    Gx                  = Gx-sum(Gx)/size(x,2);
    Gx                  = Gx/sum(0.5*x.*x.*Gx);

    % Color conversion
    I                 = double(I);
    R                   = I(:,:,1);
    G                   = I(:,:,2);
    B                   = I(:,:,3);
    gray                = 0.299*R + 0.587*G + 0.114*B;
    by                  = 0.5*R + 0.5*G - B  ;
    rg                  = R - G;    
    [row col dim]       = size(I);
    CE_gray             = double(zeros(row, col));
    CE_by               = double(zeros(row, col));
    CE_rg               = double(zeros(row, col));
    
    % CE_Gray    
    gray_temp1          = border_in(gray,border_s);
    Cx_gray             = conv2(gray_temp1, Gx, 'same');
    Cy_gray             = conv2(gray_temp1, Gx','same');
    C_gray_temp2        = sqrt(Cx_gray.^2+Cy_gray.^2);
    C_gray              = border_out(C_gray_temp2,border_s);        
    R_gray              = (C_gray.*max(C_gray(:)))./(C_gray+(max(C_gray(:))*semisaturation));
    R_gray_temp1        = R_gray - t1;
    index1              = find(R_gray_temp1 > 0.0000001);
    CE_gray(index1)     = R_gray_temp1(index1);
    
    % CE_by    
    by_temp1            = border_in(by,border_s);
    Cx_by               = conv2(by_temp1, Gx, 'same');
    Cy_by               = conv2(by_temp1, Gx','same');
    C_by_temp2          = sqrt(Cx_by.^2+Cy_by.^2);
    C_by                = border_out(C_by_temp2,border_s);        
    R_by                = (C_by.*max(C_by(:)))./(C_by+(max(C_by(:))*semisaturation));
    R_by_temp1          = R_by - t2;
    index2              = find(R_by_temp1 > 0.0000001);
    CE_by(index2)       = R_by_temp1(index2);
    
    % CE_rg   
    rg_temp1            = border_in(rg,border_s);
    Cx_rg               = conv2(rg_temp1, Gx, 'same');
    Cy_rg               = conv2(rg_temp1, Gx','same');
    C_rg_temp2          = sqrt(Cx_rg.^2+Cy_rg.^2);
    C_rg                = border_out(C_rg_temp2,border_s);    
    R_rg                = (C_rg.*max(C_rg(:)))./(C_rg+(max(C_rg(:))*semisaturation));
    R_rg_temp1          = R_rg - t3;
    index3              = find(R_rg_temp1 > 0.0000001);
    CE_rg(index3)       = R_rg_temp1(index3);