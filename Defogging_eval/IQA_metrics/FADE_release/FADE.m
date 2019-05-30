function [D D_map] = FADE(I)
    % Input: a test image, I
    % Output: perceptual fog density D and fog density map D_map
    % Detail explanation:
    % L. K. Choi, J. You, and A. C. Bovik, "Referenceless Prediction of Perceptual Fog Density and Perceptual Image Defogging",
    % IEEE Transactions on Image Processing, to appear (2015).
                                         
    %% Basic setup        
            ps                      = 8;                                            % patch size 8 x 8 pixels
        % Size of a test image for checking possilbe distinct patches 
            [row col dim]           = size(I);
            patch_row_num           = floor(row/ps);
            patch_col_num           = floor(col/ps);
            I                       = I(1:patch_row_num*ps,1:patch_col_num*ps,1:3);           
            [row col dim]           = size(I);
            patch_row_num           = floor(row/ps);
            patch_col_num           = floor(col/ps);
            I                       = I(1:patch_row_num*ps,1:patch_col_num*ps,1:3); 
        % RGB and gray channel
            R                       = double(I(:,:,1));                             % Red
            G                       = double(I(:,:,2));                             % Green
            B                       = double(I(:,:,3));                             % Blue
            Ig                      = double(rgb2gray(I));                          % Gray        
        % Dark channel prior image: Id, pixel-wise, scaled to [0 1]     
            Irn                     = R./255;
            Ign                     = G./255;
            Ibn                     = B./255;
            Id                      = min(min(Irn,Ign),Ibn);
        % HSV color space: saturation image: Is                
            I_hsv                   = rgb2hsv(I);
            Is                      = I_hsv(:,:,2);                
        % MSCN
            MSCN_window             = fspecial('gaussian',7,7/6);
            MSCN_window             = MSCN_window/sum(sum(MSCN_window));        
            warning('off');               
            mu                      = imfilter(Ig,MSCN_window,'replicate');
            mu_sq                   = mu.*mu;
            sigma                   = sqrt(abs(imfilter(Ig.*Ig,MSCN_window,'replicate') - mu_sq));
            MSCN                    = (Ig-mu)./(sigma+1);
            cv                      = sigma./mu; 
        % rg and by channel
            rg                      = R-G;
            by                      = 0.5*(R+G)-B;                 

    %% Fog aware statistical feature extraction     
        % f1        
            MSCN_var                = reshape(nanvar(im2col(MSCN, [ps ps], 'distinct')),[row col]/ps);        
        % f2,f3                    
            MSCN_V_pair_col         = im2col((MSCN.*circshift(MSCN,[1 0])),[ps ps], 'distinct'); %vertical 
            MSCN_V_pair_col_temp1   = MSCN_V_pair_col; MSCN_V_pair_col_temp1(MSCN_V_pair_col_temp1>0)=NaN;
            MSCN_V_pair_col_temp2   = MSCN_V_pair_col; MSCN_V_pair_col_temp2(MSCN_V_pair_col_temp2<0)=NaN;
            MSCN_V_pair_L_var       = reshape(nanvar(MSCN_V_pair_col_temp1),[row col]/ps);          
            MSCN_V_pair_R_var       = reshape(nanvar(MSCN_V_pair_col_temp2),[row col]/ps);          
        % f4        
            Mean_sigma              = reshape(mean(im2col(sigma, [ps ps], 'distinct')),[row col]/ps);    
        % f5        
            Mean_cv                 = reshape(mean(im2col(cv, [ps ps], 'distinct')),[row col]/ps);
        % f6, f7, f8                    
            [CE_gray CE_by CE_rg]   = CE(I);
            Mean_CE_gray            = reshape(mean(im2col(CE_gray, [ps ps], 'distinct')),[row col]/ps);
            Mean_CE_by              = reshape(mean(im2col(CE_by, [ps ps], 'distinct')),[row col]/ps);
            Mean_CE_rg              = reshape(mean(im2col(CE_rg, [ps ps], 'distinct')),[row col]/ps);     
        % f9        
            IE_temp                 = num2cell(im2col(uint8(Ig), [ps ps], 'distinct'),1);
            IE                      = reshape(cellfun(@entropy,IE_temp),[row col]/ps);
        % f10        
            Mean_Id                 = reshape(mean(im2col(Id, [ps ps], 'distinct')),[row col]/ps);
        % f11        
            Mean_Is                 = reshape(mean(im2col(Is, [ps ps], 'distinct')),[row col]/ps);
        % f12        
            CF                      = reshape(sqrt(std(im2col(rg, [ps ps], 'distinct')).^2 + std(im2col(by, [ps ps], 'distinct')).^2) + 0.3* sqrt(mean(im2col(rg, [ps ps], 'distinct')).^2 + mean(im2col(by, [ps ps], 'distinct')).^2),[row col]/ps);                                
        feat                        = [MSCN_var(:) MSCN_V_pair_R_var(:) MSCN_V_pair_L_var(:) Mean_sigma(:) Mean_cv(:) Mean_CE_gray(:) Mean_CE_by(:) Mean_CE_rg(:) IE(:) Mean_Id(:) Mean_Is(:) CF(:)]; 
        feat                        = log(1 + feat);

    %% MVG model distance                 
        %Df (foggy level distance) for each patch
        % load natural fogfree image features (mu, cov)
            load('natural_fogfree_image_features_ps8.mat');        
        % test param for each patch                
            mu_fog_param_patch      = feat;
            cov_fog_param_patch     = nanvar(feat')';
        % Distance calculation - includes intermediate steps
            feature_size            = size(feat,2);
            mu_matrix               = repmat(mu_fogfreeparam, [size(feat,1),1]) - mu_fog_param_patch;         
            cov_temp1               = [];
            cov_temp1(cumsum(feature_size.*ones(1,length(cov_fog_param_patch))))=1;
            cov_temp2               = cov_fog_param_patch(cumsum(cov_temp1)-cov_temp1+1,:);
            cov_temp3               = repmat(cov_temp2, [1,feature_size]);
            cov_temp4               = repmat(cov_fogfreeparam,[length(cov_fog_param_patch),1]);
            cov_matrix              = (cov_temp3 + cov_temp4)/2;     
        % cell computation
            mu_cell                 = num2cell(mu_matrix,2);
            cov_cell                = mat2cell(cov_matrix, feature_size*ones(1, size(mu_matrix,1)),feature_size);        
            mu_transpose_cell       = num2cell(mu_matrix',1);     
        % foggy level computation
            distance_patch          = sqrt(cell2mat(cellfun(@mtimes,cellfun(@mrdivide,mu_cell,cov_cell, 'UniformOutput',0),mu_transpose_cell', 'UniformOutput',0)));          
            Df                      = nanmean(distance_patch);   % Mean_distance_patch
            Df_map                  = reshape(distance_patch,[row,col]/ps);                                        
            clear mu_matrix cov_matrix mu_cell cov_cell mu_transpose_cell distance_patch

        %Dff
        % load natural foggy image features (mu, cov)
            load('natural_foggy_image_features_ps8.mat');            
        % calculation of distance - includes intermediate steps            
            mu_matrix               = repmat(mu_foggyparam,  [size(feat,1),1]) - mu_fog_param_patch;                    
            cov_temp5               = repmat(cov_foggyparam,[length(cov_fog_param_patch),1]);               
            cov_matrix              = (cov_temp3 + cov_temp5)/2;          
        % cell computation
            mu_cell                 = num2cell(mu_matrix,2);              
            cov_cell                = mat2cell(cov_matrix, feature_size*ones(1, size(mu_matrix,1)),feature_size);                        
            mu_transpose_cell       = num2cell(mu_matrix',1);        
        % fog-free level computation
            distance_patch          = sqrt(cell2mat(cellfun(@mtimes,cellfun(@mrdivide,mu_cell,cov_cell, 'UniformOutput',0),mu_transpose_cell', 'UniformOutput',0)));     
            Dff                     = nanmean(distance_patch);   % Mean_distance_patch                
            Dff_map                 = reshape(distance_patch,[row,col]/ps); 
            clear mu_matrix cov_matrix mu_cell cov_cell mu_transpose_cell 
    %% Perceptual fog density and density map
        D                           = Df/(Dff+1);                
        D_map                       = Df_map./(Dff_map+1);                
end

