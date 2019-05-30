function [nI]=border_out(I,ps)

% Input - input image I, patch size ps
% Output - border trimmed image

     if mod(ps,2)==0
         uc             = ps/2;      % upperside copy 
         dc             = ps/2-1;    % downside copy
       
     elseif mod(ps,2)==1
         uc             = floor(ps/2);
         dc             = uc;
     end
     
     I(:,1:uc)          = [];
     I(:,end-dc:end)    = [];
     I(1:uc,:)          = [];
     I(end-dc:end,:)    = []; 
     
     nI=I;
end
       
        
        
        