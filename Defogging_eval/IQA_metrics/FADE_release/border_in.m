function [nI]=border_in(I,ps)

% Input  - input image I, patch size ps
% Output - border added image

     if mod(ps,2)==0
         uc     = ps/2;      % upperside copy 
         dc     = ps/2-1;    % downside copy       
     elseif mod(ps,2)==1
         uc     = floor(ps/2);
         dc     = uc;
     end
     
     ucb        = I(1:uc,:);
     dcb        = I(end-dc:end,:);
     Igtemp1    = [ucb;I;dcb];

     lcb        = Igtemp1(:,1:uc);
     rcb        = Igtemp1(:,end-dc:end);
     nI         = [lcb Igtemp1 rcb]; 
      
end
       
        
        
        