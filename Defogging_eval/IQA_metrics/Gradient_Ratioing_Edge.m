function [e1,r1,ns1] = Gradient_Ratioing_Edge(fogImg, defogImg)

I1=double(fogImg);
R1=double(defogImg);

[nl,nc]=size(I1);

%%%% Sobel Gradient
Sy = double([1 2 1;0 0 0;-1 -2 -1]);
GIy = imfilter(I1,Sy,'symmetric');
GRy1 = imfilter(R1,Sy,'symmetric');

Sx = double([1 0 -1; 2 0 -2; 1 0 -1]);
GIx = imfilter(I1,Sx,'symmetric');
GRx1 = imfilter(R1,Sx,'symmetric');

GI=sqrt((GIx.^2)+(GIy.^2));
GR1=sqrt((GRx1.^2)+(GRy1.^2));

%%%% Contrast Computation at 5%
[C1, Crr1]=functionContrastAt5PerCent(R1);
[Ci, ~]=functionContrastAt5PerCent(I1);

%%%% Visible Gradients Ratio
Cratio1=zeros(nl,nc);
Cratio1(Crr1>0)=GR1(Crr1>0)./GI(Crr1>0);

%%%% Descriptor computation
% Rate of new visible edges
whitePixels1=sum(C1(:));
whitePixelsi=sum(Ci(:));
e1=(whitePixels1-whitePixelsi)/whitePixelsi;

% Number of saturated pixels after restoration
ns1=sum(R1(:)==256 |R1(:)==0);
ns1=ns1/(nl*nc);

% Restoration quality (geometric mean ratios of visibility level)
XX=log(Cratio1);
r1=exp((1/(whitePixels1))*nansum(XX(isfinite(XX))));

end