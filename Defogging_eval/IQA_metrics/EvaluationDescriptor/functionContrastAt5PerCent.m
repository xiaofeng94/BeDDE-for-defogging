% 23/02/2011
% IFSTTAR copyright
%
% The approach is described in details in 
%
% "Blind Contrast Restoration Assessment by Gradient Ratioing at Visible Edges",
% by N. Hautiere, J.-P. Tarel, D. Aubert and E. Dumont,
% in proceedings of International Congress for Stereology (ICS'07), 
% Saint Etienne, France, August 30-September 7, 2007.
% http://perso.lcpc.fr/tarel.jean-philippe/publis/ics07.html
%
function [Mask Crr]=functionContrastAt5PerCent(I1,S,percentage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I1 : Original image
% S : SubWindow size (default = 7)
% perCentage : Visibility percentage (default = 5)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mask : Mask of the visibility contrast > percentage
% Crr : valeur du contraste visibile (>percentage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
    error('not enough argument')
end
if nargin<2
    S=7;
end
if nargin<3
    percentage=5;
end

[nl,nc,dd]=size(I1);

%%% 4-neighborhood pixels
seh=[0 1 0;0 1 0;0 0 0];
seg=[0 0 0;1 1 0;0 0 0];

%%% minimum and maximum calculation
I1pad= padarray(I1,[3 3],'symmetric');

Igmin11=imerode(I1pad,seg,'full');
Ihmin11=imerode(I1pad,seh,'full');
Igmax11=imdilate(I1pad,seg,'full');
Ihmax11=imdilate(I1pad,seh,'full');

Igmin1=Igmin11(3+(1:nl),3+(1:nc));
Ihmin1=Ihmin11(3+(1:nl),3+(1:nc));
Igmax1=Igmax11(3+(1:nl),3+(1:nc));
Ihmax1=Ihmax11(3+(1:nl),3+(1:nc));

%%% pads image and minimum and maximum matrix
I1pad= padarray(I1,[S S],'symmetric');

Igmin1pad= padarray(Igmin1,[S S],'symmetric');
Ihmin1pad= padarray(Ihmin1,[S S],'symmetric');
Igmax1pad= padarray(Igmax1,[S S],'symmetric');
Ihmax1pad= padarray(Ihmax1,[S S],'symmetric');

%%% Initialization
Is=zeros(S,S);
Mask=false(nl+2*S,nc+2*S);
Crr=zeros(nl+2*S,nc+2*S);
s=1;
percentage=percentage/2;

% h = waitbar(0,'Please wait...');

%%% subwindow of size S*S
for ii=1:round(S/2):nl
    for jj=1:round(S/2):nc

        Is=double(I1pad(S+ii:2*S+ii-1,S+jj:2*S+jj-1));

        Isgmin=double(Igmin1pad(S+ii:2*S+ii-1,S+jj:2*S+jj-1));
        Ishmin=double(Ihmin1pad(S+ii:2*S+ii-1,S+jj:2*S+jj-1));
        Isgmax=double(Igmax1pad(S+ii:2*S+ii-1,S+jj:2*S+jj-1));
        Ishmax=double(Ihmax1pad(S+ii:2*S+ii-1,S+jj:2*S+jj-1));

        %%% horizontal and vertical contrasts
        Cgxx1=zeros(1,S^2);
        Chxx1=zeros(1,S^2);

        Ismin=round(min(Is(:)));
        Ismax=round(max(Is(:)));

        if(Ismin<=0) Ismin=1;end
        if(Ismin>256) Ismin=256;end
        if(Ismax<=0) Ismax=1; end
        if(Ismax>256) Ismax=256;end

        Fcube=false(S,S,Ismax-Ismin+1);

        C=zeros(1,Ismax);

        G=zeros(1,Ismax);

        for s=Ismin:Ismax  %%% we vary the threshold

            Fg=0;   %%% Pixels cardinal separated by s
            pg=1;   %%% Contrast indice
            Fh=0;
            ph=1;

            Cgxx1=zeros(1,S^2);
            Chxx1=zeros(1,S^2);

            for nn=2:S
                for mm=2:S
                    %%% Vertical contrast
                    if((Isgmin(nn,mm)<=s) && (Isgmax(nn,mm)>s))%%% test if two pixels are separated by s

                        %%%% Weber Contrast calculation
                        Cgxx1(pg)=min( abs(s-Is(nn,mm))/(max(s,Is(nn,mm))), abs(s-Is(nn,mm-1))/(max(s,Is(nn,mm-1))));

                        pg=pg+1;
                        Fg=Fg+1;

                        %% Mask calculation
                        Fcube(nn,mm,s-Ismin+1)=true;
                        Fcube(nn,mm-1,s-Ismin+1)=true;

                    end
                    %%% Horizontal contrast
                    if((Ishmin(nn,mm)<=s) && (Ishmax(nn,mm)>s))

                        %%%% Weber Contrast calculation
                        Chxx1(ph)=min( abs(s-Is(nn,mm))/(max(s,Is(nn,mm))) , abs(s-Is(nn-1,mm))/(max(s,Is(nn-1,mm))) );

                        ph=ph+1;
                        Fh=Fh+1;

                        %% Mask calculation
                        Fcube(nn,mm,s-Ismin+1)=true;
                        Fcube(nn-1,mm,s-Ismin+1)=true;

                    end

                end
            end

            if ((Fg+Fh) > 0)  %%% we verified if there is a threshold s that separates at least two pixels

                C(s)=(1/(Fg+Fh))*(sum(Cgxx1)+sum(Chxx1));   %% Contrast of threshold s

            end

        end
        %%% threshold that maximize the contrast
        [M s0] = max(C);
        s0=max(s0,Ismin);
        M=256*M;


        if(M>(256*(percentage/100)))
            %%% overlapping : we make a "logical |" because we evaluate contrast in windows of size S moving with a step of S/2
            
            Mask(S+ii:2*S+ii-1,S+jj:2*S+jj-1) =  Mask(S+ii:2*S+ii-1,S+jj:2*S+jj-1) |  Fcube(:,:,s0-Ismin+1);

            I3=Mask(S+ii:2*S+ii-1,S+jj:2*S+jj-1);

            %%% recover the contrast value of visible edge
            Crr1=zeros(S,S);

            Crr1(I3)=2*M/256;

            Crr(S+ii:2*S+ii-1,S+jj:2*S+jj-1)=Crr1;

        end
               
    end

%     waitbar(ii/nl,h)

end

Mask=Mask(S+1:nl+S,S+1:nc+S);
Crr=Crr(S+1:nl+S,S+1:nc+S);

% close(h);
