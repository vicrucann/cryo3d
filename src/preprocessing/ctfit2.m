function [P c]=ctfit2(mc,Pa,pixA,maxres,minres,test)
% function [P c]=ctfit2(mc,Pa,pixA,maxres,minres,test)
%  Derived from ppcfit2, uses the algorithm of N. Grigorieff's CTFFIND3.
% Determine the CTF parameters P from the image mc.  The cross-correlation
% coefficient c is returned too.
% Positive defocus values correspond to *underfocus*.
% Parameter ranges are defined by the structure Pa.  If vectors of values
% are given, all combinations of these values are tested in the brute-force
% search.  Afterwards, parameters are optimized using Simplex.
% The arguments pixA, maxres, minres are all in angstroms.
% The optional test argument controls display of results.  test=1 means graphical
% output (default).  test=2 gives text display of execution times.
% Here is an example:
%   Pa.lambda=EWavelength(200); % lambda in angstroms
%   Pa.defocus=0:0.5:20;  % in microns
%   Pa.deltadef=-2:.4:2;  % astigmatism in um
%   Pa.theta=0;           % astigmatism angle: all values are automatically
%   searched.
%   Pa.alpha=.07;  % If a range is given, this can also be fitted, for use with phase plate
%    % The following parameters are not fitted
%   Pa.Cs=2;
%   Pa.B=1000;  % This is not fitted
%    % Actual function call
%   P=ctfit2(m,Pa,1.4,8,50);
%   P      % print out the values

% Redundancy in references due to overfocus should be removed.

if nargin<6
    test=1;  % show graphics at least.
end;
NumRefs=2;
persistent R;  % We cache the references to allow faster execution of repeats.
overfocus=0;

mc=single(mc);
mc=mc-mean(mc(:));
[nx ny]=size(mc);
%
% maxres=6;
% minres=40;

res=pixA; % pixel size in A.
Pa.res=res;

if nx>5000
    nu=256;  % FFT block size, must be a multiple of 4.
else
    nu=128;
end;

w0=100;  % effective res of filter, in A.
disexp=0.2; % display exponent
pwExp1=4;  % pre-whitening positive exponent
pwExp4=2; % pre-whitening negative exponent (r^4)

w=nu*res/w0;  % filter radius.  Niko has it = (1/40A) in freq. space.

df=1/(res*nu);  % spatial frequency unit
freqs=(0:df:df*(nu/2-1))';

% Compute local power spectra
[nx ny]=size(mc);
nv=round(nu/4);   % overlap
tx=TileCoords(nx,nu,nv);  % no. of tiles in X direction
ty=TileCoords(ny,nu,nv);  % no. of tiles in Y direction

if test>0  % Show the micrograph
    figure(1);
    SetGrayscale;
    subplot(2,3,1);
    imacs(mc);
    LabelImage(mc,res/10,'nm');
    drawnow;
end;

sd=zeros(tx,ty);
% Make a window for spectral estiation
window=fuzzymask(nu,2,0.45*nu,0.1*nu);
winsum=sum(window(:));

sps=zeros(nu,nu,tx,ty);

if test>1
    disp('Computing spectra...');
    tic
end;
for ix=1:tx
    [x0 x1 u1]=TileCoords(nx,nu,nv,ix);
    
    for iy=1:ty
        [y0 y1 v1]=TileCoords(ny,nu,nv,iy);
        % I don't think it's necessary to remove gradients, but Niko does
        % this:
        tm=RemoveGradients(double(mc(x0+1:x0+nu,y0+1:y0+nu)));
        tm=tm.*window;
        tm=tm-sum(tm(:))*window/winsum;
        
        % subplot(2,1,1); imacs(tm); title([ix iy]);
        
        sp2=abs(fftn(tm)).^2;
        % subplot(2,1,2); imacs(fftshift(sqrt(sp2)));
        % drawnow;
        sps(:,:,ix,iy)=sp2;
        sd(ix,iy)=sqrt(sum(sp2(:)))/nu;
    end;
end;
if test>1
    toc
end;
% subplot(2,3,1);
% sumsp=sum(sum(sps,4),3);
% imacs(fftshift(sumsp.^disexp));
% LabelImage(nu,df,'A^{-1}');


% Make a histogram of the s.d. values and find tiles having s.d. near the
% mode.
nbins=tx+2;
hthresh=0.2;
% hthresh=-1;  % no selection

sdmn=min(sd(:));
sdmx=max(sd(:));
dbin=(sdmx-sdmn)/(nbins-1);
bins=sdmn:dbin:sdmx;

[h x]=histc(sd(:),bins);

[mxh j]=max(h);
% j is the mode of the distribution
% Find the first bin below hthresh on the high side
jmax=j-1+find(h(j:nbins)<hthresh*mxh,1,'first');
if numel(jmax)<1
    jmax=nbins;
end;
sdmax=x(jmax);
% find the bin below hthresh on the low side
jmin=find(h(1:j)<hthresh*mxh,1,'last');
if numel(jmin)<1
    jmin=1;
end;
sdmin=bins(jmin);
sdmax=bins(jmax)+dbin;

% Show histogram
% subplot(2,3,2);
% bar(x,h);
% hold on
% plot(sdmin,0,'w.');
% plot(sdmax,0,'w.');
% hold off

% Show regions that we used.
if test
    subplot(2,3,2);
    imacs((sd<=sdmax).*(sd>=sdmin));
    drawnow;
end;

cumsp=zeros(nu,nu);
count=0;
for ix=1:tx
    for iy=1:ty
        if (sd(ix,iy)<=sdmax) && (sd(ix,iy)>=sdmin)
            cumsp=cumsp+sps(:,:,ix,iy);
            count=count+1;
        end;
    end;
end;
if test>1
    disp([num2str(count) ' tiles of ' num2str(tx*ty) ' used']);
end;
% pre-whitening correction
r=fftshift(Radius(nu))/nu;
cumsp0=cumsp; %%%
cumsp=cumsp.*exp(r.*(pwExp1-pwExp4*r.^3));

% remove the crystal spots
% cumsp=fftshift(RemoveSpots(fftshift(cumsp),[3.5 3]));

% Filter the square root of the spectrum.
sqrsp=sqrt(cumsp);
if test
    subplot(2,3,4);
    imacs(fftshift(sqrsp.^(disexp*2)));
    LabelImage(nu,df,'A^{-1}');
end;

kernel=fuzzymask(nu,2,w,w/10);
kernel=kernel/sum(kernel(:));
% convolve with kernel
filtsp=real(ifftn(fftn(sqrsp).*fftn(fftshift(kernel))));

diffsp=fftshift(cumsp-filtsp.^2);
% diffsp=fftshift(cumsp./filtsp-filtsp);
% diffsp=fftshift(cumsp);

% look at the subtracted spectrum

% % subplot(2,3,3);  % Compare the smoothed and original spectra
% % radSpecs=[Radial(fftshift(filtsp))' Radial(fftshift(sqrt(cumsp)))'];
% % semilogy(freqs,radSpecs);
% % xlabel('A^{-1}');
% % drawnow;

radsp=Radial(fftshift(cumsp));


% determine resolution limits

rmin=1/(df*minres);
rmax=1/(df*maxres);
outerlimits=fuzzymask(nu,2,rmax,rmax/10);
limits=outerlimits-fuzzymask(nu,2,rmin,1);

if test
    subplot(2,3,5);
    imacs(limits.*diffsp);
    LabelImage(nu,df,'A^{-1}');
    
    subplot(2,3,6);
    rdiff=Radial(diffsp.*limits);
    % % plot(freqs,rdiff);
    drawnow;
end;

Pao=Pa;
Pao.defocus=-Pa.defocus;

halflimits=limits(nu/2+1:nu,:);
halflimits(1,nu/2:nu)=0;  % zero out the upper meridian

halfsp=diffsp(nu/2+1:nu,:).*halflimits;
halfspv=halfsp(:);
halfspv=halfspv/sqrt(halfspv'*halfspv);

% % Here we test for an existing set of references, and compute new ones
% % only if necessary.
% First, check whether the static variable R exists at all.
b=1;
try
    b=numel(R);
catch
    b=0;
end;

% Now see if the problem has changed from previous calls.
ind=0;
for i=1:b
    if StructsEqual(R(i).Pa0,Pa)
        ind=i;
    end;
end;

%if ind==0
    ind=min(b+1,NumRefs);
    if test>1
        disp('Making new CTF references');
    end;
    [R(ind).refs R(ind).refsz]=MakeCTFRefs(nu,res,Pa,halflimits);
    %     [R(ind).refso R(ind).refsz]=MakeCTFRefs(nu,res,Pao,halflimits);
    
    R(ind).refs=reshape(R(ind).refs,nu^2/2,prod(R(ind).refsz));
    %     R(ind).refso=reshape(R(ind).refso,nu^2/2,prod(R(ind).refsz));
    R(ind).Pa0=Pa;
    if test>1
        whos R
    end;
%end;

if test>1
    disp('Cross-correlation')
end;

cc=halfspv'*R(ind).refs;  % Cross correlation done here!

sz=size(cc);
[mxc mxi]=max(cc);
% mxi
[m l i j k]=ind2sub(R(ind).refsz,mxi);

if test>1
    Initial_c=mxc
end;

Ps=Pa;
P=Ps;

P.defocus=Ps.defocus(i);
P.deltadef=Ps.deltadef(j);
P.theta=Ps.theta(k);
P.alpha=Ps.alpha(l);
P.B=P.B(m);
% P

if test
    % Make the split-power spectrum display
    subplot(2,3,5);
    dspect=imscale(limits.*diffsp);
    dref=imscale(limits.*CTF(nu,res,P).^2);
    dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
    imac(dspect);
    LabelImage(nu,df,'A^{-1}');
    drawnow;
end;

% % Make the radial display
% subplot(2,3,6);
% plot(freqs,[Radial(imscale(limits.*diffsp))' Radial(dref)']);
% legend('data','fit');
% xlabel('A^{-1}');

ccx=reshape(cc,R(ind).refsz);

% figure(3); clf;
% colormap jet
% make a contour plot of alpha vs defocus
if test
    subplot(2,3,3);
    if numel(Ps.alpha)>1
        contourf(Ps.alpha,Ps.defocus,squeeze(ccx(m,:,:,j,k))',10);
        xlabel('Alpha');
        ylabel('Defocus');
        colorbar;
        title('Underfocus is positive');
    end;
    subplot(2,3,3);
    contourf(Ps.deltadef,Ps.defocus,squeeze(ccx(m,l,:,:,k)),10);
    xlabel('Delta-defocus');
    ylabel('Defocus');
    colorbar;
    drawnow;
end;
% Do the optimization
% disp('Final optimization');
if test>1
    disp('Final optimization.');
end;
P.theta=0;
if numel(Ps.alpha)<2  % Don't fit alpha
    p=[P.defocus; P.deltadef; P.theta];
    for is=1:2  % do extra restarts of the optimization, to get out local optima.
        %     disp(p');
        p=Simplex('init',p,[0.1 0.1 0.1]);
        for i=1:80
            P.defocus=p(1);
            P.deltadef=p(2);
            P.theta=p(3);
            c=CTFCrossCorr(diffsp,res,limits,P,0);
            p=Simplex(-c);
            %               if mod(i,10)==0
            %                 p'
            %            end;
        end;
        % disp([p(1:2)' 180/pi*p(3:4)']);
    end;
else % We are fitting alpha too
    p=[P.defocus; P.deltadef; P.theta; P.alpha];
    
    for is=1:2  % do extra restarts of the optimization, to get out local optima.
        %     disp(p');
        p=Simplex('init',p,[0.1 0.1 0.1 0.01]);
        for i=1:120
            P.defocus=p(1);
            P.deltadef=p(2);
            P.theta=p(3);
            c=CTFCrossCorr(diffsp,res,limits,P,0);
            p=Simplex(-c);
        end;
    end;
end;

if test
    subplot(2,3,4);
    dspect=imscale(fftshift(sqrsp.^(disexp*2)));
    dspect=(imscale(ImEqualize(GaussFilt(fftshift(sqrsp),0.2)).^2)*1+0);
    [dref mr ma]=imscale(limits.*CTF(nu,res,P).^2);
    dref=ma+mr*CTF(nu,res,P).^2;
    % dref=ImEqualize(CTF(nu,res,P).^2);
    dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
    imac(dspect);
    LabelImage(nu,df,'A^{-1}');
    title(['CC =' num2str(c)]);
    
    subplot(2,3,5);
    dspect=imscale(limits.*diffsp);
    dref=imscale(limits.*CTF(nu,res,P).^2);
    dspect(1:nu/2+1,:)=dref(1:nu/2+1,:);
    imac(dspect);
    LabelImage(nu,df,'A^{-1}');
    title(['\Deltaz=' num2str(-p(1)) '   Astig=' num2str(p(2))...
        '   \theta=' num2str(360/pi*p(3))]);
    
    % Make the radial display
    subplot(2,3,6);
    rs=Radial(limits.*diffsp);
    rs=150*(rs-min(rs))/(max(rs)-min(rs));
    radCorrSpecs=[rs Radial(dref)];
    plot(freqs,radCorrSpecs);
    xlabel('A^{-1}');
    drawnow;
end;
% P.defocus=-P.defocus;  % change the polarity to match scope.
if test>1
    c
    Final_values=[P.defocus P.deltadef 180/pi*P.theta P.alpha]
end;

