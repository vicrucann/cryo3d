function I_fft_ctf = ctf_filtered_img(I, defocus, pixA, lambda, B, ampContrast)

if (nargin<6)
    ampContrast=0.1;
end
if (nargin<5)
    B=30;
end

N=size(I,3);

% get FFTs of I
I_fft=get_ffts(I);

% get CTFs
M=size(defocus, 2);
ctfs=zeros(1,M);
for j=1:M
    ctfs(j)=CTF(size(I,1),pixA,lambda,defocus(j)/1e4,0,B,ampContrast);
end

assert(N<=M, 'size(I,3) must be less or eq to size(defocus,2)');

% calculate I_fft * CTF
I_fft_ctf = zeros(size(I_fft));
for k=1:N
    I_fft_ctf(:,:,k) = I_fft(:,:,k) * ctfs(k);
end
end