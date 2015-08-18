function fft_mu = get_ffts(mu)
fft_mu=zeros(size(mu));
for k=1:size(mu,3)
    fft_mu(:,:,k)=fft2(squeeze(mu(:,:,k)));
end
end