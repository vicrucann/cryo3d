function norms = get_norms(I_fft_ctf, mu_fft, map)
N=size(I,3);
norms=zeros(1,N);
mapped_k = I_fft_ctf(map);
mu_fft_mapped = mu_fft(:,:,mapped_k);
diff = I_fft_ctf - mu_fft_mapped;
diff_resh = reshape(diff, [size(diff,1)*size(diff,2), size(diff,3)]);
for k=1:N
    norms(k) = norm(diff_resh(:,k));
end
%for k=1:N
%    j = map(k);
%    t1=squeeze(I_fft_ctf(:,:,k));
%    t2=squeeze(mu_fft(:,:,j));
%    norms(k) = norm2(t1(:) - t2(:));
%end
end