function DP = compute_xy_profile2d_generalized(h,d1,d2,AP,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% h : n_subxn_ant complex channel matrix
% d1: search space for X-axis
% d2: search space for Y-axis
% AP: n_antx2 XY locations of the antennas of the AP
% opt: options containting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic implementation
% DP=zeros(length(d1),length(d2));
% freqcomp = (1j*2*pi/(3e8)).*opt.freq;
% for i = 1:length(d1)
%     for j = 1:length(d2)
%         temp = vecnorm(AP-[d1(i),d2(j)],2,2);
%         coeffecients = transpose(h).*exp(temp*freqcomp);
%         DP(i,j)=sum(coeffecients(:));
%     end
% end
% DP = abs(DP)./max(abs(DP(:)));

%% Vectorization of things
DP=zeros(length(d1),length(d2));
[D1,D2] = meshgrid(d1,d2);
AP_mat = repmat(AP,1,1,length(d2));
freqcomp = transpose((1j*2*pi/(3e8)).*opt.freq);
for i = 1:length(d1)
    d1_vec(1,:,:) = [D1(:,i).';D2(:,i).'];
    temp_mat(1,:,:) = squeeze(vecnorm(AP_mat-d1_vec,2,2));
    temp_vec = reshape(temp_mat,1,[]);
    exponent = freqcomp*temp_vec;
    exponent_mat = reshape(exponent,length(freqcomp),size(AP,1),[]);
    coeffecients = h.*exp(exponent_mat);
    sum_coeffecients = reshape(coeffecients,[],length(d2));
    DP(i,:)=sum(sum_coeffecients,1);
end
DP = abs(DP)./max(abs(DP(:)));

end