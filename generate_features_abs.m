function features = generate_features_abs(channels,ap,theta_vals,d_vals,d1,d2,opt)

n_ap=length(ap);
features = zeros(n_ap,length(d2),length(d1));

for j=1:n_ap
    P = compute_multipath_profile2d_fast_edit(squeeze(channels(:,j,:)),theta_vals,d_vals,opt);
    P_out = convert_spotfi_to_2d(P,theta_vals,d_vals,d2,d1,ap{j});
    features(j,:,:) = (abs(P_out)./max(abs(P_out(:)))).';
end
end