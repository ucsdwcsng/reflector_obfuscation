function [ch_user, ch_ap, gt_aoa_d, gt_aoa_r, gt_aod_d, gt_aod_r, gt_d_val_d, gt_d_val_r, delay_d_r, delay] = channel_generation(freq, theta_vals, d_vals, ant_sep)
% generate channel from configured environment
[ch_user, ch_ap] = ch_gen(); % channels from ch_gen is tx by rx
S = get_2dsteering_matrix(theta_vals,d_vals);

% obtain channel at user side to compute aod 
% features_user(1,1,:,:) = compute_spotfi_profile_vectorized_two_paths(fftshift(squeeze(ch_user(1,1,:,1,:)), 1),theta_vals,d_vals,opt,S); %
[features_user(1,1,:,:), ~, ~] = plot_2dfft(freq, squeeze(ch_user(1,1,:,1,:)), theta_vals, d_vals, ant_sep, [], [], [], [], []);

[gt_aods, dont_know_tofs] = get_aoa_tof_Twopaths(squeeze(features_user(1,1,:,:)),d_vals,theta_vals); % get aod at user for two paths
gt_aod_d = gt_aods(1,1)*180/pi;
gt_aod_r = gt_aods(1,2)*180/pi;


% obtain channel at ap side to compute aoa
% features_ap(1,1,:,:) = compute_spotfi_profile_vectorized_two_paths(fftshift(squeeze(ch_ap(1,1,:,1,:)), 1),theta_vals,d_vals,opt,S); %
[features_ap(1,1,:,:), ~, ~] = plot_2dfft(freq, squeeze(ch_ap(1,1,:,1,:)), theta_vals, d_vals, ant_sep, [], [], [], [], []);
[gt_aoas, gt_tofs] = get_aoa_tof_Twopaths(squeeze(features_ap(1,1,:,:)),d_vals,theta_vals); % get aod at user for two paths
gt_aoa_d = gt_aoas(1,1)*180/pi;
gt_aoa_r = gt_aoas(1,2)*180/pi;
gt_d_val_d = gt_tofs(1,1); % direct path tof
gt_d_val_r = gt_tofs(1,2); % reflected path tof
delay_d_r = gt_d_val_r-gt_d_val_d;% distfeatures_apance between the direct path and reflected path
delay = 100 ; %15.5; % how much delay you want to add

debug=false;
if debug
 figure(2)
 n_sub=48;
 subplot(2,2,1)
 temps = squeeze(ch_ap(1,1,:,1,:));
 [P, aoa_pred_d, tof_pred_d]=plot_spotfi(temps, theta_vals, d_vals, n_sub);
 plot_transform_profile((P), theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], []);
% plot_transform_profile(squeeze(features_user(1,1,:,:)), theta_vals, d_vals, [], [], [], [], [], []);
title('spotfi: channel at AP')

subplot(2,2,2)
temps_gt = squeeze(ch_ap(1,1,:,1,:));%
[DP1, aoa_pred_d, tof_pred_d]= plot_2dfft(freq, temps_gt, theta_vals, d_vals, ant_sep, gt_aoa_d, gt_aoa_r, gt_d_val_d, gt_d_val_r, delay);
plot_transform_profile(DP1, theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], [])
title('2dfft: channel at AP')

end
end
