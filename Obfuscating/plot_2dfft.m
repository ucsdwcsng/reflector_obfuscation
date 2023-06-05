function [DP, aoa_pred_d, tof_pred_d] = plot_2dfft(freq, temps_gt, theta_vals, d_vals, ant_sep, gt_aoa_d, gt_aoa_r, gt_d_val_d, gt_d_val_r, delay)
assert(length(size(temps_gt)) == 2, "h must be 2d matrix");
freq_cent = freq(1); % assumes freq is fftshifted, with dc subcarrier at index 1
const = 1j*2*pi/(3e8);
const2 = 1j*2*pi*ant_sep*freq_cent/(3e8);
temps_gt = temps_gt.';
d_rep = const*(freq'.*repmat(d_vals,length(freq),1));
temp = temps_gt*exp(d_rep);
theta_rep = const2*((1:size(temps_gt,1)).*repmat(sin(theta_vals'),1,size(temps_gt,1)));
DP = exp(theta_rep)*(temp);
DP = abs(DP);
[aoa_pred_d,tof_pred_d]=get_aoa_for_least_tof(DP,d_vals,theta_vals);

% figure(1)
% plot_transform_profile(DP, theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], [])
% title(sprintf('without mirage-2dfft: direct-link AoA is %0.1f deg, ToF is %0.1f m; reflected-link AoA is %0.1f deg, ToF is %0.1f m; delayed %0.1f m', rad2deg(gt_aoa_d), gt_d_val_d, rad2deg(gt_aoa_r), gt_d_val_r, delay))
end