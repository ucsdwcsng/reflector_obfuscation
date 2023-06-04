function [bf_d, bf_r, bf_d_pilots, bf_r_pilots]=bf_v(freq, freq_pilots, gt_aod_d, gt_aod_r, ant_sep)
n_ant=4;
n_sub=48;
bf_d = ones(48, 4);
bf_r = ones(48, 4);

bf_d_pilots = ones(4,4);
bf_r_pilots = ones(4,4);

bf_t_d = ones(48, 4);
bf_t_r = ones(48, 4);
for ant_dix = 1:n_ant
  % for data sc
  for subcarrier_idx =1:n_sub
      f = freq(subcarrier_idx);
      bf_d(subcarrier_idx,ant_dix) = exp(1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aod_d*pi/180))/3e8)/2; 
      bf_r(subcarrier_idx,ant_dix) = exp(1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aod_r*pi/180))/3e8)/2;
  end
  % for pilot sc
  for subcarrier_idx = 1:4
      f = freq_pilots(subcarrier_idx);
      bf_d_pilots(subcarrier_idx,ant_dix) = exp(1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aod_d*pi/180))/3e8)/2;
      bf_r_pilots(subcarrier_idx,ant_dix) = exp(1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aod_r*pi/180))/3e8)/2;
  end
end
end