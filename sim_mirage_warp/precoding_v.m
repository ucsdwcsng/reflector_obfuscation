function [apply_ch, ch_all, bf_null_delay] =precoding_v(freq, freq_pilots, ch_ap, SC_IND_DATA, SC_IND_PILOTS, ant_sep, gt_aod_d, gt_aod_r, gt_aoa_d, gt_aoa_r, gt_d_val_d, gt_d_val_r)
n_sub=48;

% for data sc
bf_null_delay = ones(48, 4, 4);
nulling = ones(48, 4, 4);
for subcarrier_idx =1:n_sub
    f = freq(subcarrier_idx);
    % beamforming&nulling&delaying weight
%     bf_null_delay(subcarrier_idx, :,:) = (null_nd(subcarrier_idx, :).')*bf_d(subcarrier_idx, :).*exp(-1j*2*pi*f*delay/3e8) ...
%                                         + (null_d(subcarrier_idx,:).')*bf_r(subcarrier_idx,:);
%      bf_null_delay(subcarrier_idx, :,:) = (null_nd(subcarrier_idx, :).')*bf_d(subcarrier_idx, :).*exp(-1j*2*pi*f*delay/3e8);
    % beamforming&nulling weight
%     bf_null_delay(subcarrier_idx, :,:) = (null_d(subcarrier_idx, :).')*transpose(bf_r(subcarrier_idx, :).');
    bf_null_delay(subcarrier_idx, :,:) = eye(4);

    %beamforming and delaying wight
%     bf_null_delay(subcarrier_idx,:,:) = conj(bf_d(subcarrier_idx, :).')*bf_d(subcarrier_idx, :).*exp(-1j*2*pi*f*delay/3e8) ...
%                                         + conj(bf_r(subcarrier_idx,:).')*bf_r(subcarrier_idx,:);
%     bf_null_delay(subcarrier_idx,:,:) = conj(bf_r(subcarrier_idx,:).')*bf_r(subcarrier_idx,:);

    %nulling weight weight
%     bf_null_delay(subcarrier_idx, :,:) = conj(null_d_space(:,1))*transpose(null_d_space(:,1));
end

% for pilot sc
bf_null_delay_pilots = ones(4, 4, 4);
nulling_pilots = ones(4, 4, 4);
for subcarrier_idx =1:4
    f = freq_pilots(subcarrier_idx);
%     bf_null_delay_pilots(subcarrier_idx, :,:) = (null_nd_pilots(subcarrier_idx, :).')*bf_d_pilots(subcarrier_idx, :).*exp(-1j*2*pi*f*delay/3e8) + ...
%                                             (null_d_pilots(subcarrier_idx,:).')*bf_r_pilots(subcarrier_idx,:);
%     bf_null_delay_pilots(subcarrier_idx, :,:) = (null_nd_pilots(subcarrier_idx, :).')*bf_d_pilots(subcarrier_idx, :).*exp(-1j*2*pi*f*delay/3e8);

%   bf_null_delay_pilots(subcarrier_idx, :,:) = (null_d_pilots(subcarrier_idx, :).')*transpose(bf_r(subcarrier_idx, :).');
    bf_null_delay_pilots(subcarrier_idx, :,:) = eye(4);
%   bf_null_delay_pilots(subcarrier_idx,:,:) = conj(bf_d_pilots(subcarrier_idx, :).')*bf_d_pilots(subcarrier_idx, :).*exp(-1j*2*pi*f*delay/3e8) + ...
%                                               conj(bf_r_pilots(subcarrier_idx,:).')*bf_r_pilots(subcarrier_idx,:);
%   bf_null_delay_pilots(subcarrier_idx,:,:) = conj(bf_r_pilots(subcarrier_idx,:).')*bf_r_pilots(subcarrier_idx,:);
%   bf_null_delay_pilots(subcarrier_idx, :,:) = (null_d_space_pilots(:,1))*transpose(null_d_space_pilots(:,1));
end

%% Combine all the weights together 
ch_all = ones(64, 4, 4);
ch_all(SC_IND_DATA, :, :) = bf_null_delay;  %48, 4, 4
ch_all(SC_IND_PILOTS,:, :) = bf_null_delay_pilots;



%% define wireless channel using aod at tx, aoa at rx 
% h = zeros(48,4,4);
apply_ch  = zeros(64,4,4);

% for pilot sc
h2 = zeros(4, 4, 4);
for subcarrier_idx =1:4
    for rx=1:4 %rx
        for tx=1:4 %tx

            f = freq_pilots(subcarrier_idx);

            t_t_d = ((tx-1) * ant_sep * sin(gt_aod_d*pi/180))/3e8;
            t_t_r = ((tx-1) * ant_sep * sin(gt_aod_r*pi/180))/3e8;

            t_r_d = ((rx-1) * ant_sep * sin(gt_aoa_d*pi/180))/3e8;
            t_r_r = ((rx-1) * ant_sep * sin(gt_aoa_r*pi/180))/3e8;

            h_d = exp(-1j*2*pi*f*(t_t_d+t_r_d+gt_d_val_d/3e8));

            h_r = exp(-1j*2*pi*f*(t_t_r+t_r_r+gt_d_val_r/3e8));

            h2(subcarrier_idx, rx, tx) = h_d + h_r;
        end
    end
end

% apply_ch(SC_IND_DATA, :,:) = h;

apply_ch(SC_IND_DATA, :,:) = permute(squeeze(ch_ap(1,1,:,:,:)), [1, 3, 2]);

apply_ch(SC_IND_PILOTS, :,:) = h2;


end