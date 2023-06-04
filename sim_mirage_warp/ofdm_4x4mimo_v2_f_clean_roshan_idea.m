clear
clc
close all

% Params:
NUM_PACKETS             = 10;
USE_WARPLAB_TXRX        = 0;           % Enable WARPLab-in-the-loop (otherwise sim-only)
WRITE_PNG_FILES         = 0;           % Enable writing plots to PNG
CHANNEL                 = 11;          % Channel to tune Tx and Rx radios

% Waveform params
N_OFDM_SYMS             = 4000;        % Number of OFDM symbols (must be even valued)
MOD_ORDER               = 2;         % Modulation order (2/4/16/64/128 = BSPK/QPSK/16-QAM/64-QAM/128-QAM)
TX_SCALE                = 1;         % Scale for Tx waveform ([0:1])
INTERP_RATE             = 2;           % Interpolation rate (must be 2)
TX_SPATIAL_STREAM_SHIFT = 3;           % Number of samples to shift the transmission from RFB

% OFDM params
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)

% Rx processing params
FFT_OFFSET                    = 4;           % Number of CP samples to use in FFT (on average)
LTS_CORR_THRESH               = 0.8;         % Normalized threshold for LTS correlation
DO_APPLY_CFO_CORRECTION       = 1;           % Enable CFO estimation/correction
DO_APPLY_PHASE_ERR_CORRECTION = 1;           % Enable Residual CFO estimation/correction0
DO_APPLY_SFO_CORRECTION       = 1;           % Enable SFO estimation/correction
DECIMATE_RATE                 = INTERP_RATE;

% WARPLab experiment params
USE_AGC                 = true;        % Use the AGC if running on WARP hardware
MAX_TX_LEN              = 2^19;        % Maximum number of samples to use for this experiment
SAMP_PADDING            = 100;         % Extra samples to receive to ensure both start and end of waveform visible


% Use sane defaults for hardware-dependent params in sim-only version
maximum_buffer_len  = min(MAX_TX_LEN, 2^20);
SAMP_FREQ           = 500e3;
example_mode_string = 'sim';

%% Define a half-band 2x interpolation filter response
interp_filt2 = zeros(1,43);
interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
interp_filt2(22) = 16384;
interp_filt2 = interp_filt2./max(abs(interp_filt2));

% Define the preamble
% Note: The STS symbols in the preamble meet the requirements needed by the
% AGC core at the receiver. Details on the operation of the AGC are
% available on the wiki: http://warpproject.org/trac/wiki/WARPLab/AGC
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];


sts_t = ifft(sqrt(13/6).*sts_f, 64);
sts_t = sts_t(1:16);


% LTS for CFO and channel estimation
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];

sts_t_rep = repmat(sts_t, 1, 30);

%% parameters for array and mirage
n_ant = 4;%number of antennas
n_sub = length(SC_IND_DATA);%48 data subcarriers
c = 3e8;%speed of light
ant_sep = 0.0259;% antenna separation
theta_vals = -pi/2:0.01:pi/2; % AoA search space
d_vals = -30:0.1:150; %tof search space in distance (meter)
BW = 20e6;%bandwidth
subcarrier_indices = SC_IND_DATA-32; %symmetric data subcarriers 
freq = fftshift(5.18e9 + subcarrier_indices.*BW./64);%frequency of each data subcarrier
subcarrier_indices_pilots = SC_IND_PILOTS - 32; % symmetric pilot subcarriers
freq_pilots = fftshift(5.18e9 + subcarrier_indices_pilots.*BW./64);%frequency of each pilot subcarrier
gt_d_val_d = 10;% ground-truth direct path distance
delay_d_r = 5;
gt_d_val_r = gt_d_val_d+delay_d_r;%  ground truth reflected path distance
delay = 25;% delay for direct path which should be larger than delay_d_r
gt_aoa_d = pi/6; % aoa for direct path
gt_aoa_r = -pi/4;% aoa for reflected path
gt_aod_d = gt_aoa_d;% aod for direct path
gt_aod_r = -gt_aoa_r;% aod for reflected path


%% define beamforming vector based on aod
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
%       bf_d(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aod_d))/3e8)/2; 
%       bf_r(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aod_r))/3e8)/2;
      %roshan's aoa 
      bf_d(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aoa_d))/3e8)/2; 
      bf_r(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aoa_r))/3e8)/2;

  end
  % for pilot sc
  for subcarrier_idx = 1:4
      f = freq_pilots(subcarrier_idx);
%       bf_d_pilots(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aod_d))/3e8)/2;
%       bf_r_pilots(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aod_r))/3e8)/2;
       %roshan's aoa 
      bf_d_pilots(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aoa_d))/3e8)/2; 
      bf_r_pilots(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aoa_r))/3e8)/2;
  end
end


%% beamforming and nulling direction
% For data sc
null_d = zeros(48, 4);
null_nd = zeros(48, 4);
for i=1:48
 null_d_space = null(repmat((bf_d(i,:)), [4,1]));
 temp_w_1 = conj(bf_r(i,:))*null_d_space(:,1);
 temp_w_2 = conj(bf_r(i,:))*null_d_space(:,2);
 temp_w_3 = conj(bf_r(i,:))*null_d_space(:,3);
 temp_w = temp_w_1*null_d_space(:,1)+temp_w_2*null_d_space(:,2)+temp_w_3*null_d_space(:,3);
 null_d(i,:) = temp_w/vecnorm(temp_w, 2,1);


 null_nd_space = null(repmat((bf_r(i,:)), [4, 1]));  % conjugate nullspace   
 temp_w_1 = conj(bf_d(i,:))*null_nd_space(:,1);
 temp_w_2 = conj(bf_d(i,:))*null_nd_space(:,2);
 temp_w_3 = conj(bf_d(i,:))*null_nd_space(:,3);
 temp_w = temp_w_1*null_nd_space(:,1)+temp_w_2*null_nd_space(:,2)+temp_w_3*null_nd_space(:,3);
 null_nd(i,:) = temp_w/vecnorm(temp_w, 2,1);
end

% For pilots sc
null_d_pilots = zeros(4, 4);
null_nd_pilots = zeros(4, 4);
for i=1:4
 null_d_space_pilots = null(repmat((bf_d_pilots(i,:)), [4,1]));
 temp_w_1 = conj(bf_r_pilots(i,:))*null_d_space_pilots(:,1);
 temp_w_2 = conj(bf_r_pilots(i,:))*null_d_space_pilots(:,2);
 temp_w_3 = conj(bf_r_pilots(i,:))*null_d_space_pilots(:,3);
 temp_w = temp_w_1*null_d_space_pilots(:,1)+temp_w_2*null_d_space_pilots(:,2)+temp_w_3*null_d_space_pilots(:,3);
 null_d_pilots(i,:) = temp_w/vecnorm(temp_w, 2,1);


 null_nd_space_pilots = null(repmat((bf_r_pilots(i,:)), [4, 1]));  % conjugate nullspace   
 temp_w_1 = conj(bf_d_pilots(i,:))*null_nd_space_pilots(:,1);
 temp_w_2 = conj(bf_d_pilots(i,:))*null_nd_space_pilots(:,2);
 temp_w_3 = conj(bf_d_pilots(i,:))*null_nd_space_pilots(:,3);
 temp_w = temp_w_1*null_nd_space_pilots(:,1)+temp_w_2*null_nd_space_pilots(:,2)+temp_w_3*null_nd_space_pilots(:,3);
 null_nd_pilots(i,:) = temp_w/vecnorm(temp_w, 2,1);
end


lts_t = ifft(lts_f, 64);


%% define wireless channel using aod at tx, aoa at rx 
% for data sc
h = zeros(48,4,4);
apply_ch  = zeros(64,4,4);
for subcarrier_idx =1:48
    for rx=1:4 %rx
        for tx=1:4 %tx
            f = freq(subcarrier_idx);

            t_t_d = ((tx-1) * ant_sep * sin(gt_aod_d))/3e8;
            t_t_r = ((tx-1) * ant_sep * sin(gt_aod_r))/3e8;

            t_r_d = ((rx-1) * ant_sep * sin(gt_aoa_d))/3e8;
            t_r_r = ((rx-1) * ant_sep * sin(gt_aoa_r))/3e8;

            h_d = exp(-1j*2*pi*f*(t_t_d+t_r_d+gt_d_val_d/3e8));

            h_r = exp(-1j*2*pi*f*(t_t_r+t_r_r+gt_d_val_r/3e8));

            h(subcarrier_idx, rx, tx) = h_d + 1/2*h_r;
        end
    end
end
% for pilot sc
h2 = zeros(4, 4, 4);
for subcarrier_idx =1:4
    for rx=1:4 %rx
        for tx=1:4 %tx

            f = freq_pilots(subcarrier_idx);

            t_t_d = ((tx-1) * ant_sep * sin(gt_aod_d))/3e8;
            t_t_r = ((tx-1) * ant_sep * sin(gt_aod_r))/3e8;

            t_r_d = ((rx-1) * ant_sep * sin(gt_aoa_d))/3e8;
            t_r_r = ((rx-1) * ant_sep * sin(gt_aoa_r))/3e8;

            h_d = exp(-1j*2*pi*f*(t_t_d+t_r_d+gt_d_val_d/3e8));

            h_r = exp(-1j*2*pi*f*(t_t_r+t_r_r+gt_d_val_r/3e8));

            h2(subcarrier_idx, rx, tx) = h_d + h_r;
        end
    end
end

apply_ch(SC_IND_DATA, :,:) = h;

apply_ch(SC_IND_PILOTS, :,:) = h2;





% for data sc
bf_null_delay = ones(48, 4, 4);
nulling = ones(48, 4, 4);
for subcarrier_idx =1:n_sub
    f = freq(subcarrier_idx);
    % beamforming&nulling&delaying weight
%     bf_null_delay(subcarrier_idx, :,:) = (null_nd(subcarrier_idx, :).')*bf_d(subcarrier_idx, :).*exp(-1j*2*pi*f*delay/3e8) + (null_d(subcarrier_idx,:).')*bf_r(subcarrier_idx,:);
    % beamforming&nulling weight
%     bf_null_delay(subcarrier_idx, :,:) = (null_d(subcarrier_idx, :).')*transpose(bf_r(subcarrier_idx, :).');
%     bf_null_delay(subcarrier_idx, :,:) = eye(4);

    %beamforming and delaying wight
    %bf_null_delay(subcarrier_idx,:,:) = conj(bf_d(subcarrier_idx, :).')*bf_d(subcarrier_idx, :).*exp(1j*2*pi*f*delay/3e8)+ conj(bf_r(subcarrier_idx,:).')*bf_r(subcarrier_idx,:);
    %nulling weight weight
%     nulling(subcarrier_idx, :,:) = conj(null_d_space(:,1))*transpose(null_d_space(:,1));
    %roshan's idea
    bf_null_delay(subcarrier_idx,:,:)=conj(h(subcarrier_idx,:,:))*((transpose(bf_d(subcarrier_idx,:))*conj(bf_d(subcarrier_idx,:)))*exp(-1j*2*pi*f*delay/3e8)+transpose(bf_r(subcarrier_idx,:))*conj(bf_r(subcarrier_idx,:)));
end

% for pilot sc
bf_null_delay_pilots = ones(4, 4, 4);
nulling_pilots = ones(4, 4, 4);
for subcarrier_idx =1:4
    f = freq_pilots(subcarrier_idx);
%     bf_null_delay_pilots(subcarrier_idx, :,:) = (null_nd_pilots(subcarrier_idx, :).')*bf_d_pilots(subcarrier_idx, :).*exp(-1j*2*pi*f*delay/3e8) + ...
%                                             (null_d_pilots(subcarrier_idx,:).')*bf_r_pilots(subcarrier_idx,:);
%   bf_null_delay_pilots(subcarrier_idx, :,:) = (null_d_pilots(subcarrier_idx, :).')*transpose(bf_r(subcarrier_idx, :).');
%     bf_null_delay_pilots(subcarrier_idx, :,:) = eye(4);
%   bf_null_delay_pilots(subcarrier_idx,:,:) = conj(bf_d_pilots(subcarrier_idx, :).')*bf_d_pilots(subcarrier_idx, :).*exp(1j*2*pi*f*delay/3e8)+ ...
%                                               conj(bf_r_pilots(subcarrier_idx,:).')*bf_r_pilots(subcarrier_idx,:);
%   nulling_pilots(subcarrier_idx, :,:) = (null_d_space_pilots(:,1))*transpose(null_d_space_pilots(:,1));
   % roshan's idea
  bf_null_delay_pilots(subcarrier_idx,:,:)=conj(h(subcarrier_idx,:,:))*((transpose(bf_d(subcarrier_idx,:))*conj(bf_d(subcarrier_idx,:)))*exp(-1j*2*pi*f*delay/3e8)+transpose(bf_r(subcarrier_idx,:))*conj(bf_r(subcarrier_idx,:)));

end

% Combine all the weights together 
ch_all = ones(64, 4, 4);
ch_all(SC_IND_DATA, :, :) = bf_null_delay;  %48, 4, 4
ch_all(SC_IND_PILOTS,:, :) = bf_null_delay_pilots;








%% debug the aoa-tof profile using the precoding weight we generate
debug = true;
if debug
    ch_f = zeros(48, 4,4);
    for sub =1:48
        ch_f(sub, :,:)  = squeeze(h(sub,:,:))*squeeze(bf_null_delay(sub,:,:));
    end


    figure(100)
    temps = squeeze(ch_f(:,:,1));
    [P, aoa_pred_d, tof_pred_d]=plot_spotfi(fftshift(temps, 1), theta_vals, d_vals, n_sub);

    plot_transform_profile(db(P), theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], []);
    title('spotfi: ground-truth channel multiplying precoding weight')

    figure(101)
    [DP1, aoa_pred_d, tof_pred_d]= plot_2dfft(freq, temps, theta_vals, d_vals, ant_sep, gt_aoa_d, gt_aoa_r, gt_d_val_d, gt_d_val_r, delay);
    plot_transform_profile(DP1, theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], [])
    title('2dfft: ground-truth channel multiplying precoding weight')


end


%% applying the precoding weight (BUG IS HERE: HOW WILL WE APPLY THE PRECODING WEIGHT ON MIMO PREAMBLE)
sts_f_across_antenna = zeros(64, 4,1);
lts_f_across_antenna = zeros(64, 4,4);

for sub=1:64
    sts_f_across_antenna(sub,:,:) = transpose(squeeze(ch_all(sub, :,:)))*repmat(sts_f(1, sub), 4, 1);
end

for rx = 1:4
    for tx = 1:4
        lts_f_across_antenna(:,rx,tx) = squeeze(ch_all(:, rx,tx)).*lts_f.';
    end
end

sts_pre_A = sts_f_across_antenna(:,1,:).';% 64, 1
sts_pre_B = sts_f_across_antenna(:,2,:).';
sts_pre_C = sts_f_across_antenna(:,3,:).';
sts_pre_D = sts_f_across_antenna(:,4,:).';



sts_t_w_A_temp = ifft(sqrt(13/6).*sts_pre_A, 64); %%% checking the sts_f to be 1, 4, 64, ch_all_A is 64, 4
sts_t_w_A = sts_t_w_A_temp(1:16);
sts_t_rep_A = repmat(sts_t_w_A, 1, 30);


lts_t_AA = ifft(lts_f_across_antenna(:, 1, 1), 64).';
lts_t_AB = ifft(lts_f_across_antenna(:, 2, 1), 64).';
lts_t_AC = ifft(lts_f_across_antenna(:, 3, 1), 64).';
lts_t_AD = ifft(lts_f_across_antenna(:, 4, 1), 64).';

preamble_legacy_A = [sts_t_rep_A, lts_t_AA(33:64), lts_t_AA, lts_t_AA];
preamble_mimo_A = [lts_t_AA(33:64), lts_t_AA, ...
                   lts_t_AB(33:64), lts_t_AB, ... 
                   lts_t_AC(33:64), lts_t_AC, ...
                   lts_t_AD(33:64), lts_t_AD];


%apply weight on frequency domain for preamble on antenna B
sts_t_w_B_temp = ifft(sqrt(13/6).*sts_pre_B, 64);
sts_t_w_B = sts_t_w_B_temp(1:16);
sts_t_rep_B = repmat(sts_t_w_B, 1, 30);


lts_t_BA = ifft(lts_f_across_antenna(:, 1, 2), 64).';
lts_t_BB = ifft(lts_f_across_antenna(:, 2, 2), 64).';
lts_t_BC = ifft(lts_f_across_antenna(:, 3, 2), 64).';
lts_t_BD = ifft(lts_f_across_antenna(:, 4, 2), 64).';

preamble_legacy_B = [circshift(sts_t_rep_B, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];
preamble_mimo_B = [lts_t_BA(33:64), lts_t_BA, ...
                   lts_t_BB(33:64), lts_t_BB, ... 
                   lts_t_BC(33:64), lts_t_BC, ...
                   lts_t_BD(33:64), lts_t_BD];

%apply weight on frequency domain for preamble on antenna C
sts_t_w_C_temp = ifft(sqrt(13/6).*sts_pre_C, 64);
sts_t_w_C = sts_t_w_C_temp(1:16);
sts_t_rep_C = repmat(sts_t_w_C, 1, 30);
lts_t_CA = ifft(lts_f_across_antenna(:, 1, 3), 64).';
lts_t_CB = ifft(lts_f_across_antenna(:, 2, 3), 64).';
lts_t_CC = ifft(lts_f_across_antenna(:, 3, 3), 64).';
lts_t_CD = ifft(lts_f_across_antenna(:, 4, 3), 64).';

preamble_legacy_C = [circshift(sts_t_rep_C, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];
preamble_mimo_C = [lts_t_CA(33:64), lts_t_CA, ...
                   lts_t_CB(33:64), lts_t_CB, ... 
                   lts_t_CC(33:64), lts_t_CC, ...
                   lts_t_CD(33:64), lts_t_CD];

%apply weight on frequency domain for preamble on antenna D
sts_t_w_D_temp = ifft(sqrt(13/6).*sts_pre_D, 64);
sts_t_w_D = sts_t_w_D_temp(1:16);
sts_t_rep_D = repmat(sts_t_w_D, 1, 30);
lts_t_DA = ifft(lts_f_across_antenna(:, 1, 4), 64).';
lts_t_DB = ifft(lts_f_across_antenna(:, 2, 4), 64).';
lts_t_DC = ifft(lts_f_across_antenna(:, 3, 4), 64).';
lts_t_DD = ifft(lts_f_across_antenna(:, 4, 4), 64).';

preamble_legacy_D = [circshift(sts_t_rep_D, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];
preamble_mimo_D = [lts_t_DA(33:64), lts_t_DA, ...
                   lts_t_DB(33:64), lts_t_DB, ... 
                   lts_t_DC(33:64), lts_t_DC, ...
                   lts_t_DD(33:64), lts_t_DD];



% We break the construction of our preamble into two pieces. First, the
% legacy portion, is used for CFO recovery and timing synchronization at
% the receiver. The processing of this portion of the preamble is SISO.
% Second, we include explicit MIMO channel training symbols.

% Legacy Preamble

% Use 30 copies of the 16-sample STS for extra AGC settling margin
% To avoid accidentally beamforming the preamble transmissions, we will
% let RFA be dominant and handle the STS and first set of LTS. We will
% append an extra LTS sequence from RFB so that we can build out the
% channel matrix at the receiver


preamble_A = [preamble_legacy_A, preamble_mimo_A];
preamble_B = [preamble_legacy_B, preamble_mimo_B];
preamble_C = [preamble_legacy_C, preamble_mimo_C];
preamble_D = [preamble_legacy_D, preamble_mimo_D];

% Sanity check variables that affect the number of Tx samples
if(SAMP_PADDING + INTERP_RATE*((N_OFDM_SYMS/4 * (N_SC + CP_LEN)) + length(preamble_A) + 100) > maximum_buffer_len)
    fprintf('Too many OFDM symbols for TX_NUM_SAMPS!\n');
    fprintf('Raise TX_NUM_SAMPS to %d, or \n', SAMP_PADDING + INTERP_RATE*((N_OFDM_SYMS/4 * (N_SC + CP_LEN)) + length(preamble_A) + 100));
    fprintf('Reduce N_OFDM_SYMS to %d\n',  2*(floor(( (maximum_buffer_len/INTERP_RATE)-length(preamble_A)-100-SAMP_PADDING )/( N_SC + CP_LEN )) - 1));
    return;
end


ch_sc = zeros(4, 4, N_SC);

%% Generate a payload of random integers
 tx_data = randi(MOD_ORDER, 1, N_DATA_SYMS) - 1;


% Functions for data -> complex symbol mapping (like qammod, avoids comm toolbox requirement)
% These anonymous functions implement the modulation mapping from IEEE 802.11-2012 Section 18.3.5.8
modvec_bpsk   =  (1/sqrt(2))  .* [-1 1];
modvec_16qam  = (1/sqrt(10))  .* [-3 -1 +3 +1];
modvec_64qam  =  (1/sqrt(43)) .* [-7 -5 -1 -3 +7 +5 +1 +3];


mod_fcn_bpsk  = @(x) complex(modvec_bpsk(1+x),0);
mod_fcn_qpsk  = @(x) complex(modvec_bpsk(1+bitshift(x, -1)), modvec_bpsk(1+mod(x, 2)));
mod_fcn_16qam = @(x) complex(modvec_16qam(1+bitshift(x, -2)), modvec_16qam(1+mod(x,4)));
mod_fcn_64qam = @(x) complex(modvec_64qam(1+bitshift(x, -3)), modvec_64qam(1+mod(x,8)));


% Map the data values on to complex symbols
switch MOD_ORDER
    case 2         % BPSK
        tx_syms = arrayfun(mod_fcn_bpsk, tx_data);
    case 4         % QPSK
        tx_syms = arrayfun(mod_fcn_qpsk, tx_data);
    case 8
        refconst = qammod(0:7, 8);
        y = qammod(tx_data, 8);
        nf = modnorm(refconst, 'peakpow', 1)*1.34;
        tx_syms = nf*y;         
    case 16        % 16-QAM
        tx_syms = arrayfun(mod_fcn_16qam, tx_data);
    case 64        %64-QAM
          tx_syms = arrayfun(mod_fcn_64qam, tx_data);
    case 128       %128-QAM
        refconst = qammod(0:127, 128);
        y = qammod(tx_data, 128);
        nf = modnorm(refconst,'peakpow',1)*1.34;
        tx_syms = nf*y;
    case 256      %256-QAM
        refconst = qammod(0:255, 256);
        y = qammod(tx_data, 256);
        nf = modnorm(refconst,'peakpow',1)*1.34;
        tx_syms = nf*y;
    case 512      %512-QAM
        refconst = qammod(0:511, 512);
        y = qammod(tx_data, 512);
        nf = modnorm(refconst,'peakpow',1)*1.34;
        tx_syms = nf*y;
    case 1024      %1024-QAM
        refconst = qammod(0:1023, 1024);
        y = qammod(tx_data, 1024);
        nf = modnorm(refconst,'peakpow',1)*1.34;
        tx_syms = nf*y;     
    otherwise
        fprintf('Invalid MOD_ORDER (%d)!  Must be in [2, 4, 16]\n', MOD_ORDER);
        return;
end



% Reshape the symbol vector into two different spatial streams
tx_syms_space_mat = reshape(tx_syms, 4, length(tx_syms)/4);


% Break up the matrix into a vector for each antenna
tx_syms_A = tx_syms_space_mat(1,:);
tx_syms_B = tx_syms_space_mat(2,:);
tx_syms_C = tx_syms_space_mat(3,:);
tx_syms_D = tx_syms_space_mat(4,:);

% Reshape the symbol vector to a matrix with one column per OFDM symbol
tx_syms_mat_A = reshape(tx_syms_A, length(SC_IND_DATA), N_OFDM_SYMS/4);
tx_syms_mat_B = reshape(tx_syms_B, length(SC_IND_DATA), N_OFDM_SYMS/4);
tx_syms_mat_C = reshape(tx_syms_C, length(SC_IND_DATA), N_OFDM_SYMS/4);
tx_syms_mat_D = reshape(tx_syms_D, length(SC_IND_DATA), N_OFDM_SYMS/4);

% Define the pilot tone values as BPSK symbols
%  We will transmit pilots only on RF A
pilots_A = [1 1 -1 1].';
pilots_B = [0 0 0 0].';
pilots_C = [0 0 0 0].';
pilots_D = [0 0 0 0].';


% Repeat the pilots across all OFDM symbols
pilots_mat_A = repmat(pilots_A, 1, N_OFDM_SYMS/4);
pilots_mat_B = repmat(pilots_B, 1, N_OFDM_SYMS/4);
pilots_mat_C = repmat(pilots_C, 1, N_OFDM_SYMS/4);
pilots_mat_D = repmat(pilots_D, 1, N_OFDM_SYMS/4);


%% IFFT

% Construct the IFFT input matrix
ifft_in_mat_A = zeros(N_SC, N_OFDM_SYMS/4);
ifft_in_mat_B = zeros(N_SC, N_OFDM_SYMS/4);
ifft_in_mat_C = zeros(N_SC, N_OFDM_SYMS/4);
ifft_in_mat_D = zeros(N_SC, N_OFDM_SYMS/4);



% Insert the data and pilot values; other subcarriers will remain at 0lts_t
ifft_in_mat_A(SC_IND_DATA, :)   = tx_syms_mat_A;
ifft_in_mat_A(SC_IND_PILOTS, :) = pilots_mat_A;

ifft_in_mat_B(SC_IND_DATA, :)   = tx_syms_mat_B;
ifft_in_mat_B(SC_IND_PILOTS, :) = pilots_mat_B;


ifft_in_mat_C(SC_IND_DATA, :)   = tx_syms_mat_C;
ifft_in_mat_C(SC_IND_PILOTS, :) = pilots_mat_C;


ifft_in_mat_D(SC_IND_DATA, :)   = tx_syms_mat_D;
ifft_in_mat_D(SC_IND_PILOTS, :) = pilots_mat_D;

%Perform the IFFT and applying precoding weight to data symbols

pre_data = zeros(64, 4, N_OFDM_SYMS/4);
for sub =1:64
      temp_A = ifft_in_mat_A(sub,:);
      temp_B = ifft_in_mat_B(sub,:);
      temp_C = ifft_in_mat_C(sub,:);
      temp_D = ifft_in_mat_D(sub,:);
      temp = [temp_A; temp_B; temp_C; temp_D];
      pre_data(sub,:,:) = squeeze(ch_all(sub, :,:))*temp;%4,1000
end
pre_ifft_in_mat_A = squeeze(pre_data(:,1,:));
pre_ifft_in_mat_B = squeeze(pre_data(:,2,:));
pre_ifft_in_mat_C = squeeze(pre_data(:,3,:));
pre_ifft_in_mat_D = squeeze(pre_data(:,4,:));



tx_payload_mat_A = ifft(pre_ifft_in_mat_A, N_SC, 1);
tx_payload_mat_B = ifft(pre_ifft_in_mat_B, N_SC, 1);
tx_payload_mat_C = ifft(pre_ifft_in_mat_C, N_SC, 1);
tx_payload_mat_D = ifft(pre_ifft_in_mat_D, N_SC, 1);



% Insert the cyclic prefix
if(CP_LEN > 0)
    tx_cp = tx_payload_mat_A((end-CP_LEN+1 : end), :);
    tx_payload_mat_A = [tx_cp; tx_payload_mat_A];

    tx_cp = tx_payload_mat_B((end-CP_LEN+1 : end), :);
    tx_payload_mat_B = [tx_cp; tx_payload_mat_B];
    
    tx_cp = tx_payload_mat_C((end-CP_LEN+1 : end), :);
    tx_payload_mat_C = [tx_cp; tx_payload_mat_C];  

    tx_cp = tx_payload_mat_D((end-CP_LEN+1 : end), :);
    tx_payload_mat_D = [tx_cp; tx_payload_mat_D]; 
end

% Reshape to a vector
tx_payload_vec_A = reshape(tx_payload_mat_A, 1, numel(tx_payload_mat_A));
tx_payload_vec_B = reshape(tx_payload_mat_B, 1, numel(tx_payload_mat_B));
tx_payload_vec_C = reshape(tx_payload_mat_C, 1, numel(tx_payload_mat_C));
tx_payload_vec_D = reshape(tx_payload_mat_D, 1, numel(tx_payload_mat_D));



% Construct the full time-domain OFDM waveform
tx_vec_A = [preamble_A tx_payload_vec_A];
tx_vec_B = [preamble_B tx_payload_vec_B];
tx_vec_C = [preamble_C tx_payload_vec_C];
tx_vec_D = [preamble_D tx_payload_vec_D];


% Pad with zeros for transmission
tx_vec_padded_A = [tx_vec_A zeros(1,50)];
tx_vec_padded_B = [tx_vec_B zeros(1,50)];
tx_vec_padded_C = [tx_vec_C zeros(1,50)];
tx_vec_padded_D = [tx_vec_D zeros(1,50)];


%% Interpolate
if(INTERP_RATE == 1)
    tx_vec_air_A = tx_vec_padded_A;
    tx_vec_air_B = tx_vec_padded_B;
    tx_vec_air_C = tx_vec_padded_C;
    tx_vec_air_D = tx_vec_padded_D;
elseif(INTERP_RATE == 2)
	% Zero pad then filter (same as interp or upfirdn without signal processing toolbox)
    tx_vec_2x_A = zeros(1, 2*numel(tx_vec_padded_A));
    tx_vec_2x_A(1:2:end) = tx_vec_padded_A;
    tx_vec_air_A = filter(interp_filt2, 1, tx_vec_2x_A);
    
    tx_vec_2x_B = zeros(1, 2*numel(tx_vec_padded_B));
    tx_vec_2x_B(1:2:end) = tx_vec_padded_B;
    tx_vec_air_B = filter(interp_filt2, 1, tx_vec_2x_B);
    
    tx_vec_2x_C = zeros(1, 2*numel(tx_vec_padded_C));
    tx_vec_2x_C(1:2:end) = tx_vec_padded_C;
    tx_vec_air_C = filter(interp_filt2, 1, tx_vec_2x_C); 

    tx_vec_2x_D = zeros(1, 2*numel(tx_vec_padded_D));
    tx_vec_2x_D(1:2:end) = tx_vec_padded_D;
    tx_vec_air_D = filter(interp_filt2, 1, tx_vec_2x_D); 

end

% Scale the Tx vector to +/- 1
tx_vec_air_A = TX_SCALE .* tx_vec_air_A ./ max(abs(tx_vec_air_A));
tx_vec_air_B = TX_SCALE .* tx_vec_air_B ./ max(abs(tx_vec_air_B));
tx_vec_air_C = TX_SCALE .* tx_vec_air_C ./ max(abs(tx_vec_air_C));
tx_vec_air_D = TX_SCALE .* tx_vec_air_D ./ max(abs(tx_vec_air_D));


%% apply_ch 64 4 4
subs = 1;
subt = 3;
temp_ch = zeros(subt,1);
weig_ch = zeros(4,4,subt*2);
weig_ch_64 = zeros(4,4,64);
weig_ch_64_f = zeros(4,4,64);

for rx=1:4
    for tx=1:4
        temp=ifft(apply_ch(:,rx,tx));
        tep_ch = temp(subs:subt,1);%only take first 16 taps of time domain channel
        weig_ch(rx,tx,:) = upsample(tep_ch,2); % upsample channel
        weig_ch_64_f(rx,tx,:) = fft(tep_ch, 64);
    end
end

%% applying wireless channel in time domain
rx_vec_air_A = conv(tx_vec_air_A, squeeze(weig_ch(1,1,:)), 'same')+conv(tx_vec_air_B, squeeze(weig_ch(1,2,:)), 'same')+conv(tx_vec_air_C, squeeze(weig_ch(1,3,:)), 'same')+conv(tx_vec_air_D, squeeze(weig_ch(1,4,:)), 'same');
rx_vec_air_B = conv(tx_vec_air_A, squeeze(weig_ch(2,1,:)), 'same')+conv(tx_vec_air_B, squeeze(weig_ch(2,2,:)), 'same')+conv(tx_vec_air_C, squeeze(weig_ch(2,3,:)), 'same')+conv(tx_vec_air_D, squeeze(weig_ch(2,4,:)), 'same');
rx_vec_air_C = conv(tx_vec_air_A, squeeze(weig_ch(3,1,:)), 'same')+conv(tx_vec_air_B, squeeze(weig_ch(3,2,:)), 'same')+conv(tx_vec_air_C, squeeze(weig_ch(3,3,:)), 'same')+conv(tx_vec_air_D, squeeze(weig_ch(3,4,:)), 'same');
rx_vec_air_D = conv(tx_vec_air_A, squeeze(weig_ch(4,1,:)), 'same')+conv(tx_vec_air_B, squeeze(weig_ch(4,2,:)), 'same')+conv(tx_vec_air_C, squeeze(weig_ch(4,3,:)), 'same')+conv(tx_vec_air_D, squeeze(weig_ch(4,4,:)), 'same');

   

%% Decimate
if(DECIMATE_RATE == 1)
    raw_rx_dec_A = rx_vec_air_A;
    raw_rx_dec_B = rx_vec_air_B;
    raw_rx_dec_C = rx_vec_air_C;
    raw_rx_dec_D = rx_vec_air_D;
elseif(DECIMATE_RATE == 2)
    raw_rx_dec_A = filter(interp_filt2, 1, rx_vec_air_A);
    raw_rx_dec_A = raw_rx_dec_A(1:2:end);
    raw_rx_dec_B = filter(interp_filt2, 1, rx_vec_air_B);
    raw_rx_dec_B = raw_rx_dec_B(1:2:end);
    raw_rx_dec_C = filter(interp_filt2, 1, rx_vec_air_C);
    raw_rx_dec_C = raw_rx_dec_C(1:2:end);
    raw_rx_dec_D = filter(interp_filt2, 1, rx_vec_air_D);
    raw_rx_dec_D = raw_rx_dec_D(1:2:end);

end

%% Correlate for LTS

% For simplicity, we'll only use RFA for LTS correlation and peak
% discovery. A straightforward addition would be to repeat this process for
% RFB and combine the results for detection diversity.

% Complex cross correlation of Rx waveform with time-domain LTS
lts_corr = abs(conv(conj(fliplr(lts_t)), sign(raw_rx_dec_A)));

% Skip early and late samples - avoids occasional false positives from pre-AGC samples
lts_corr = lts_corr(32:end-32);

% Find all correlation peaks
lts_peaks = find(lts_corr > LTS_CORR_THRESH*max(lts_corr));

% Select best candidate correlation peak as LTS-payload boundary
% In this MIMO example, we actually have 3 LTS symbols sent in a row.
% The first two are sent by RFA on the TX node and the last one was sent
% by RFB. We will actually look for the separation between the first and the
% last for synchronizing our starting index.

[LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
[lts_last_peak_index,y] = find(LTS2-LTS1 == length(lts_t));

% Stop if no valid correlation peak was found
if(isempty(lts_last_peak_index))
    fprintf('No LTS Correlation Peaks Found!\n');
    return;
end

% Set the sample indices of the payload symbols and preamble
% The "+32" here corresponds to the 32-sample cyclic prefix on the preamble LTS
% The "+192" corresponds to the length of the extra training symbols for MIMO channel estimation
mimo_training_ind = lts_peaks(max(lts_last_peak_index))+32;


payload_ind = mimo_training_ind + 192 + 96 + 96;


% Subtract of 2 full LTS sequences and one cyclic prefixes
% The "-160" corresponds to the length of the preamble LTS (2.5 copies of 64-sample LTS)
lts_ind = mimo_training_ind-160;

if(DO_APPLY_CFO_CORRECTION)
    %Extract LTS (not yet CFO corrected)
    rx_lts = raw_rx_dec_A(lts_ind : lts_ind+159); %Extract the first two LTS for CFO
    rx_lts1 = rx_lts(-64+-FFT_OFFSET + [97:160]);
    rx_lts2 = rx_lts(-FFT_OFFSET + [97:160]);

    %Calculate coarse CFO est
    rx_cfo_est_lts = mean(unwrap(angle(rx_lts2 .* conj(rx_lts1))));
    rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*64);
else
    rx_cfo_est_lts = 0;
end

% Apply CFO correction to raw Rx waveforms
rx_cfo_corr_t = exp(-1i*2*pi*rx_cfo_est_lts*[0:length(raw_rx_dec_A)-1]);
rx_dec_cfo_corr_A = raw_rx_dec_A .* rx_cfo_corr_t;
rx_dec_cfo_corr_B = raw_rx_dec_B .* rx_cfo_corr_t;
rx_dec_cfo_corr_C = raw_rx_dec_C .* rx_cfo_corr_t;
rx_dec_cfo_corr_D = raw_rx_dec_D .* rx_cfo_corr_t;



%MIMO Channel Estimatation
% lts_ind_TXA_start = mimo_training_ind + 32 - FFT_OFFSET;
% lts_ind_TXA_end = lts_ind_TXA_start + 64 - 1;
% 
% lts_ind_TXB_start = mimo_training_ind + 32 + 64 + 32 - FFT_OFFSET;
% lts_ind_TXB_end = lts_ind_TXB_start + 64 - 1;
% 
% lts_ind_TXC_start = mimo_training_ind +32 + 64 + 32 - FFT_OFFSET;
% lts_ind_TXC_end = lts_ind_TXC_start + 64 - 1;


lts_ind_TXA_start = mimo_training_ind + 32 - FFT_OFFSET;
lts_ind_TXA_end = lts_ind_TXA_start + 64 - 1;

lts_ind_TXB_start = mimo_training_ind + 32 + 64 + 32 - FFT_OFFSET;
lts_ind_TXB_end = lts_ind_TXB_start + 64 - 1;

lts_ind_TXC_start = mimo_training_ind + 32 + 64 + 32 + 64 + 32 - FFT_OFFSET;
lts_ind_TXC_end = lts_ind_TXC_start + 64 - 1;


lts_ind_TXD_start = mimo_training_ind + 32 + 64 + 32 + 64 + 32 + 64 + 32 - FFT_OFFSET;
lts_ind_TXD_end = lts_ind_TXD_start + 64 - 1;


rx_lts_AA = rx_dec_cfo_corr_A( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BA = rx_dec_cfo_corr_A( lts_ind_TXB_start:lts_ind_TXB_end );
rx_lts_CA = rx_dec_cfo_corr_A( lts_ind_TXC_start:lts_ind_TXC_end );
rx_lts_DA = rx_dec_cfo_corr_A( lts_ind_TXD_start:lts_ind_TXD_end );



rx_lts_AB = rx_dec_cfo_corr_B( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BB = rx_dec_cfo_corr_B( lts_ind_TXB_start:lts_ind_TXB_end );
rx_lts_CB = rx_dec_cfo_corr_B( lts_ind_TXC_start:lts_ind_TXC_end );
rx_lts_DB = rx_dec_cfo_corr_B( lts_ind_TXD_start:lts_ind_TXD_end );


rx_lts_AC = rx_dec_cfo_corr_C( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BC = rx_dec_cfo_corr_C( lts_ind_TXB_start:lts_ind_TXB_end );
rx_lts_CC = rx_dec_cfo_corr_C( lts_ind_TXC_start:lts_ind_TXC_end );
rx_lts_DC = rx_dec_cfo_corr_C( lts_ind_TXD_start:lts_ind_TXD_end );



rx_lts_AD = rx_dec_cfo_corr_D( lts_ind_TXA_start:lts_ind_TXA_end );
rx_lts_BD = rx_dec_cfo_corr_D( lts_ind_TXB_start:lts_ind_TXB_end );
rx_lts_CD = rx_dec_cfo_corr_D( lts_ind_TXC_start:lts_ind_TXC_end );
rx_lts_DD = rx_dec_cfo_corr_D( lts_ind_TXD_start:lts_ind_TXD_end );



rx_lts_AA_f = fft(rx_lts_AA);
rx_lts_BA_f = fft(rx_lts_BA);
rx_lts_CA_f = fft(rx_lts_CA);
rx_lts_DA_f = fft(rx_lts_DA);


rx_lts_AB_f = fft(rx_lts_AB);
rx_lts_BB_f = fft(rx_lts_BB);
rx_lts_CB_f = fft(rx_lts_CB);
rx_lts_DB_f = fft(rx_lts_DB);


rx_lts_AC_f = fft(rx_lts_AC);
rx_lts_BC_f = fft(rx_lts_BC);
rx_lts_CC_f = fft(rx_lts_CC);
rx_lts_DC_f = fft(rx_lts_DC);


rx_lts_AD_f = fft(rx_lts_AD);
rx_lts_BD_f = fft(rx_lts_BD);
rx_lts_CD_f = fft(rx_lts_CD);
rx_lts_DD_f = fft(rx_lts_DD);


% Calculate channel estimate
rx_H_est_AA = lts_f .* rx_lts_AA_f;
rx_H_est_BA = lts_f .* rx_lts_BA_f;
rx_H_est_CA = lts_f .* rx_lts_CA_f;
rx_H_est_DA = lts_f .* rx_lts_DA_f;



rx_H_est_AB = lts_f .* rx_lts_AB_f;
rx_H_est_BB = lts_f .* rx_lts_BB_f;
rx_H_est_CB = lts_f .* rx_lts_CB_f;
rx_H_est_DB = lts_f .* rx_lts_DB_f;



rx_H_est_AC = lts_f .* rx_lts_AC_f;
rx_H_est_BC = lts_f .* rx_lts_BC_f;
rx_H_est_CC = lts_f .* rx_lts_CC_f;
rx_H_est_DC = lts_f .* rx_lts_DC_f;


rx_H_est_AD = lts_f .* rx_lts_AD_f;
rx_H_est_BD = lts_f .* rx_lts_BD_f;
rx_H_est_CD = lts_f .* rx_lts_CD_f;
rx_H_est_DD = lts_f .* rx_lts_DD_f;


%% Rx payload processing

% Extract the payload samples (integral number of OFDM symbols following preamble)
payload_vec_A = rx_dec_cfo_corr_A(payload_ind : payload_ind+(N_OFDM_SYMS/4)*(N_SC+CP_LEN)-1);
payload_mat_A = reshape(payload_vec_A, (N_SC+CP_LEN), (N_OFDM_SYMS/4));

payload_vec_B = rx_dec_cfo_corr_B(payload_ind : payload_ind+(N_OFDM_SYMS/4)*(N_SC+CP_LEN)-1);
payload_mat_B = reshape(payload_vec_B, (N_SC+CP_LEN), (N_OFDM_SYMS/4));

payload_vec_C = rx_dec_cfo_corr_C(payload_ind : payload_ind+(N_OFDM_SYMS/4)*(N_SC+CP_LEN)-1);
payload_mat_C = reshape(payload_vec_C, (N_SC+CP_LEN), (N_OFDM_SYMS/4));


payload_vec_D = rx_dec_cfo_corr_D(payload_ind : payload_ind+(N_OFDM_SYMS/4)*(N_SC+CP_LEN)-1);
payload_mat_D = reshape(payload_vec_D, (N_SC+CP_LEN), (N_OFDM_SYMS/4));


% Remove the cyclic prefix, keeping FFT_OFFSET samples of CP (on average)
payload_mat_noCP_A = payload_mat_A(CP_LEN-FFT_OFFSET+[1:N_SC], :);
payload_mat_noCP_B = payload_mat_B(CP_LEN-FFT_OFFSET+[1:N_SC], :);
payload_mat_noCP_C = payload_mat_C(CP_LEN-FFT_OFFSET+[1:N_SC], :);
payload_mat_noCP_D = payload_mat_D(CP_LEN-FFT_OFFSET+[1:N_SC], :);


% Take the FFT
syms_f_mat_A = fft(payload_mat_noCP_A, N_SC, 1);
syms_f_mat_B = fft(payload_mat_noCP_B, N_SC, 1);
syms_f_mat_C = fft(payload_mat_noCP_C, N_SC, 1);
syms_f_mat_D = fft(payload_mat_noCP_D, N_SC, 1);



% Equalize pilots
% Because we only used Tx RFA to send pilots, we can do SISO equalization
% here. This is zero-forcing (just divide by chan estimates)
syms_eq_mat_pilots = syms_f_mat_A ./ repmat(rx_H_est_AA.', 1, N_OFDM_SYMS/4);

if DO_APPLY_SFO_CORRECTION
    % SFO manifests as a frequency-dependent phase whose slope increases
    % over time as the Tx and Rx sample streams drift apart from one
    % another. To correct for this effect, we calculate this phase slope at
    % each OFDM symbol using the pilot tones and use this slope to
    % interpolate a phase correction for each data-bearing subcarrier.

	% Extract the pilot tones and "equalize" them by their nominal Tx values
    pilots_f_mat = syms_eq_mat_pilots(SC_IND_PILOTS,:);
    pilots_f_mat_comp = pilots_f_mat.*pilots_mat_A;

	% Calculate the phases of every Rx pilot tone
    pilot_phases = unwrap(angle(fftshift(pilots_f_mat_comp, 1)), [], 1);

	pilot_spacing_mat = repmat(mod(diff(fftshift(SC_IND_PILOTS)), 64).', 1, N_OFDM_SYMS/4);
    pilot_slope_mat = mean(diff(pilot_phases) ./ pilot_spacing_mat);

	% Calculate the SFO correction phases for each OFDM symbol
    pilot_phase_sfo_corr = fftshift((-32:31).' * pilot_slope_mat, 1);
    pilot_phase_corr = exp(-1i*(pilot_phase_sfo_corr));

    % Apply the pilot phase correction per symbol
    syms_f_mat_A = syms_f_mat_A .* pilot_phase_corr;
    syms_f_mat_B = syms_f_mat_B .* pilot_phase_corr;
    syms_f_mat_C = syms_f_mat_C .* pilot_phase_corr;
    syms_f_mat_D = syms_f_mat_D .* pilot_phase_corr;

    
else
	% Define an empty SFO correction matrix (used by plotting code below)
    pilot_phase_sfo_corr = zeros(N_SC, N_OFDM_SYMS);
end

% Extract the pilots and calculate per-symbol phase error
if DO_APPLY_PHASE_ERR_CORRECTION
    pilots_f_mat = syms_eq_mat_pilots(SC_IND_PILOTS, :);
    pilot_phase_err = angle(mean(pilots_f_mat.*pilots_mat_A));
else
	% Define an empty phase correction vector (used by plotting code below)
    pilot_phase_err = zeros(1, N_OFDM_SYMS/4);
end
pilot_phase_corr = repmat(exp(-1i*pilot_phase_err), N_SC, 1);

% Apply pilot phase correction to both received streams
syms_f_mat_pc_A = syms_f_mat_A .* pilot_phase_corr;
syms_f_mat_pc_B = syms_f_mat_B .* pilot_phase_corr;
syms_f_mat_pc_C = syms_f_mat_C .* pilot_phase_corr;
syms_f_mat_pc_D = syms_f_mat_D .* pilot_phase_corr;


% MIMO Equalization
% We need to apply the MIMO equalization to each subcarrier separately.
% There, unfortunately, is no great vector-y solution to do this, so we
% reluctantly employ a FOR loop.

syms_eq_mat_A = zeros(N_SC, N_OFDM_SYMS/4);
syms_eq_mat_B = zeros(N_SC, N_OFDM_SYMS/4);
syms_eq_mat_C = zeros(N_SC, N_OFDM_SYMS/4);
syms_eq_mat_D = zeros(N_SC, N_OFDM_SYMS/4);


channel_condition_mat = zeros(1,N_SC);
for sc_idx = [SC_IND_DATA, SC_IND_PILOTS]
   y = [syms_f_mat_pc_A(sc_idx,:) ; syms_f_mat_pc_B(sc_idx,:); syms_f_mat_pc_C(sc_idx,:); syms_f_mat_pc_D(sc_idx,:)];
   H = [rx_H_est_AA(sc_idx), rx_H_est_BA(sc_idx), rx_H_est_CA(sc_idx), rx_H_est_DA(sc_idx);
        rx_H_est_AB(sc_idx), rx_H_est_BB(sc_idx), rx_H_est_CB(sc_idx), rx_H_est_DB(sc_idx);  
        rx_H_est_AC(sc_idx), rx_H_est_BC(sc_idx), rx_H_est_CC(sc_idx), rx_H_est_DC(sc_idx);
        rx_H_est_AD(sc_idx), rx_H_est_BD(sc_idx), rx_H_est_CD(sc_idx), rx_H_est_DD(sc_idx)
        ];

   ch_sc(:,:, sc_idx) = H;
%    x = inv(H'*H)*(H')*y;
   x = inv(H)*y;
   syms_eq_mat_A(sc_idx, :) = x(1,:);
   syms_eq_mat_B(sc_idx, :) = x(2,:);
   syms_eq_mat_C(sc_idx, :) = x(3,:);
   syms_eq_mat_D(sc_idx, :) = x(4,:);


   channel_condition_mat(sc_idx) = rcond(H);

end


payload_syms_mat_A = syms_eq_mat_A(SC_IND_DATA, :);
payload_syms_mat_B = syms_eq_mat_B(SC_IND_DATA, :);
payload_syms_mat_C = syms_eq_mat_C(SC_IND_DATA, :);
payload_syms_mat_D = syms_eq_mat_D(SC_IND_DATA, :);


%% Demodulate
rx_syms_A = reshape(payload_syms_mat_A, 1, N_DATA_SYMS/4);
rx_syms_B = reshape(payload_syms_mat_B, 1, N_DATA_SYMS/4);
rx_syms_C = reshape(payload_syms_mat_C, 1, N_DATA_SYMS/4);
rx_syms_D = reshape(payload_syms_mat_D, 1, N_DATA_SYMS/4);


% Combine both streams to a single vector of symbols
rx_syms_space_mat = [rx_syms_A; rx_syms_B; rx_syms_C; rx_syms_D];
rx_syms = reshape(rx_syms_space_mat, 1, length(rx_syms_A)*4);

demod_fcn_bpsk = @(x) double(real(x)>0);
demod_fcn_qpsk = @(x) double(2*(real(x)>0) + 1*(imag(x)>0));
demod_fcn_16qam = @(x) (8*(real(x)>0)) + (4*(abs(real(x))<0.6325)) + (2*(imag(x)>0)) + (1*(abs(imag(x))<0.6325));
demod_fcn_64qam = @(x) (32*(real(x)>0)) + (16*(abs(real(x))<0.6172)) + (8*((abs(real(x))<(0.9258))&&((abs(real(x))>(0.3086))))) + (4*(imag(x)>0)) + (2*(abs(imag(x))<0.6172)) + (1*((abs(imag(x))<(0.9258))&&((abs(imag(x))>(0.3086)))));

switch(MOD_ORDER)
    case 2         % BPSK
        rx_data = arrayfun(demod_fcn_bpsk, rx_syms);
    case 4         % QPSK
        rx_data = arrayfun(demod_fcn_qpsk, rx_syms);
    case 8 
        refconst = qammod(0:7, 8);
        nf = modnorm(refconst, 'peakpow', 1)*1.34;
        rx_data = qamdemod(rx_syms./nf, 8);      
    case 16        % 16-QAM
        rx_data = arrayfun(demod_fcn_16qam, rx_syms);
    case 64        % 64-QAM
          rx_data = arrayfun(demod_fcn_64qam, rx_syms); 
    case 128       %128-QAM
        refconst = qammod(0:127, 128);
        nf = modnorm(refconst,'peakpow',1)*1.34;
        rx_data = qamdemod(rx_syms./nf, 128);
    case 256       %256-QAM
        refconst = qammod(0:255, 256);
        nf = modnorm(refconst,'peakpow',1)*1.34;
        rx_data = qamdemod(rx_syms./nf, 256); 
     case 512       %512-QAM
        refconst = qammod(0:511, 512);
        nf = modnorm(refconst,'peakpow',1)*1.34;
        rx_data = qamdemod(rx_syms./nf, 512);  
     case 1024       %1024-QAM
        refconst = qammod(0:1023, 1024);
        nf = modnorm(refconst,'peakpow',1)*1.34;
        rx_data = qamdemod(rx_syms./nf, 1024);     
end

ch_f = ch_sc(:,:, SC_IND_DATA); % 4 4 48 



%% without mirage
figure(1)
temps_gt = squeeze(apply_ch(SC_IND_DATA,:,1));%
[DP1, aoa_pred_d, tof_pred_d]= plot_2dfft(freq, temps_gt, theta_vals, d_vals, ant_sep, gt_aoa_d, gt_aoa_r, gt_d_val_d, gt_d_val_r, delay);
plot_transform_profile(DP1, theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], [])
%title(sprintf('without mirage-2dfft: direct-link AoA is %0.1f deg, ToF is %0.1f m; reflected-link AoA is %0.1f deg, ToF is %0.1f m; delayed %0.1f m', rad2deg(gt_aoa_d), gt_d_val_d, rad2deg(gt_aoa_r), gt_d_val_r, delay))
title('2dfft without mirage ground truth')




figure(2)
[P, aoa_pred_d, tof_pred_d]=plot_spotfi(fftshift(temps_gt, 1), theta_vals, d_vals, n_sub);
plot_transform_profile(db(P), theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], []);
title('spotfi without mirage ground truth')


%% with mirage
figure(3)
temps = squeeze(ch_f(:,1,:));
[DP1, aoa_pred_d, tof_pred_d]= plot_2dfft(freq, temps.', theta_vals, d_vals, ant_sep, gt_aoa_d, gt_aoa_r, gt_d_val_d, gt_d_val_r, delay);
plot_transform_profile(DP1, theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], [])
title('2d fft with mirage')

figure(4)
temps = squeeze(ch_f(:,1,:));
[P, aoa_pred_d, tof_pred_d]=plot_spotfi(fftshift(temps.', 1), theta_vals, d_vals, n_sub);
plot_transform_profile(db(P), theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], []);
title('spotfi with mirage')


figure(5)
temps = squeeze(ch_f(:,1,:));
[DP1, aoa_pred_d, tof_pred_d]= plot_2dfft(freq, temps.', theta_vals, d_vals, ant_sep, gt_aoa_d, gt_aoa_r, gt_d_val_d, gt_d_val_r, delay);
plot_transform_profile(DP1, theta_vals, d_vals, [], [], aoa_pred_d, tof_pred_d, [], [])
title('2d fft without mirage - mirage')




debug = false;
if debug

  save('workspace.mat');
  mimo_4x4_plot();


%% Calculate Rx stats

sym_errs = sum(tx_data ~= rx_data);
bit_errs = length(find(dec2bin(bitxor(tx_data, rx_data),8) == '1'));
rx_evm   = sqrt(sum((real(rx_syms) - real(tx_syms)).^2 + (imag(rx_syms) - imag(tx_syms)).^2)/(length(SC_IND_DATA) * N_OFDM_SYMS));


fprintf('\nResults:\n');
fprintf('Num Bytes:   %d\n', N_DATA_SYMS * log2(MOD_ORDER) / 8);
fprintf('Sym Errors:  %d (of %d total symbols)\n', sym_errs, N_DATA_SYMS);

fprintf('Bit Errors:  %d (of %d total bits)\n', bit_errs, N_DATA_SYMS * log2(MOD_ORDER));

cfo_est_lts = rx_cfo_est_lts*(SAMP_FREQ/INTERP_RATE);
cfo_est_phaseErr = mean(diff(unwrap(pilot_phase_err)))/(4e-6*2*pi);
cfo_total_ppm = ((cfo_est_lts + cfo_est_phaseErr) /  ((2.412+(.005*(CHANNEL-1)))*1e9)) * 1e6;

fprintf('CFO Est:     %3.2f kHz (%3.2f ppm)\n', (cfo_est_lts + cfo_est_phaseErr)*1e-3, cfo_total_ppm);
fprintf('     LTS CFO Est:                  %3.2f kHz\n', cfo_est_lts*1e-3);
fprintf('     Phase Error Residual CFO Est: %3.2f kHz\n', cfo_est_phaseErr*1e-3);

if DO_APPLY_SFO_CORRECTION
    drift_sec = pilot_slope_mat / (2*pi*312500);
    sfo_est_ppm =  1e6*mean((diff(drift_sec) / 4e-6));
    sfo_est = sfo_est_ppm*20;
    fprintf('SFO Est:     %3.2f Hz (%3.2f ppm)\n', sfo_est, sfo_est_ppm);

end
end
