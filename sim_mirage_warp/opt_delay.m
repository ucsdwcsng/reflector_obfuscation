%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wl_example_mimo_ofdm_txrx.m
% 2x2 MIMO OFDM Example
% A detailed write-up of this example is available on the wiki:
% http://warpproject.org/trac/wiki/WARPLab/Examples/MIMO_OFDM
%
% Copyright (c) 2015 Mango Communications - All Rights Reserved
% Distributed under the WARP License (http://warpproject.org/license)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear
% clc
% close all


function  avg_snr = opt_delay(delay)

% Params:
NUM_PACKETS             = 1;
USE_WARPLAB_TXRX        = 0;           % Enable WARPLab-in-the-loop (otherwise sim-only)
WRITE_PNG_FILES         = 0;           % Enable writing plots to PNG
CHANNEL                 = 11;          % Channel to tune Tx and Rx radios

% Waveform params
N_OFDM_SYMS             = 1000;        % Number of OFDM symbols (must be even valued)
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

%applied weight
% weight_f = channel();
% weight_f = rand(52, 4);
% save('random_weight.mat', 'weight_f');
% weight = load('random_weight.mat');
% weight_f = weight.weight_f;
% weights_f = ones(52, 4);

%% beamforming, nulling and delaying
n_ant = 4;
n_sub = length(SC_IND_DATA);
c = 3e8;
ant_sep = 0.0259;
theta_vals = -pi/2:0.01:pi/2; % AoA search space
d_vals = 0:0.1:100; 
BW = 20e6;
subcarrier_indices = SC_IND_DATA; %48 subcarrier
% CHANNEL=35;
freq = double(5e9 + 5*34.5*1e6) + subcarrier_indices.*BW./n_sub;
% freq = 5.1692e9+subcarrier_indices.*BW./n_sub;

% aoa = pi/2;
d_val = 10;% distance 25
r_delay = 5; 

% 
gt_aoa_d = pi/6;
gt_aoa_r = -pi/5;

gt_aod_d = gt_aoa_d;
gt_aod_r = gt_aoa_r;

gt_d_val_d = d_val;
gt_d_val_r = d_val+r_delay;

delay = 15; %10

bf_d = ones(48, 4);
bf_r = ones(48, 4);

for subcarrier_idx =1:n_sub
  for ant_dix = 1:n_ant
      f = freq(subcarrier_idx);
      bf_d(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aoa_d))/3e8);
      bf_r(subcarrier_idx,ant_dix) = exp(-1j*2*pi*f*((ant_dix-1)*ant_sep*sin(gt_aoa_r))/3e8);
  end
end

h = zeros(48,4,4);
for subcarrier_idx =1:48
    for j=1:4
        for k=1:4
            f = freq(subcarrier_idx);

            t_t_d = ((j-1) * ant_sep * sin(gt_aod_d))/3e8;
            t_t_r = ((j-1) * ant_sep * sin(gt_aod_r))/3e8;

            t_r_d = ((k-1) * ant_sep * sin(gt_aoa_d))/3e8;
            t_r_r = ((k-1) * ant_sep * sin(gt_aoa_r))/3e8;

            h_d = exp(-1j*2*pi*f*(t_t_d+t_r_d+gt_d_val_d/3e8));
            h_r = exp(-1j*2*pi*f*(t_t_r+t_r_r+gt_d_val_r/3e8));

            h(subcarrier_idx, j, k) = h_d + h_r;
        end
    end
end



%% beamforming and nulling direction
null_d = zeros(48, 4);
null_nd = zeros(48, 4);
for i=1:48
 null_d_space = null(bf_d(i,:));
 temp_w_1 = conj(bf_r(i,:))*null_d_space(:,1);
 temp_w_2 = conj(bf_r(i,:))*null_d_space(:,2);
 temp_w_3 = conj(bf_r(i,:))*null_d_space(:,3);
 null_d(i,:) = (temp_w_1*null_d_space(:,1)+temp_w_2*null_d_space(:,2)+temp_w_3*null_d_space(:,3))/(temp_w_1+temp_w_2+temp_w_3);

 null_nd_space = null(bf_r(i,:));    
 temp_w_1 = conj(bf_d(i,:))*null_nd_space(:,1);
 temp_w_2 = conj(bf_d(i,:))*null_nd_space(:,2);
 temp_w_3 = conj(bf_d(i,:))*null_nd_space(:,3);
 null_nd(i,:) = (temp_w_1*null_nd_space(:,1)+temp_w_2*null_nd_space(:,2)+temp_w_3*null_nd_space(:,3))/(temp_w_1+temp_w_2+temp_w_3);
end

bf_null_delay = ones(48, 4);
opt_ch = zeros(48,4, 4);
for subcarrier_idx =1:n_sub
  for ant_dix = 1:n_ant
      f = freq(subcarrier_idx);
      bf_null_delay(subcarrier_idx,ant_dix) = exp(1j*2*pi*f*delay/3e8)*conj(null_nd(subcarrier_idx, ant_dix))*null_nd(subcarrier_idx,ant_dix)+conj(null_d(subcarrier_idx, ant_dix))*null_d(subcarrier_idx,ant_dix);
      opt_ch(subcarrier_idx, ant_dix, :) = h(subcarrier_idx, ant_dix,:).*bf_null_delay(subcarrier_idx, ant_dix);
  end
end

%% average snr
avg_snr = 0;
for subcarrier_idx = 1:n_sub
    for j=1:4
        for k=1:4
             avg_snr = avg_snr + mag2db(abs(opt_ch(subcarrier_idx, j,k)));
        end
    end
end
avg_snr = avg_snr/(n_sub*4*4);

end


