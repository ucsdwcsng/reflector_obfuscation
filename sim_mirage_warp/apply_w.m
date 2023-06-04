function [preamble_A, preamble_B, preamble_C, preamble_D]=apply_w(precoding_w, lts_f, sts_f, TX_SPATIAL_STREAM_SHIFT, sts_t, lts_t)
f_domain=false;
if f_domain
sts_f_across_antenna = zeros(64, 4,1);
lts_f_across_antenna = zeros(64, 4,4);

for sub=1:64
    sts_f_across_antenna(sub,:,:) = transpose(squeeze(precoding_w(sub, :,:)))*repmat(sts_f(1, sub), 4, 1);
%     sts_f_across_antenna(sub,:,:) = repmat(sts_f(1, sub), 4, 1);
end

for rx = 1:4
    for tx = 1:4
        lts_f_across_antenna(:,rx,tx) = squeeze(precoding_w(:, rx,tx)).*lts_f.';
%         lts_f_across_antenna(:,rx,tx) = lts_f.';
    end
end

sts_pre_A = sts_f_across_antenna(:,1,:).';% 64, 1
sts_pre_B = sts_f_across_antenna(:,2,:).';
sts_pre_C = sts_f_across_antenna(:,3,:).';
sts_pre_D = sts_f_across_antenna(:,4,:).';



sts_t_w_A_temp = ifft(sqrt(13/6).*sts_pre_A, 64); %%% checking the sts_f to be 1, 4, 64, ch_all_A is 64, 4
sts_t_w_A = sts_t_w_A_temp(1:16);
sts_t_rep_A = repmat(sts_t_w_A, 1, 30);

scale = 1;
lts_t_AA = ifft(lts_f_across_antenna(:, 1, 1), 64).';
lts_t_AB = ifft(lts_f_across_antenna(:, 2, 1), 64).';
lts_t_AC = ifft(lts_f_across_antenna(:, 3, 1), 64).';
lts_t_AD = ifft(lts_f_across_antenna(:, 4, 1), 64).';

preamble_legacy_A = [sts_t_rep_A, lts_t_AA(33:64), lts_t_AA, lts_t_AA];
preamble_mimo_A = [lts_t_AA(33:64), lts_t_AA, ...
                   scale*lts_t_AB(33:64),  scale*lts_t_AB, ... 
                    scale*lts_t_AC(33:64),  scale*lts_t_AC, ...
                    scale*lts_t_AD(33:64),  scale*lts_t_AD];


%apply weight on frequency domain for preamble on antenna B
sts_t_w_B_temp = ifft(sqrt(13/6).*sts_pre_B, 64);
sts_t_w_B = sts_t_w_B_temp(1:16);
sts_t_rep_B = repmat(sts_t_w_B, 1, 30);


lts_t_BA = ifft(lts_f_across_antenna(:, 1, 2), 64).';
lts_t_BB = ifft(lts_f_across_antenna(:, 2, 2), 64).';
lts_t_BC = ifft(lts_f_across_antenna(:, 3, 2), 64).';
lts_t_BD = ifft(lts_f_across_antenna(:, 4, 2), 64).';

preamble_legacy_B = [circshift(sts_t_rep_B, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];
preamble_mimo_B = [ scale*lts_t_BA(33:64),  scale*lts_t_BA, ...
                   lts_t_BB(33:64), lts_t_BB, ... 
                    scale*lts_t_BC(33:64),  scale*lts_t_BC, ...
                    scale*lts_t_BD(33:64),  scale*lts_t_BD];

%apply weight on frequency domain for preamble on antenna C
sts_t_w_C_temp = ifft(sqrt(13/6).*sts_pre_C, 64);
sts_t_w_C = sts_t_w_C_temp(1:16);
sts_t_rep_C = repmat(sts_t_w_C, 1, 30);
lts_t_CA = ifft(lts_f_across_antenna(:, 1, 3), 64).';
lts_t_CB = ifft(lts_f_across_antenna(:, 2, 3), 64).';
lts_t_CC = ifft(lts_f_across_antenna(:, 3, 3), 64).';
lts_t_CD = ifft(lts_f_across_antenna(:, 4, 3), 64).';

preamble_legacy_C = [circshift(sts_t_rep_C, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];
preamble_mimo_C = [ scale*lts_t_CA(33:64),  scale*lts_t_CA, ...
                    scale*lts_t_CB(33:64),  scale*lts_t_CB, ... 
                   lts_t_CC(33:64), lts_t_CC, ...
                    scale*lts_t_CD(33:64), scale* lts_t_CD];

%apply weight on frequency domain for preamble on antenna D
sts_t_w_D_temp = ifft(sqrt(13/6).*sts_pre_D, 64);
sts_t_w_D = sts_t_w_D_temp(1:16);
sts_t_rep_D = repmat(sts_t_w_D, 1, 30);
lts_t_DA = ifft(lts_f_across_antenna(:, 1, 4), 64).';
lts_t_DB = ifft(lts_f_across_antenna(:, 2, 4), 64).';
lts_t_DC = ifft(lts_f_across_antenna(:, 3, 4), 64).';
lts_t_DD = ifft(lts_f_across_antenna(:, 4, 4), 64).';

preamble_legacy_D = [circshift(sts_t_rep_D, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];
preamble_mimo_D = [ scale*lts_t_DA(33:64),  scale*lts_t_DA, ...
                    scale*lts_t_DB(33:64),  scale*lts_t_DB, ... 
                    scale*lts_t_DC(33:64),  scale*lts_t_DC, ...
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
else
sts_t_rep = repmat(sts_t, 1, 30);

preamble_legacy_A = [sts_t_rep, lts_t(33:64), lts_t, lts_t];
preamble_legacy_B = [circshift(sts_t_rep, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];
preamble_legacy_C = [circshift(sts_t_rep, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];
preamble_legacy_D = [circshift(sts_t_rep, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];


preamble_mimo_A = [lts_t(33:64), lts_t, zeros(1,96), zeros(1,96), zeros(1,96)];
preamble_mimo_B = [zeros(1,96), lts_t(33:64), lts_t, zeros(1,96), zeros(1,96)];
preamble_mimo_C = [zeros(1,96), zeros(1,96), lts_t(33:64), lts_t, zeros(1,96)];
preamble_mimo_D = [zeros(1,96), zeros(1,96), zeros(1,96), lts_t(33:64), lts_t];

end
preamble_A = [preamble_legacy_A, preamble_mimo_A];
preamble_B = [preamble_legacy_B, preamble_mimo_B];
preamble_C = [preamble_legacy_C, preamble_mimo_C];
preamble_D = [preamble_legacy_D, preamble_mimo_D];

end