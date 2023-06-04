function [rx_vec_air_A, rx_vec_air_B, rx_vec_air_C, rx_vec_air_D]=ch_cov(tx_vec_air_A, tx_vec_air_B, tx_vec_air_C, tx_vec_air_D, wireless_ch)

subs = 1;
subt = 16;
temp_ch = zeros(subt,1);
weig_ch = zeros(4,4,subt);
% weig_ch_64 = zeros(4,4,64);
% weig_ch_64_f = zeros(4,4,64);

for rx=1:4
    for tx=1:4
        temp=ifft(squeeze(wireless_ch(:,rx,tx)));
        tep_ch = temp(subs:subt,1);%only take first 16 taps of time domain channel
        weig_ch(rx,tx,:) = tep_ch;
%         weig_ch(rx,tx,:) = upsample(tep_ch,2); % upsample channel
%         weig_ch_64_f(rx,tx,:) = fft(tep_ch, 64);
    end
end

%% applying wireless channel in time domain
rx_vec_air_A = conv(tx_vec_air_A, squeeze(weig_ch(1,1,:)), 'same')+conv(tx_vec_air_B, squeeze(weig_ch(1,2,:)), 'same')+conv(tx_vec_air_C, squeeze(weig_ch(1,3,:)), 'same')+conv(tx_vec_air_D, squeeze(weig_ch(1,4,:)), 'same');
rx_vec_air_B = conv(tx_vec_air_A, squeeze(weig_ch(2,1,:)), 'same')+conv(tx_vec_air_B, squeeze(weig_ch(2,2,:)), 'same')+conv(tx_vec_air_C, squeeze(weig_ch(2,3,:)), 'same')+conv(tx_vec_air_D, squeeze(weig_ch(2,4,:)), 'same');
rx_vec_air_C = conv(tx_vec_air_A, squeeze(weig_ch(3,1,:)), 'same')+conv(tx_vec_air_B, squeeze(weig_ch(3,2,:)), 'same')+conv(tx_vec_air_C, squeeze(weig_ch(3,3,:)), 'same')+conv(tx_vec_air_D, squeeze(weig_ch(3,4,:)), 'same');
rx_vec_air_D = conv(tx_vec_air_A, squeeze(weig_ch(4,1,:)), 'same')+conv(tx_vec_air_B, squeeze(weig_ch(4,2,:)), 'same')+conv(tx_vec_air_C, squeeze(weig_ch(4,3,:)), 'same')+conv(tx_vec_air_D, squeeze(weig_ch(4,4,:)), 'same');
end