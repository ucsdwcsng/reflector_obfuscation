clearvars
close all;
%% Description: 
% The goal of this code is to get channels given a physical space. 
% The definition of physical space can include finite reflectors and 
% bounded obstacles.
% Reflectors are obstacles by default
% Color code: Red is for walls, blue is for obstacles and green is for
% reflectors
% For paths, red path is a non-blocked path and black path is a blocked
% path
%% parameter definitions
ANGLE_STEP = deg2rad(1);%0.01;
D_STEP = 0.1;
THETA_VALS = -deg2rad(90):ANGLE_STEP:deg2rad(90); % aoa values for FFT
PHI_VALS = -deg2rad(90):ANGLE_STEP:deg2rad(90);
D_VALS = -10:D_STEP:80;
C = 3e8;
BW = 80e6;
Ts = 1/BW;
CHAN_NUM = 155;
N_SUB = 256;
SUB_INDICES = [-122:-104,-102:-76,-74:-40,-38:-12,-10:-2,2:10,12:38,40:74,76:102,104:122];% 80MHz
NOISE_VAR = 0.26;
OUTPUT_SIGMA = 0.25;
% define opt
opt.freq = double(5e9 + 5e6*CHAN_NUM) + (BW/N_SUB).*SUB_INDICES;%6489.6e6;
opt.lambda = C./opt.freq;
opt.ant_sep = min(opt.lambda)/2;
opt.threshold = 0.2;
zeroFreq = mean(opt.freq);
opt.ant_sep = 0.026;
t_res = 1;
n_freq = length(SUB_INDICES);
%% Define the Model of the Space
side_length = 15;

walls = get_rectangle([0,0],40,30);  % Walls are defined just for plotting. They are non-reflecting, non-blocking                              

ap_centers = [20,0.5;39.5,15;20,29.5;0.5,15];
ap_axis = [1,2,1,2];
ap_orient = [1,1,-1,-1];
ap_aoas = [pi/2,pi,-pi/2,0];
n_aps = size(ap_centers,1);
n_ant_per_ap = 4;
aps = zeros(n_aps,n_ant_per_ap,2);
ant_array = (-(n_ant_per_ap-1)/2:1:(n_ant_per_ap-1)/2).'.*opt.ant_sep;
for i = 1:n_aps
    if ap_axis(i) == 1
        aps(i,:,:) = ap_centers(i,:) + ap_orient(i).*[ant_array, zeros(n_ant_per_ap,1)];
    elseif ap_axis(i) ==2
        aps(i,:,:) = ap_centers(i,:) + ap_orient(i).*[zeros(n_ant_per_ap,1), ant_array];
    end
end

obstacles{1} = get_rectangle([2,22],8,7); %Obstacles are bounded and opaque to wireless signals
% obstacles{2} = get_rectangle([5,5],30,21);
reflectors{1} = [0,0;40,0]; % Reflectors are linear for ease 
reflectors{2} = [0,30;40,30]; % Reflectors are linear for ease 
reflectors{3} = [0,0;0,30]; % Reflectors are linear for ease 
reflectors{4} = [40,0;40,30]; % Reflectors are linear for ease 
reflectors{5} = [30,27;35,20]; % Reflectors are linear for ease 

model.walls = walls;
model.reflectors =  reflectors;
model.obstacles = obstacles;
model.obs_attenuation = 0.7;
model.ref_attenuation = 1-model.obs_attenuation;
model.lambda = opt.lambda;
model.freq = opt.freq;
model.amps = ones(n_freq,1)*50; % Lambda and amps can be arrays

figure(1), display_model(model), hold on
for i = 1:n_aps
    scatter(aps(i,:,1), aps(i,:,2),'k+');
end
figure(1), hold off;
%% Get variables
% n_total = size(labels,1);
labels_sup = zeros(1000,2);
aoas_gt = zeros(1000,2);
aoas_pred = zeros(1000,2);
aoas_pred_new = zeros(1000,2);
for n_exp = 1:100

n_total = 10;
user_bbox = get_rectangle([5,5],30,21);

%% Generate the points. Ap 1 antenna 1 serves as the reference
channels_user = zeros(n_total,n_aps, n_freq,n_ant_per_ap, n_ant_per_ap);
channels_ap = zeros(n_total,n_aps, n_freq,n_ant_per_ap, n_ant_per_ap);
channels_new = zeros(n_total,n_aps, n_freq,n_ant_per_ap, n_ant_per_ap);
desired_channel = zeros(n_total,n_aps, n_freq,n_ant_per_ap);
bfmg_mat = zeros(n_total,n_aps, n_freq,n_ant_per_ap, n_ant_per_ap);
channels_ap_new = zeros(n_total,n_aps, n_freq,n_ant_per_ap, n_ant_per_ap);

S = get_2dsteering_matrix(THETA_VALS,D_VALS,length(SUB_INDICES),mean(opt.freq), mean(diff(opt.freq)), opt.ant_sep);

% dst = zeros(n_total,n_ant_per_ap,2);
disp('Generating channels');
figure(2), subplot(4,4,9.5), display_model(model), hold on,...
    scatter(src(:,1),src(:,2),'r*');
xlim([-2 side_length+3]), ylim([-2 side_length+3]);
title('Setup')
n_curr = 0;

%%
labels = zeros(n_total,2);
labels_all = zeros(n_total,n_ant_per_ap, 2);
features_ap = zeros(n_total, n_aps, length(THETA_VALS), length(D_VALS));
features_ap_new = zeros(n_total, n_aps, length(THETA_VALS), length(D_VALS));

aoa_gt = zeros(n_total,n_aps);

delay = 4*(3e8/BW);
aoas_ap = zeros(n_total,n_aps);
aoas_ap_new = zeros(n_total,n_aps);
parfor n_curr =1:n_total
    curr_loc = rand(1,2).*([25,15])+[5,5];
%     if (is_in_obstacle(obstacles,curr_loc))
%         continue;
%     end
    offset = rand(1);
    labels(n_curr,:) = curr_loc;
    theta_curr = (rand()-0.5)*pi;
    labels_all(n_curr,:, :) = curr_loc + [ant_array, zeros(n_ant_per_ap,1)]*[cos(theta_curr), sin(theta_curr);-1*sin(theta_curr), cos(theta_curr)];
    [channels_user(n_curr,:,:,:,:), channels_ap(n_curr,:,:,:,:)] = get_ap_user_channels_from_model(model, aps, squeeze(labels_all(n_curr,:,:)), offset);
    [features_ap(n_curr,:,:,:), features_ap_new(n_curr,:,:,:), aoas_ap(n_curr,:), aoas_ap_new(n_curr,:)] = ...
        get_features_from_channels_mirage(squeeze(channels_user(n_curr,:,:,:,:)), squeeze(channels_ap(n_curr,:,:,:,:)), opt, delay, THETA_VALS, D_VALS, S);
%     for j = 1:n_aps
%         precoding_matrix(n_curr,j,:,:,:) = create_mirage_roshan(squeeze(channels_user(n_curr,j,:,:,:)),opt, delay);
% %         bfmg_mat(i,:,j,:) = squeeze(conj(channels_user(i,:,j,:))).*squeeze(squeeze(desired_channel(i,:,:)));
% %         for fq = 1:n_freq
%         channels_ap_new(n_curr,j,:,:,:) = squeeze(channels_ap(n_curr,j,:,:,:)).*squeeze(precoding_matrix(n_curr,j,:,:,:));
% %         end
%         features_user(n_curr,j,:,:) = compute_spotfi_profile_vectorized(squeeze(channels_user(n_curr,j,:,1,:)),THETA_VALS,D_VALS,opt,S);
%         features_ap(n_curr,j,:,:) = compute_spotfi_profile_vectorized(squeeze(channels_ap(n_curr,j,:,1,:)),THETA_VALS,D_VALS,opt,S);
%         features_ap_new(n_curr,j,:,:) = compute_spotfi_profile_vectorized(squeeze(channels_ap_new(n_curr,j,:,1,:)),THETA_VALS,D_VALS,opt,S);
%         
%         aoas_ap(n_curr,j) = rad2deg(get_aoa_for_least_tof(squeeze(features_ap(n_curr,j,:,:)),D_VALS,THETA_VALS));
%         aoas_ap_new(n_curr,j) = rad2deg(get_aoa_for_least_tof(squeeze(features_ap_new(n_curr,j,:,:)),D_VALS,THETA_VALS));
% %         features_desired(i,:,:) = compute_spotfi_profile_vectorized(squeeze(desired_channel(i,:,:)),THETA_VALS,D_VALS,opt,S);
% 
% %         features_user(i,j,:,:) = compute_multipath_profile2d_fast_edit(squeeze(channels_user(i,:,:,1)),THETA_VALS,D_VALS,opt);
% %         features_ap(i,j,:,:) = compute_multipath_profile2d_fast_edit(squeeze(channels_ap(i,:,1,:)),THETA_VALS,D_VALS,opt);
% %         features_ap_new(i,j,:,:) = compute_multipath_profile2d_fast_edit(squeeze(channels_ap_new(i,:,1,:)),THETA_VALS,D_VALS,opt);
% %         features_desired(i,:,:) = compute_multipath_profile2d_fast_edit(squeeze(desired_channel(i,:,:)),THETA_VALS,D_VALS,opt);
%     end
%     for j = 1:4
%         figure(2),
%         subplot(3,4,j), imagesc(D_VALS,THETA_VALS*180/pi,squeeze(features_user(n_curr,j,:,:))),...
%             axis xy, xlabel('ToF (m)'), ylabel('AoA (\circ)'),...
%             title(['User end profiles for AP Ant', num2str(j)]);
%         subplot(3,4,4+j), imagesc(D_VALS,THETA_VALS*180/pi,squeeze(features_ap(n_curr,j,:,:))),...
%             axis xy, xlabel('ToF (m)'), ylabel('AoA (\circ)'),...
%             title(['AP end profiles for User Ant', num2str(j)]);
%         subplot(3,4,8+j), imagesc(D_VALS,THETA_VALS*180/pi,squeeze(features_ap_new(n_curr,j,:,:))),...
%             axis xy, xlabel('ToF (m)'), ylabel('AoA (\circ)'),...
%             title(['Modified AP end profiles for user Ant', num2str(j)]);
%         subplot(4,4,11.5), imagesc(D_VALS,THETA_VALS*180/pi,squeeze(features_ap_new(i,:,:))),...
%             axis xy, xlabel('ToF (m)'), ylabel('AoA (\circ)'),...
%             title('Profile for the phase that is added at the user end');
%     end
%     for j = 1:4
%         figure(3),
%         subplot(4,4,j), imagesc(D_VALS,THETA_VALS*180/pi,squeeze(features_user_s(i,j,:,:))),...
%             axis xy, xlabel('ToF (m)'), ylabel('AoA (\circ)'),...
%             title(['User end profiles for AP Ant', num2str(j)]);
%         subplot(4,4,4+j), imagesc(D_VALS,THETA_VALS*180/pi,squeeze(features_ap_s(i,j,:,:))),...
%             axis xy, xlabel('ToF (m)'), ylabel('AoA (\circ)'),...
%             title(['AP end profiles for User Ant', num2str(j)]);
%         subplot(4,4,12+j), imagesc(D_VALS,THETA_VALS*180/pi,squeeze(features_ap_new_s(i,j,:,:))),...
%             axis xy, xlabel('ToF (m)'), ylabel('AoA (\circ)'),...
%             title(['Modified AP end profiles for user Ant', num2str(j)]);
%         subplot(4,4,[10,11]), imagesc(D_VALS,THETA_VALS*180/pi,squeeze(features_ap_new_s(i,:,:))),...
%  aoa_gt           axis xy, xlabel('ToF (m)'), ylabel('AoA (\circ)'),...
%             title('Profile for the phase that is added at the user end');
%     end
%     waitforbuttonpress;
    if(mod(n_curr,10))
        print([num2str(n_curr),'points have been processed']);
    end
end

for n_curr = 1:n_total
    for ap_idx = 1:n_aps
        ap_pos = squeeze(mean(aps(ap_idx,:,:)))';
        ap_vec=squeeze(aps(ap_idx,1,:)-aps(ap_idx,end,:))';
        X=labels(n_curr,1)-ap_pos(1);
        Y=labels(n_curr,2)-ap_pos(2);
        aoa_gt(n_curr,ap_idx)= rad2deg(sign(sum([X,Y].*ap_vec))*(pi/2-acos(abs(sum([X,Y].*ap_vec))/norm([X,Y])/norm(ap_vec))));
    end
end
aoas_gt((n_exp-1)*n_total+(1:n_total)) = aoa_gt;
labels_sup((n_exp-1)*n_total+(1:n_total),:) = labels;
aoas_pred((n_exp-1)*n_total+(1:n_total)) = aoas_ap;
aoas_pred_new((n_exp-1)*n_total+(1:n_total)) = aoas_ap_new;
end

% figure(2), subplot(4,4,9.5), hold off;
% disp('Loading Real data features')






