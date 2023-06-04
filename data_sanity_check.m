ANGLE_STEP = deg2rad(1);%0.01;
D_STEP = 0.1;
THETA_VALS = -deg2rad(90):ANGLE_STEP:deg2rad(90); % aoa values for FFT
PHI_VALS = -deg2rad(90):ANGLE_STEP:deg2rad(90);
D_VALS = -10:D_STEP:50;
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

%% Define the Model of the Space
datasets = {'jacobs_July28','jacobs_July28_2','jacobs_Aug16_1','jacobs_Aug16_2','jacobs_Aug16_3','jacobs_Aug16_4_ref'};
dataset_idx = 2;
dataset_name = ['dataset_edit_',datasets{dataset_idx},'.mat'];
load(fullfile('/media/datadisk/Roshan/datasets/quantenna/features',dataset_name),...
    'labels','features_w_offset','features_wo_offset','labels_gaussian_2d','opt');
sim_data_folder = fullfile('/media/datadisk/Roshan/datasets/DLoc_sim_data/features/',datasets{dataset_idx});

DATA_SAVE_TOP = '/media/datadisk/Roshan/datasets/DLoc_sim_data';
chan_save_name = ['channels/channels_sim_',datasets{dataset_idx},'.mat'];
load(fullfile(DATA_SAVE_TOP,chan_save_name),...
    'd1','d2');


n_points = size(labels,1);
data_points = 1:100:n_points;
data_points_sim = 1:400:4*n_points;
[~,n_ap,n_d2,n_d1] = size(features_w_offset);

real_features_w_offset = features_w_offset(data_points,:,:,:);
real_features_wo_offset = features_wo_offset(data_points,:,:,:);
real_labels_gaussian_2d = labels_gaussian_2d(data_points,:,:);
real_labels = labels(data_points,:);

clear labels features_w_offset features_wo_offset labels_gaussian_2d

sim_features_w_offset = zeros(length(data_points),n_ap,n_d2,n_d1);
sim_features_wo_offset = zeros(length(data_points),n_ap,n_d2,n_d1);
sim_labels_gaussian_2d = zeros(length(data_points),n_d2,n_d1);
sim_labels = zeros(length(data_points),2);

n = 1; 

for i = data_points
    sim_filename = fullfile(sim_data_folder,[num2str(i),'.h5']);
    sim_features_w_offset(n,:,:,:) = squeeze(h5read(sim_filename,'/features_w_offset'));
    sim_features_wo_offset(n,:,:,:) = squeeze(h5read(sim_filename,'/features_wo_offset'));
    sim_labels_gaussian_2d(n,:,:) = squeeze(h5read(sim_filename,'/labels_gaussian_2d'));
    sim_labels(n,:) = squeeze(h5read(sim_filename,'/labels'));
    n = n+1;
end


%% 

figure(1)
for i = 1:length(data_points)
    for j=1:4
        subplot(4,4,j), imagesc(d1,d2,squeeze(real_features_w_offset(i,j,:,:))),axis xy,...
            hold on, scatter(real_labels(i,1), real_labels(i,2),'r','filled'), hold off,...
            xlabel('X (m)'), xlabel('Y (m)'),...
            title(['Real features with offsets for AP',num2str(j)]);
        subplot(4,4,4+j), imagesc(d1,d2,squeeze(real_features_wo_offset(i,j,:,:))),axis xy,...
            hold on, scatter(real_labels(i,1), real_labels(i,2),'r','filled'), hold off,...
            xlabel('X (m)'), xlabel('Y (m)'),...
            title(['Real features without offsets for AP',num2str(j)]);
        subplot(4,4,8+j), imagesc(d1,d2,squeeze(sim_features_w_offset(i,j,:,:))),axis xy,...
            hold on, scatter(real_labels(i,1), real_labels(i,2),'r','filled'), hold off,...
            xlabel('X (m)'), xlabel('Y (m)'),...
            title(['Sim features with offsets for AP',num2str(j)]);
        subplot(4,4,12+j), imagesc(d1,d2,squeeze(sim_features_wo_offset(i,j,:,:))),axis xy,...
            hold on, scatter(real_labels(i,1), real_labels(i,2),'r','filled'), hold off,...
            xlabel('X (m)'), xlabel('Y (m)'),...
            title(['Sim features without offsets for AP',num2str(j)]);
    end
    waitforbuttonpress;
end
        
  

