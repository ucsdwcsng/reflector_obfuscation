


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
opt.ant_sep = 0.026;


datasets = {'jacobs_July28','jacobs_July28_2','jacobs_Aug16_1','jacobs_Aug16_2','jacobs_Aug16_3','jacobs_Aug16_4_ref'};
dataset_idx = 1;

DATA_SAVE_TOP = '/media/datadisk/Roshan/datasets/DLoc_sim_data';
chan_save_name = ['channels/channels_sim_',datasets{dataset_idx},'.mat'];
save(fullfile(DATA_SAVE_TOP,chan_save_name),...
    'channels_w_offset','channels_wo_offset','labels',...
    'offsets','ap','d1','d2');

n_total = size(channels_w_offset,1);
features_w_offset_all = zeros(n_total, length(ap), length(d2), length(d1));
features_wo_offset_all = zeros(n_total, length(ap), length(d2), length(d1));
parfor i = 1:n_total
    features_w_offset_all(i,:,:,:) = generate_features_abs(squeeze(channels_w_offset(i,:,:,:)),ap,THETA_VALS,D_VALS,d1,d2,opt);
    features_wo_offset_all(i,:,:,:) = generate_features_abs(squeeze(channels_wo_offset(i,:,:,:)),ap,THETA_VALS,D_VALS,d1,d2,opt);
    if(mod(i,1000)==0)
        hrss = num2str(floor(rem(now -start_time,1)*24));
        minss = num2str(rem((now - start_time)*24,1)*60);
        disp(['Time taken so far for ',num2str(i),' iterations is ',hrss,'hrs ',minss,'mins']);
    end
end
labels_gaussian_2d_all = get_gaussian_labels(labels,...
        OUTPUT_SIGMA,...
        d1,...
        d2);

    %% Saving individual feature files
dataset_name = datasets{dataset_idx};
if ~exist(fullfile(DATA_SAVE_TOP,'features'), 'dir')
   mkdir(fullfile(DATA_SAVE_TOP,'features'));
end
if ~exist(fullfile(DATA_SAVE_TOP,'features',dataset_name), 'dir')
    mkdir(fullfile(DATA_SAVE_TOP,'features',dataset_name))
end
n_start = 0;
for i=1:n_total
    features_wo_offset = features_wo_offset_all(i,:,:,:);
    features_w_offset = features_w_offset_all(i,:,:,:);
    labels_gaussian_2d = labels_gaussian_2d_all(i,:,:);
    labels_discrete = labels(i,:,:);
    if (mod(i,1000)==0)
        fprintf('Saving....%d.h5\n',i+n_start);
    end
    fname = [num2str(i+n_start),'.h5'];
    
    if(~exist(fullfile(DATA_SAVE_TOP,'features',dataset_name,fname), 'file'))
        h5create(fullfile(DATA_SAVE_TOP,'features',dataset_name,fname),...
            '/features_w_offset',size(features_w_offset));
        h5create(fullfile(DATA_SAVE_TOP,'features',dataset_name,fname),...
            '/features_wo_offset',size(features_wo_offset));
        h5create(fullfile(DATA_SAVE_TOP,'features',dataset_name,fname),...
            '/labels',size(labels_discrete));
        h5create(fullfile(DATA_SAVE_TOP,'features',dataset_name,fname),...
            '/labels_gaussian_2d',size(labels_gaussian_2d));
    end

    h5write(fullfile(DATA_SAVE_TOP,'features',dataset_name,fname),...
        '/features_w_offset',features_w_offset);
    h5write(fullfile(DATA_SAVE_TOP,'features',dataset_name,fname),...
        '/features_wo_offset',features_wo_offset);
    h5write(fullfile(DATA_SAVE_TOP,'features',dataset_name,fname),...
        '/labels',labels_discrete);
    h5write(fullfile(DATA_SAVE_TOP,'features',dataset_name,fname),...
        '/labels_gaussian_2d',labels_gaussian_2d);

end
