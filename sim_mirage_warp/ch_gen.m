% clearvars
% close all;
%% Description: 
% The goal of this code is to get channels given a physical space. 
% The definition of physical space can include finite reflectors and 
% bounded obstacles.
% Reflectors are obstacles by default
% Color code: Red is for walls, blue is for obstacles and green is for
% reflectors
% For paths, red path is a non-blocked path and black path is a blocked
% path
function [channels_user, channels_ap]=ch_gen()
    %% parameter definitions
    ANGLE_STEP = deg2rad(1);%0.01;
    D_STEP = 0.1;
    THETA_VALS = -deg2rad(90):ANGLE_STEP:deg2rad(90); % aoa values for FFT
    PHI_VALS = -deg2rad(90):ANGLE_STEP:deg2rad(90);
    D_VALS = -10:D_STEP:80;
    C = 3e8;
    BW = 20e6;
    Ts = 1/BW;
    NOISE_VAR = 0.26;
    OUTPUT_SIGMA = 0.25;
    
    
    opt.ant_sep = 0.026;
    %%
    SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
    SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
    subcarrier_indices = SC_IND_DATA-32; %symmetric data subcarriers 
    freq = fftshift(5.18e9 + subcarrier_indices.*BW./64);%frequency of each data subcarrier
    subcarrier_indices_pilots = SC_IND_PILOTS - 32; % symmetric pilot subcarriers
    
    opt.freq = freq;
    opt.lambda = C./opt.freq;
    opt.ant_sep = min(opt.lambda)/2;
    opt.threshold = 0.2;
    opt.ant_sep = 0.026;
    n_freq = length(SC_IND_DATA);
    
    %% Define the Model of the Space
    walls = get_rectangle([0,0],70,70);  % Walls are defined just for plotting. They are non-reflecting, non-blocking                              
    
    ap_centers = [20,10]; % position of center of ap
    ap_orient = [0*pi/180]; % orienation of antenna array wrt positive x axis
    ap_ant_ordering = [1]; % swap the antenna ordering, ie place 1st antenna in 4th antenna pos
    n_aps = size(ap_centers,1);
    n_ant_per_ap = 4;
    aps = zeros(n_aps,n_ant_per_ap,2);
    ant_array = (-(n_ant_per_ap-1)/2:1:(n_ant_per_ap-1)/2).'.*opt.ant_sep;
    for i = 1:n_aps
        pos_curr = ap_centers(i,:);
        theta_curr = ap_orient(i);
        ordering_curr = ap_ant_ordering(i);
        aps(i,:, :) = pos_curr + [ordering_curr*ant_array, zeros(n_ant_per_ap,1)]*[cos(theta_curr), sin(theta_curr);-1*sin(theta_curr), cos(theta_curr)];
    end

%     reflectors{1} = [60-31.3848/2,0;60-31.3848/2,70]; % Reflectors are linear for ease 
%     reflectors{1} = [50, 0; 50, 70]; % Reflectors are linear for ease 
   reflectors{1} = [10, 0; 10, 70];
%    reflectors{2} = [40, 0; 40, 70];
%    reflectors{3} = [0, 20; 10, 70];


    model.walls = walls;
    model.reflectors = reflectors;
    model.obstacles = {};
    model.obs_attenuation = 1;
    model.ref_attenuation = 1; %1-model.obs_attenuation;
    model.lambda = opt.lambda;
    model.freq = opt.freq;
    model.amps = ones(n_freq,1)*50;%50; % Lambda and amps can be arrays
    
    
    figure(300)
    display_model(model), hold on
    for i = 1:n_aps
        scatter(aps(i,:,1), aps(i,:,2),'k+');
    end

    for n_exp = 1:1
        n_total = 1;%10;        
        %% Generate the points. Ap 1 antenna 1 serves as the reference
        channels_user = zeros(n_total,n_aps, n_freq,n_ant_per_ap, n_ant_per_ap);
        channels_ap = zeros(n_total,n_aps, n_freq,n_ant_per_ap, n_ant_per_ap);
        
        labels = zeros(n_total,2);
        labels_all = zeros(n_total,n_ant_per_ap, 2);
        
        for n_curr =1:n_total
            curr_loc = [28,40]; % rand(1,2).*([25,15])+[5,5];
            offset = 0;
            labels(n_curr,:) = curr_loc;
            theta_curr = 0; % (rand()-0.5)*pi;
            labels_all(n_curr,:, :) = curr_loc + [ant_array, zeros(n_ant_per_ap,1)]*[cos(theta_curr), sin(theta_curr);-1*sin(theta_curr), cos(theta_curr)];
            scatter(labels_all(n_curr,:,1), labels_all(n_curr,:,2),'ko');
            [channels_user(n_curr,:,:,:,:), channels_ap(n_curr,:,:,:,:)] = get_ap_user_channels_from_model(model, aps, squeeze(labels_all(n_curr,:,:)), offset);
        end
    end
end


