function [channels_user, channels_ap] = get_ap_user_channels_from_model(model, aps, labels, offset)

[n_aps, n_ant_per_ap, ~] = size(aps);
n_freq = length(model.lambda);
n_ant_per_user = size(labels,1);

channels_user = zeros(n_aps,n_freq,n_ant_per_ap,n_ant_per_user);
channels_ap = zeros(n_aps,n_freq,n_ant_per_user,n_ant_per_ap);

for ap_num = 1:n_aps
    for j=1:n_ant_per_ap
        for k=1:n_ant_per_user
            channels_user(ap_num,:,j,k) = get_channels_from_model(model,squeeze(aps(ap_num, j,:)).',labels(k,:),false,offset);
            channels_user(ap_num,:,j,k) = squeeze(awgn(squeeze(channels_user(ap_num,:,j,k)),15));
            channels_ap(ap_num,:,k,j) = get_channels_from_model(model,labels(k,:),squeeze(aps(ap_num, j,:)).',false,offset);
            channels_ap(ap_num,:,k,j) = squeeze(awgn(squeeze(channels_ap(ap_num,:,k,j)),15));
        end
    end
end
    
end