function [channels_w_offset,channels_wo_offset,offset] = ...
    call_get_channels_from_model(offset_model,no_offset_model,pos,ap, offset)

n_ant_per_ap = size(ap{1},1);

channels_wo_offset = zeros(length(no_offset_model.lambda),length(ap),n_ant_per_ap);
channels_w_offset = zeros(length(offset_model.lambda),length(ap),n_ant_per_ap);

for j=1:length(ap)
    
    for k=1:n_ant_per_ap
        
        channels_w_offset(:,j,k)=get_channels_from_model_edit(no_offset_model, ...
                                            pos, ap{j}(k,:), false, offset(j));
        channels_wo_offset(:,j,k)=get_channels_from_model_edit(offset_model, ...
                                            pos, ap{j}(k,:), false, 0);

        mean_power_w_offset = db(mean(abs(squeeze(channels_w_offset(:,j,k)))))./2;
        mean_power_wo_offset = db(mean(abs(squeeze(channels_wo_offset(:,j,k)))))./2;
    
        channels_w_offset(:,j,k) = awgn(squeeze(channels_w_offset(:,j,k)), ...
                                        20, mean_power_w_offset, 'db');
        channels_wo_offset(:,j,k) = awgn(squeeze(channels_wo_offset(:,j,k)), ...
                                        20, mean_power_wo_offset, 'db');
    end

end