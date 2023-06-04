function [channels, rays, is_ray_blocked]  = get_noisy_channels_from_model(model,src,dst,t_res,enable_printing,offset,ant_noise)
    dir_path_len = norm(src-dst);
    d_res = t_res*3e8;
    lambda = model.lambda;
    amps = model.amps;
    if(enable_printing)
%         figure; 
        display_model(model);
    end
    %% Check if there are reflections
    rays{1} = [src;dst]; %Direct path
    % Assuming Direct Path Exists
%     dp_ang = rad2deg(atan2(dst(2)-src(2),dst(1)-src(1)));
%     if dp_ang < 0
%         dp_ang = dp_ang + 360;
%     end
%     src_atten_idx = ceil(src_orient-dp_ang);
%     dst_atten_idx = ceil(dst_orient-dp_ang);
%     if src_atten_idx <= 0
%         src_atten_idx = src_atten_idx + 360;
%     end
%     if dst_atten_idx <= 0
%         dst_atten_idx = dst_atten_idx + 360;
%     end
%     src_atten = db2mag(gain_pattern(src_atten_idx));
%     dst_atten = db2mag(gain_pattern(dst_atten_idx));
    
    ray_status=[0,0];
    for i=1:length(model.reflectors)
        %[p,t]=is_reflect_edit_fast2(src,dst, model.reflectors{i});  % TODO, Add 3D  
        [p,t]=is_reflect_edit_2d(src,dst, model.reflectors{i});  % TODO, Add 3D  
        if(t)
            rays{end+1} =[src;p;dst]; % Add a reflected path
            ray_status(end+1,:)=[1,i];
        end    
    end

    %% Check if a ray is blocked

    % Collect all blocking line segments into one array
    all_blocking_rays ={};
    blocking_ray_status=[];
    for j=1:length(model.reflectors)
        all_blocking_rays{end+1}=model.reflectors{j};
        blocking_ray_status(end+1,:)=[1,j];
    end
    for j=1:length(model.obstacles)
        for k=1:size(model.obstacles{j},1)
            start_point = model.obstacles{j}(k,:);
            if(k==size(model.obstacles{j},1))
                end_point = model.obstacles{j}(1,:);
            else
                end_point = model.obstacles{j}(k+1,:);
            end
            all_blocking_rays{end+1}=[start_point;end_point];
            blocking_ray_status(end+1,:)=[0,j];
        end
    end 

    is_ray_clear = ones(length(rays),1);
    for i=1:length(rays)
        % Check if these line segments are blocking the signal

        for j=1:length(all_blocking_rays)
            if(sum(abs(ray_status(i,:)-blocking_ray_status(j,:)))==0)
                is_ray_clear(i) = model.obs_attenuation*ant_noise;
                if(enable_printing)
                    fprintf('Skipping because a reflector cannot block itself\n');
                end
            else
                for k=1:size(rays{i},1)-1
                    t = is_blocked(rays{i}(k:k+1,:),all_blocking_rays{j}); % Main function
                    if(t && blocking_ray_status(j, 1) == 1)
                        is_ray_clear(i) = model.obs_attenuation*ant_noise;
                    elseif(t && blocking_ray_status(j, 1) == 0) 
                        is_ray_clear(i) = model.ref_attenuation*ant_noise;
                    end
                    if(is_ray_clear(i) ~= 1)
                        break;
                    end
                end
                if(is_ray_clear(i) ~= 1)
                    break;
                end
            end
        end


    end

    if(enable_printing)
        for i=1:length(rays)
            if(is_ray_clear(i))
                clr = 'k';
            else
                clr = 'r';
            end
            for j=1:size(rays{i},1)-1
                plot(rays{i}(j:j+1,1),rays{i}(j:j+1,2),'-','color',clr);
            end
        end
    end

    model.rays = rays;
    % TODO: convert is_ray_blocked to boolean before passing it into model
    is_ray_blocked = is_ray_clear == 0;
    %% Define th e channels
    channels = zeros(length(lambda),1);
    for i=1:length(rays)
        if(is_ray_clear(i) ~= 0)
            cur_d = 0;
            for j=1:size(rays{i},1)-1
                cur_d = cur_d + norm(rays{i}(j+1,[1,2])-rays{i}(j,[1,2]));
            end
            if abs(cur_d-dir_path_len) > d_res % If 
                continue
            end
            for j=1:length(lambda)
%                 if(i==1)
%                     REF_ATTEN = 1;
%                 else
%                     REF_ATTEN = model.ref_attenuation;
%                 end
                channels(j) = channels(j) + is_ray_clear(i)*amps(j)*1/(cur_d/lambda(j))*exp(-1j*2*pi/lambda(j)*(cur_d+offset));% * src_atten*dst_atten;
            end
        end
    end    
end
