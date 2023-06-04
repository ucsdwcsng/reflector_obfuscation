function aoas = get_aoa_for_least_tof(profile,d_vals,theta_vals)
    
    profiles = abs(profile);
    [m,n] = size(profiles);
    if(m==length(d_vals))
        profiles = profiles.';
    elseif(n~=length(d_vals))
        error('Check your channel matrix dimensions')
    end
    bw = FastPeakFind(profiles);
    aoa = bw(2:2:end);
    tof = bw(1:2:end);
    if(isempty(tof) || isempty(aoa))
         [~,idx] =  max(abs(profile(:)));
         [aoa,tof] = ind2sub(size(profiles),idx);
    end
        
    if(length(aoa)<2)
        aoas(:) = theta_vals(aoa);
%         tofs(:) = d_vals(tof);
    else
        [~,tof_ids] = sort(tof,'ascend');
        aoas = theta_vals(aoa(tof_ids(1)));
%         tofs(:) = d_vals(tof(tof_ids(1:2)));
    end

end
