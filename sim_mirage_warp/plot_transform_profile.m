function plot_transform_profile(P,theta_vals,d_vals,gt_angle, gt_angle_r, peak_aoa, peak_tof, track_aoa, track_tof)    
    % plot AoA-ToF
    %{
    Plot transform profile for music or fft. Display AoA/AoD-ToF plot.
    
    param:
        P {array, dim length(theta_val) x length(d_val)} -- transform profile.
        theta_val {array, dim 1 x n}                     -- all theta values in radian
        d_val {array, dim 1 x n}                         -- all d values in meters
        gt_angle {double, optional}                      -- ground truth angle in deg
        peak_aoa {array, 1 x n_proposals, optional}      -- profile peak angle(s)
        peak_tof {array, 1 x n_proposals, optional}      -- profile peak tof(s)
        track_aoa {double, optional}                     -- aoa tracking results
        track_tof {double, optional}                     -- tof tracking results
    %}
    assert(isreal(P), "P must be real value array")
    display_legend = false;
    imagesc(d_vals,theta_vals*180/pi,P);
    hold on; 
    
    % plot ground truth line if exist
    if exist('gt_angle', 'var') && ~isempty(gt_angle)
        y = ones(1,length(d_vals))*gt_angle;
        
        y_r = ones(1,length(d_vals))*gt_angle_r;


        plot(d_vals,y,'g','LineWidth',2,'DisplayName','Ground truth-direct path');
        hold on;
        plot(d_vals,y_r,'g:','LineWidth',2,'DisplayName', 'Ground truth-reflected path')
        display_legend = true;
    end
    
    % plot prediction/ peak candidates
    if exist('peak_aoa', 'var') && exist('peak_tof', 'var') && ~isempty(peak_aoa) && ~isempty(peak_tof)
        scatter(peak_tof, peak_aoa, '*r','DisplayName','peaks');
%         display_legend = true;
    end
    
%     % plot tracking results
%     if exist('track_aoa', 'var') && exist('track_tof', 'var') %&& ~isempty(track_aoa) && ~isempty(track_tof)
%         scatter(track_tof, track_aoa, 40, '*g', 'DisplayName', 'peak - tracking');
%         display_legend = true;
%     end
    
    set(gca,'YDir','normal');    % reverse y axis
    xlabel("ToF (m)");
    ylabel("Angle (deg)");
    if display_legend
        legend('Location', 'best');
    end
    hold off;
end

