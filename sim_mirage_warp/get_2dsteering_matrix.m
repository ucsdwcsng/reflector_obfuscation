function S = get_2dsteering_matrix(theta_vals,d_vals)
    %{
    Defines the Steering Matrices given the required search space in theta_val
    and d_vals. Refers to spotfi https://web.stanford.edu/~skatti/pubs/sigcomm15-spotfi.pdf
    for details about formulation.
    
    Args:
        theta_vals {vector, dim 1 x n} -- all possible AoA/AoD values.
        d_vals {vector, dim 1 x n}     -- search space for all possible ToF
        opt {struct}                   -- contains all constants.
    
    Return:
        S {array} -- dim (n_antenna * n_sub/2) x (len(theta_vals) * len(d_vals))
        steering matrix for all possible theta and d in search space.
    %}
    
    % get all constants from opt
    n_sub = 48; %length(opt.subcarrier_indices); % number of subcarriers.
    f = 5.18e9;                       % OFDM central frequency


BW = 20e6;
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
subcarrier_indices = SC_IND_DATA; %48 subcarrier
% CHANNEL = 1;
% freq =f + subcarrier_indices.*BW./n_sub;

%     df =312500;% median(diff(freq));            % subcarrier spacing.
%     df = 300000;
     freq = fftshift(5.18e9 + subcarrier_indices.*BW./48);

    df = median(diff(freq));
    ant_sep =0.0259;                  % antenna separation.
    
    % steering vector -> AoA, spotfi Eq. 2
    Phi = zeros(2,length(theta_vals));
    for i = 1:length(theta_vals)
        phi_theta = exp(-1j*2*pi*f/3e8*sin(theta_vals(i))*ant_sep);
        Phi(1,i) = 1;
        Phi(2,i) = phi_theta;
        Phi(3,i) = phi_theta.^2; % assumes there are three antennas
    end
    
    % Omega -> ToF, spotFi Eq. 6
    Omega = zeros(n_sub/2,length(d_vals));
    for i = 1: length(d_vals)
        omega_t = exp(-1j*2*pi*df*d_vals(i)/3e8);
        Omega(:,i) = omega_t.^((1:n_sub/2)');
    end
    
    S = kron(Phi,Omega);
end