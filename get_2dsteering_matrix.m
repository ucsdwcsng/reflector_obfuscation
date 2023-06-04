function S = get_2dsteering_matrix(theta_vals,d_vals,n_sub,f,df,ant_sep)    
    Phi = zeros(3,length(theta_vals));
    for i = 1:length(theta_vals)
        phi_theta = exp(-1j*2*pi*f/3e8*sin(theta_vals(i))*ant_sep);
        Phi(1,i) = 1;
        Phi(2,i) = phi_theta;
        Phi(3,i) = phi_theta.^2;
    end
    Omega = zeros(n_sub/2,length(d_vals));
    for i = 1: length(d_vals)
        omega_t = exp(-1j*2*pi*df*d_vals(i)/3e8);
        Omega(1:n_sub/2,i) = omega_t.^((1:n_sub/2)');
    end
    
    S = kron(Phi,Omega);
end