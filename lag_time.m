function [tau] = lag_time(n_casc, B)
    %LAG_TIME creates an array of lag times tau
    % for use with symmetric ACF
    jmax = n_casc * B;
    tau = ones(jmax,1);
    
    for j = 2 : jmax
        exponent = floor( (j-1) / B );
        tau(j) = tau(j-1) + 2.^exponent;
    end
    
    tau = uint32(tau(1:length(tau)-1));
    
end