function [ps_times] = create_abs_ps_times(macrot, microt)
    
    macrot = double(macrot);
    microt = double(microt);
    
    n = length(macrot);
    
    ps_times = zeros(n,1);
    ps_macrot = zeros(n,1);
    ps_microt = zeros(n,1);

    
    ps_macrot = macrot .* 25000;
    ps_microt = microt .* 25;
    ps_times = ps_macrot + ps_microt;

end

