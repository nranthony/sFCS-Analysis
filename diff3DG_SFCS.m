function Gtau = diff3DG_SFCS(tau, C, wxy, wz, D, Ginf, R, f) % for fit fittype
%% 3D Gaussian PSF, single component, normal diffusion, circular scan

Veff = pi.^(1.5) .* wxy .* wxy .* wz; % V effective
ttd = 4 .* D .* tau ./ (wxy .* wxy); % tau over taud
a = Veff .* C;
b = 1 + ttd;
c = (1 + ttd .* (wxy/wz).^2).^(0.5);

Gtau = Ginf + 1 ./ (a .* b .* c);

expon = ( R.^2 .* (sind(pi .* f .* tau)).^2 ) ./ ( (wz/wxy).^2 + D .* tau );

Gtau = Gtau .* exp(-expon);

end

