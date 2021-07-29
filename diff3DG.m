function Gtau = diff3DG(x, xdata, err)
%% 3D Gaussian PSF, single component, normal diffusion
C = x(1);
wxy = x(2);
wz = x(3);
D = x(4);
Ginf = x(5);

Veff = pi.^(1.5) .* wxy .* wxy .* wz; % V effective
ttd = 4 .* D .* xdata ./ (wxy .* wxy); % tau over taud
a = Veff .* C;
b = 1 + ttd;
c = (1 + ttd .* (wxy/wz).^2).^(0.5);

Gtau = Ginf + 1 ./ (a .* b .* c);

Gtau = Gtau + (Gtau .* err);
% Gtau = Gtau .* err;

end

