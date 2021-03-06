%% Sequential Data Split - Stationary Check - Error Estimation
% required for fitting - needs to be consistent across analyses

% shorten tau to remove lag times < ###
tau = lag_time(23,14);
tauus = double(tau) .* 0.025;  % converts to µs
shortest_lag = 0.3; % 300 ns
tau = tau(tauus >= shortest_lag);
tauus = tauus(tauus >= shortest_lag);

% setup sequential spliting of photon counts
split = 5; % aim for 50000 to 100000 photons, or more, per split
split_len = floor(size(macrot,1)/split);
ACFs = zeros(size(tau,1)-1,split);

% decay histogram scaled; used to normalize in symmetricACF calc
split_decay = double(decay) ./ split;
for i = 1:split
    tStart = tic;
    sidx = ((i - 1) * split_len) + 1;
    eidx = (i * split_len);
    [cor, ~] = correctedFCS(split_decay, macrot(sidx:eidx), microt(sidx:eidx), tau);
    ACFs(:,i) = cor;
    tEnd = toc(tStart);
    fprintf('ACF(:,%d), photons counts %d -> %d; %.3f seconds.\n', i, sidx, eidx, tEnd);
end

clear tStart tEnd sidx eidx shortest_lag i;

% generate mean and std across all splits for each tau value
ACFmean = mean(ACFs,2);
ACFstd = std(ACFs,0,2);
ACFvar = ACFstd.^2;

% generate ACF for all data
[ACFall, ~] = correctedFCS(decay, macrot, microt, tau);
% auto scale y axis testing
ACFall_var = movvar(ACFall, [0 6]); % moving variance 0 backwards, 6 + 1 forward
ylim_upper = mean(ACFall(ACFall_var < 0.01 & ACFall_var > 0.0007)) .* 1.85;

%% Fiting Using Standard fit(...) Command


tau = tauus(1:size(ACFs,1));

% point FCS
ft = fittype('diff3DG(tau, C, wxy, wz, D, Ginf)', ...
    'coefficients', {'C', 'wxy', 'wz', 'D', 'Ginf'}, ...
    'dependent', {'Gtau'}, ...
    'independent', {'tau'});

wxy = 0.2;
wz = 1.8;
coef_0 = [3, wxy, wz, 0.002, 0.00001]; % see x() parameters in fit functions
lb = [0.001, wxy, wz, 0.000001, 0];
ub = [1000, wxy, wz, 10, 0.05];
fitOptions = fitoptions( ...
    'Method', 'NonlinearLeastSquares', ...
    'Lower', lb, 'Upper', ub, ...
    'Display', 'iter', ...
    'Weights', ACFstd);

f = fit( tau, ACFmean, ft, 'StartPoint', coef_0);

% scanning FCS circular
ft = fittype('diff3DG_SFCS(tau, C, wxy, wz, D, Ginf, R, f)', ...
    'coefficients', {'C', 'wxy', 'wz', 'D', 'Ginf', 'R', 'f'}, ...
    'dependent', {'Gtau'}, ...
    'independent', {'tau'});

wxy = 0.2;
wz = 1.8;
R = 0.5; % ExpControl X 500nm in Abberior htm meta export
f = 1 / (double(Nc) .* meta{1,5} .* 1000); % frequency
coef_0 = [6, wxy, wz, 0.002, 0.00001, R, f]; % see x() parameters in fit functions
lb = [0.001, wxy, wz, 0.000001, 0, R, f];
ub = [1000, wxy, wz, 10, 0.05, R, f];
fitOptions = fitoptions( ...
    'Method', 'NonlinearLeastSquares', ...
    'Lower', lb, 'Upper', ub, ...
    'Display', 'iter', ...
    'Weights', ACFstd);

f = fit( tau, ACFmean, ft, 'StartPoint', coef_0);

Gtau = diff3DG_SFCS(tau, 3, wxy, wz, 0.002, 0, R, f);


figure
hAx = axes;
hAx.XScale = 'log';
xlim([0.07 1e7]);
ylim([-0.01 ylim_upper]);
hold all
plot( tau, Gtau );

plot( f, tau, ACFmean );


%% SIDE NOTE __Error estimation__:
%
% Can use moving var:
% e.g. we'd add a fraction of the array to itself in the fit function to
% weight the loss/error during fitting; gtau = gtau + (gtau .* ACFall_var)
%
% Can generate error from difference between smoothed ACFall and ACFmean:
% e.g. simply add the difference between the SG smoothed ACFall and the
% ACFmean; gtau = gtau + (ACFmean - ACFsmthSG)
% TODO - check smoothing window size for sFCS (not carpet) data
%
% Use square of ACFstd array; ACFstd.^2

%  NOTE  --  Needs to be included in the fit function; therefore, needs to
%  be positive and negative around the 'mean'
%  - look to delta of smoothed ACFall and mean of sequential splits
%  - err = ACFmean - ACFsmthSG;

% Follow Up - Smoothing too smooth for sFCS decay

%% Plot All Splits and Mean +- Std - requires above variables

% IMPORTANT!
% plot all sequentially split ACFs to look for significant y offset
% y offset between ACFs of different time windows can be due to aggregates
% - note look for offsets at longer lag times, but ignore everything past a
% lag time of ~ 1 second, 10^6 us
figure
hAx = axes;
hAx.XScale = 'log';
xlim([0.07 1e7]);
ylim([-0.01 ylim_upper]);
hold all
tauplot = tauus(1:size(ACFs,1));
for i = 1:split
    semilogx(tauplot, ACFs(:,i));
end

% plots the mean with std error bars
figure
hAx=axes;
hAx.XScale = 'log';
xlim([0.07 1e7]);
ylim([-0.01 ylim_upper]);
hold all
tauplot = tauus(1:length(ACFmean));
errorbar(tauplot, ACFmean, ACFstd,'x');


%% Use Smoothed ACFall - requires above variables

% set smoothing window size
smth_window = 45;
% smoothing types testing; I think SG will be best for changes seen in sfcs
ACFsmthG = smoothdata(ACFall, 'gaussian', smth_window);
ACFsmthM = smoothdata(ACFall, 'movmedian', smth_window);
ACFsmthSG = smoothdata(ACFall, 'sgolay', smth_window);


figure
hAx=axes;
hAx.XScale = 'log';
xlim([0.07 1e7]);
ylim([-0.01 ylim_upper]);
hold all
tauplot = tauus(1:length(ACFall));
semilogx(tauplot, ACFall);
semilogx(tauplot, ACFsmthG);
semilogx(tauplot, ACFsmthM);
semilogx(tauplot, ACFsmthSG);
legend('ACF', 'Gauss', 'Median', 'Savitzky-Golay');



%% Find Approx Fit Parameters

x0_nDif = [3, 0.2, 1.8, 0.0001, 0.00001];
Gtau = diff3DG(x0_nDif, tauplot, zeros(length(tauplot),1));
 
figure
hAx=axes;
hAx.XScale = 'log';
xlim([0.07 1e7]);
ylim([-0.02 ylim_upper]);
hold all
errorbar(tauplot, ACFmean, ACFstd,'x');
semilogx(tauplot,Gtau);

%% Fit Functions - 3D Gauss, normal diff
% C = x(1);
% wxy = x(2);
% wz = x(3);
% D = x(4);
% Ginf = x(5);

% wxy and wz are the 1/e point of the assumed 3D Gaussian volume
wxy = 0.2;
wz = 1.8;
x0_nDif = [3, wxy, wz, 0.0002, 0.00001]; % see x() parameters in fit functions

lb = [0.001, wxy, wz, 0.000001, -0.05];
ub = [100000, wxy, wz, 1000000, 0.5];

xdata = tauplot; % use microsecond tau array for xdata
ydata = zeros(size(xdata,1),1); % create empty ydata

% err = ACFvar; % does not work as all positive
err = ACFmean - ACFsmthSG; % change to be 'normalized' to some degree; now multiplying in diff3DG
errAbs = abs(err);
err = err ./ max(errAbs);

% function reference below 'fills' the call with current xdata and err;
% call function reference again if/when arrays change
func = @(x) diff3DG(x, xdata, err);

% set options
options.Algorithm = 'levenberg-marquardt';
% options.Algorithm = 'trust-region-reflective';
options.Display = 'iter'; % comment out if annoying
% fit with diff3DG fit function
% [x,resnorm,residual,exitflag,output] = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
%[x_nDif,resnorm,residual,exitflag,output] = lsqcurvefit(func,x0_nDif,xdata,ydata,lb,ub,options);

% [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,options)
[x_nDif,resnorm,residual,exitflag,output] = lsqnonlin(func,x0_nDif,lb,ub,options);


Gtau = diff3DG(x_nDif, xdata, zeros(length(tauplot),1));

% Gtau = diff3DG(x_nDif, xdata, err);


tiledlayout(4,1);
hAx = nexttile(1,[3 1]);
hAx.XScale = 'log';
xlim([0.07 1e7]);
ylim([-0.02 ylim_upper]);
hold on
errorbar(tauplot, ACFmean, ACFstd,'x');
semilogx(tauplot, Gtau);


%% Fit Functions - 3D Gauss, anomalous diff
% C = x(1);
% wxy = x(2);
% wz = x(3);
% D = x(4);
% Ginf = x(5);
% alpha = x(6);
x0_anom = [100, 0.4, 0.7, 0.02, 0.00001, 0.9]; % see x() parameters in fit functions

xdata = tauus; % use microsecond tau array for xdata
ydata = zeros(size(xdata,1),1); % create empty ydata
% fit with diff3DG fit function
x_anom = lsqcurvefit(@anomalousDiff3DG, x0_anom, xdata, ydata);

%% Bin pixCounts Matrix

factor = 50; % number off 'bins' to add together
newLen = floor(size(pixCounts,1)/factor);
pixCounts_b = zeros(newLen, size(pixCounts,2));
for i = 1:newLen
    sidx = ((i - 1) * factor) + 1;
    eidx = (i * factor);
    pixCounts_b(i,:) = sum(pixCounts(sidx:eidx,:));
end

%% Create test macrot and markertime data
clc
tst_cyc_dt = 200;
tst_fb_pixels = 1;
tst_nc = 20;
tst_c = 9;

[tst_macrot, tst_mrks, tst_mrk_type, tst_dt] = generate_tst_data(tst_nc, tst_c, tst_cyc_dt, tst_fb_pixels);

%% Run sfcsCarpet on Test Data
[tst_pixCounts, tst_photonCarpet] = sfcsCarpet(uint32(tst_dt), uint32(tst_nc), uint32(tst_c), uint16(tst_macrot), uint64(tst_macrot), uint64(tst_mrks), uint8(tst_mrks));

%%  Check Total Counts in Photon Count Carpets
s=0;
for i = 1:Nc
    s = s + length(photonCarpet{1,i});
end

%% use cross_corr TTTR Matlab Exhange code
[ps_times] = create_abs_ps_times(macrot, microt);
[corr_nDif, lags] = cross_corr(ps_times, ps_times, 1000000, 10000000000000, 10, 0);

figure
lagsus = lags .* 1E-6; % create microsecond lags
corr_nDif = corr_nDif - 1; % given we're using intensities/counts, not delta int/counts

semilogx(lagsus, corr_nDif);
ylim([1 1.1])
xlim([0.1 1e6])

%% plot all corr funcs in correlation matrix created from cross_corr_subsets.m

figure;
set(gca, 'XScale', 'log')
hold on;
for i = 1:size(A,2)
    semilogx(lags,A(:,i));
end

ylim([0.75 1.25]);
xlim([0.1 1e6]);

%% create delta marker time
sizemrkr = length(marktime);
mktimedelta = zeros(sizemrkr,1);
for i=2:sizemrkr
    mktimedelta(i) = marktime(i) - marktime(i-1);
end

%% loop through resolution options in cross_corr

corr_cells = cell(10,2);
%corr_cells(1,:) = {corr_nDif, lags};

parfor i = 1 : 10
    [corr_nDif, lags] = cross_corr(ps_times, ps_times, 1000000, 10000000000000, i, 0);
    corr_cells(i,:) = {corr_nDif, lags .* 1E-6};
end

%% plot different resolution options
figure
semilogx(corr_cells{1,2}, corr_cells{1,1});
hold on

for i = 2 : 10
    semilogx(corr_cells{i,2}, corr_cells{i,1});
end

ylim([1 1.1])
xlim([0.1 1e6])

hold off

