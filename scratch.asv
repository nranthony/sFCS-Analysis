%% Sequential Data Split - Stationary Check

split = 5; % aim for 50000 to 100000 photons, or more, per split
split_len = floor(size(macrot,1)/split);
cx = lag_time(23,14);
ACFs = zeros(size(cx,1)-1,split);

split_decay = double(decay) ./ split;
for i = 1:split
    tStart = tic;
    sidx = ((i - 1) * split_len) + 1;
    eidx = (i * split_len);
    [cor, ~] = correctedFCS(split_decay, macrot(sidx:eidx), microt(sidx:eidx), cx);
    cor = cor - 1;
    ACFs(:,i) = cor;
    tEnd = toc(tStart);
    fprintf('ACF(:,%d), photons counts %d -> %d; %.3f seconds.\n', i, sidx, eidx, tEnd);
end

ACFmean = mean(ACFs,2);
ACFstd = std(ACFs,0,2);

figure
hAx=axes;
hAx.XScale = 'log';
xlim([0.07 1e7]);
ylim([-0.02 (mean(ACFmean(ACFmean > 0.1 & ACFmean < 10))*1.85)]);
hold all
errorbar(cxus(1:length(cor)), ACFmean, ACFstd,'x');

figure
hAx = axes;
hAx.XScale = 'log';
xlim([0.07 1e7]);
ylim([-0.02 (mean(ACFmean(ACFmean > 0.1 & ACFmean < 10))*1.85)]);
hold all
cxplot = cxus(1:length(cor));
for i = 1:split
    semilogx(cxplot, ACFs(:,i));
end

%% Find Approx Parameters

x0_norm = [3, 0.2, 1.8, 0.0001, 0.00001];
Gtau = diff3DG(x0_norm, cxus);
 
figure
hAx=axes;
hAx.XScale = 'log';
xlim([0.07 1e7]);
ylim([-0.02 (mean(ACFmean(ACFmean > 0.1 & ACFmean < 10))*1.85)]);
hold all
errorbar(cxus(1:length(cor)), ACFmean, ACFstd,'x');
semilogx(cxus,Gtau);

%% Fit Functions - 3D Gauss, normal diff
% C = x(1);
% wxy = x(2);
% wz = x(3);
% D = x(4);
% Ginf = x(5);
x0_norm = [5, 0.4, 0.7, 0.02, 0.00001]; % see x() parameters in fit functions

lb = [0.001, 0.39999999999999, 0.6999999999999, 0.000001, -0.5];
ub = [100000, 0.4000000000001, 0.7000000000001, 1000000, 0.5];

func = @(x) diff3DG(x, xdata, err);

xdata = cxus; % use microsecond cx array for xdata
ydata = zeros(size(xdata,1),1); % create empty ydata
% fit with diff3DG fit function
x_norm = lsqcurvefit(@diff3DG, x0_norm, ub, lb, xdata, ydata);

figure
semilogx(cxus(1:length(cor)),cor);
% y limits will likely need to be updated
ylim([0.0 1.5]);
xlim([0.07 1e7]);
hold on
semilogx(cxus, Gtau);


%% Fit Functions - 3D Gauss, anomalous diff
% C = x(1);
% wxy = x(2);
% wz = x(3);
% D = x(4);
% Ginf = x(5);
% alpha = x(6);
x0_anom = [100, 0.4, 0.7, 0.02, 0.00001, 0.9]; % see x() parameters in fit functions

xdata = cxus; % use microsecond cx array for xdata
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
[corr_norm, lags] = cross_corr(ps_times, ps_times, 1000000, 10000000000000, 10, 0);

figure
lagsus = lags .* 1E-6; % create microsecond lags
corr_norm = corr_norm - 1; % given we're using intensities/counts, not delta int/counts

semilogx(lagsus, corr_norm);
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
%corr_cells(1,:) = {corr_norm, lags};

parfor i = 1 : 10
    [corr_norm, lags] = cross_corr(ps_times, ps_times, 1000000, 10000000000000, i, 0);
    corr_cells(i,:) = {corr_norm, lags .* 1E-6};
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

