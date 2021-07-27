%% Bin pixCounts Matrix

factor = 50; % number off 'bins' to add together
newLen = floor(size(pixCounts,1)/factor);
pixCounts_b = zeros(newLen, size(pixCounts,2));
for i = 1:newLen
    sidx = ((i - 1) * factor) + 1;
    eidx = (i * factor);
%     for j = 1:size(pixCounts,2)
%         pixCounts_b(i,j) = sum(pixCounts(sidx:eidx,j));
%     end
    pixCounts_b(i,:) = sum(pixCounts(sidx:eidx,:));
    %disp(sprintf('%d -> %d', sidx, eidx));
end


%% Create test macrot and markertime data
clc

tst_cyc_dt = 200;
tst_fb_pixels = 1;
tst_nc = 20;
tst_c = 9;

[tst_macrot, tst_mrks, tst_mrk_type, tst_dt] = generate_tst_data(tst_nc, tst_c, tst_cyc_dt, tst_fb_pixels);

%%

[tst_pixCounts, tst_photonCarpet] = sfcsCarpet(uint32(tst_dt), uint32(tst_nc), uint32(tst_c), uint16(tst_macrot), uint64(tst_macrot), uint64(tst_mrks), uint8(tst_mrks));

%%
s=0;
for i = 1:Nc
    s = s + length(photonCarpet{1,i});
end

%% Check sfcs_carpet_.cpp code
mrktimeArr = uint64(zeros(C,1));
pmrkrstime = 1;
for c = 1:C-1
    mrktimeArr(pmrkrstime + c) = marktime(pmrkrstime + (2 * c) + 2) - marktime(pmrkrstime + (2 * c));
end


%% use cross_corr TTTR Matlab Exhange code
[ps_times] = create_abs_ps_times(macrot, microt);
[corr_norm, lags] = cross_corr(ps_times, ps_times, 1000000, 10000000000000, 10, 0);

figure

lagsus = lags .* 1E-6;
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

%% get modal 4-2 and 2-4 times
t42 = [];
for i=2:2:min([2000 length(mktimedelta)])
    t42 = [t42, mktimedelta(i)];
end
mode(t42)

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

