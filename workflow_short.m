%%  Import TTTR Data

% note '' creates char arrays, while "" creates Matlab strings
% for below map ICI cloud and replace C:\ici-cloud-sections with Z:
% see https://www.cores.emory.edu/ici/resources/ici-cloud.html

%% note to self - change permissions on below folders

%% Point Example
fpath = 'C:\ici-cloud-sections\WBRB Abberior STED\2021\Neil\2021-07-13 - FCS JF cali\TimeHarp_2021-07-13_11-17-19_pnt30s60pc.ptu';
% use for point measurements
[macrot, microt, decay, meta] = import_point_ptu(fpath);


%% Scanning Example
fpath = 'C:\ici-cloud-sections\WBRB Abberior STED\2021\Neil\2021-07-13 - FCS JF cali\TimeHarp_2021-07-13_12-05-48_sFCS60pc_b_fast.ptu';
% use for scanning measurements
sfcs = 1;
[macrot, microt, decay, marktime, marktype, meta, deltat, Nc, C] = import_sfcs_ptu(fpath, sfcs);



%% compile if needed - mexw64 file should work as is on most Windows machines
% needed in correctedFCS.m
% requires compiler; recommend VS2019
% uncomment below as needed and run
% clc
% mex calcSymmetricACF.c
% mex -g calcSymmetricACF.c % for debugging in VS2019

%% Time Symmetric ACF

% generate log scaled tau, see lag_time.m
% shorten tau to remove lag times < ###; here < 300ns
tau = lag_time(23,14);
tauus = double(tau) .* 0.025;  % converts to µs
shortest_lag = 0.3; % 300 ns
tau = tau(tauus >= shortest_lag);
tauus = tauus(tauus >= shortest_lag);
% optional: view lag times
%semilogy(tau);

% note correctedFCS uses parfor, and needs parallel processing kit
% change to 'for' if not installed
% first run will take longer, and display error (sometimes)
% subsequent runs will be significantly quicker depending on core count, and
% not error
tic;
[cor, pwave] = correctedFCS(decay, macrot, microt, tau);
toc;


%% compile if needed - mexw64 file should work as is on most Windows machines
% needed for below sfcsCarpet function call
% requires compiler; recommend VS2019
% uncomment below as needed and run
% clc
% mex sfcsCarpet.cpp 
% mex -g sfcsCarpet.cpp % for debugging in VS2019

%% Generate Carpet
% __Inputs:
% deltat: width of each pixel in laser pulses -> cycleTime / Nc  -> different with circular vs line scans
% Nc: number of points in each 'cycle'
% C: number of cycles, i.e. how many complete revolutions in entire dataset
% microt: ps collection window of photon between laser pulses
% macrot: number of laser pulses from beginning of dataset - note, time zero is modified to be the start of the first line/circle
% marktime: number of pulses count for the frame and start/stop markers
% marktype: index of marker type: 1-frame, 4-line start, 2-line end
% __Outputs:
% pixCounts: 2D array of photon counts binned into scan pixels
% photonCarpet: cell array of corresponding pairs of microt and macrot arrays for each pixel
[pixCounts, photonCarpet] = sfcsCarpet(uint32(deltat), uint32(Nc), uint32(C), uint16(microt), uint64(macrot), uint64(marktime), uint8(marktype));


%% plot results - change units and edit limits for your needs

% convert number of pulses to time scaled abcissa (here 25ns pulse to pulse
% window
tauus = double(tau) .* 0.025;  % converts to µs
semilogx(tauus(1:length(cor)),cor);
% y limits will likely need to be updated
ylim([0.0 1.5]);
xlim([0.07 1e7]);







