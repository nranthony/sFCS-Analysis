%%  Set path

% note '' creates char arrays, while "" creates Matlab strings
%fpath = 'C:\ici-cloud-sections\WBRB Abberior STED\2021\Neil\2021-07-13 - FCS JF cali\TimeHarp_2021-07-13_12-05-48_sFCS60pc_b_fast.ptu';
fpath = 'C:\ici-cloud-sections\WBRB Abberior STED\2021\Neil\2021-07-13 - FCS JF cali\TimeHarp_2021-07-13_11-17-19_pnt30s60pc.ptu';
%fpath = 'Z:\WBRB Abberior STED\2021\ICI\Neil\2021-07-13 - FCS JF cali\TimeHarp_2021-07-13_12-05-48_sFCS60pc_b_fast.ptu';


% convert backslash to forward slash;  other option would be to convert to
% \\ to avoid \ being an escape character
fpath = strrep(fpath,'\','/');

%% compile if needed - mexw64 file should work as is on most Windows machines
% needed for below ptu_data function call
% requires compiler; recommend VS2019
% uncomment below as needed and run
% clc
% mex ptu_data.cpp
% mex -g ptu_data.cpp % for debugging in VS2019

%%  Load Data - Requires files in https://github.com/nranthony/TTTR/tree/master/PTU/MatLab
% not needed if using other method to access count data

% pull data from ptu file
% macrot: counts in units of pulse/sync
% microt: tcspc bin in units of time resolution
% decay: histogram of all counts
% marktime: marker time in units of pulse/sync
% marktype: marker time at above time - NOTE please double check,
% especially after hardware or service changes
% - Abberior FL PQ TimeHarp 260 P
% - Marker 1: Frame start
% - Marker 4: Line start
% - Marker 2: Line end
% meta: cell array of related metadata from ptu header
% fpath: full filepath as char, not Matlab string
% boolean flag: 1,0 for display header in command window
[macrot, microt, decay, marktime, marktype, meta] = ptu_data(fpath, 1);

%% Remove Points Prior to Frame Marker and Set First Point to Time Zero

[row,] = find(macrot > marktime(1)); % assumes first marker is frame marker
% pre measurement noise removed from macro and micro - must match perfectly
% Note: do this only once
macrot = macrot(row(1):end);
microt = microt(row(1):end);

macrot = macrot - marktime(1); % this uses marker 1 as the t zero
marktime = marktime - marktime(1); % also remove from marker time

%% compile if needed - mexw64 file should work as is on most Windows machines
% needed in correctedFCS.m
% requires compiler; recommend VS2019
% uncomment below as needed and run
% clc
% mex calcSymmetricACF.c
% mex -g calcSymmetricACF.c % for debugging in VS2019

%% Time Symmetric ACF

% generate log scaled tau
cx = lag_time(23,14);
% optional: view lag times
%semilogy(cx);

% note correctedFCS uses parfor, and needs parallel processing kit
% change to for if not installed
% first run will take longer, but subsequent runs will be significantly
% quicker depending on core count
tic;
[cor, pwave] = correctedFCS(decay, macrot, microt, cx);
cor = cor - 1;
toc;

%% Prep for sFCS/lineFCS Carpet

% sFCS or lineFCS
sFCS = true;

% create delta marker time
sizemrkr = length(marktime);
mktimedelta = zeros(sizemrkr,1);
for i=2:sizemrkr
    mktimedelta(i) = marktime(i) - marktime(i-1);
end
clear sizemrkr;

% get number of pixels along/around each cycle: Nc in C++ code
% use PixX + 1 for sFCS and pixX for lineFCS
PixX = meta{1,3};
% get time between cycles
% use line start, mrk#4, to line end, mrk#2, for lineFCS
% use line start, to line start for sFCS, mrk#4->mrk#4; here the line end
% to line start gap #2->#4 is 1 x pixel, hence above pixX + 1
% assumes 1,4,2,4,2... markers, i.e. single frame marker at the start
tarr = [];
if (sFCS)
    for i=2:2:min([4000 length(mktimedelta)])
        tarr = [tarr, (mktimedelta(i)+mktimedelta(i+1))];
    end
    Nc = PixX + 1;
else % lineFCS
    for i=3:2:min([4000 length(mktimedelta)])
        tarr = [tarr, mktimedelta(i)];
    end
    Nc = PixX;
end
cycleTime = mode(tarr);
deltat = round(cycleTime/Nc);
clear PixX;
clear tarr;
C = sum(marktype == 4) - 1;

%% compile if needed - mexw64 file should work as is on most Windows machines
% needed for below sfcsCarpet function call
% requires compiler; recommend VS2019
% uncomment below as needed and run
% clc
% mex sfcsCarpet.cpp 
% mex -g sfcsCarpet.cpp % for debugging in VS2019

%% Generate Carpet

[pixCounts, photonCarpet] = sfcsCarpet(uint32(deltat), uint32(Nc), uint32(C), uint16(microt), uint64(macrot), uint64(marktime), uint8(marktype));



%% plot results - change units and edit limits for your needs

% convert number of pulses to time scaled abcissa (here 25ns pulse to pulse
% window
cxus = double(cx) .* 0.025;  % converts to µs
semilogx(cxus(1:length(cor)),cor);
% y limits will likely need to be updated
ylim([0.0 0.9]);
xlim([0.07 1e7]);







