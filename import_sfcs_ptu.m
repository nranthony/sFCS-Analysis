function [macrot, microt, decay, marktime, marktype, meta, deltat, Nc, C] = import_sfcs_ptu(fpath, sfcs)
% fpath: full path with standard Windows backslash '\' and using '', not ""
% sfcs: boolean for using sfcs, 1, or line fcs, 0
% sfcs will use the one pixel between lineend and linestart as an
% additional pixel in the circular scan
% line will only use the pixels in the line between linestart and lineend


    %%  Set path
    % note '' creates char arrays, while "" creates Matlab strings
    % example paths:
    %fpath = 'C:\ici-cloud-sections\WBRB Abberior STED\2021\Neil\2021-07-13 - FCS JF cali\TimeHarp_2021-07-13_12-05-48_sFCS60pc_b_fast.ptu';
    %fpath = 'Z:\WBRB Abberior STED\2021\ICI\Neil\2021-07-13 - FCS JF cali\TimeHarp_2021-07-13_12-05-48_sFCS60pc_b_fast.ptu';

    % convert backslash to forward slash;  other option would be to convert to
    % \\ to avoid \ being an escape character
    fpath = strrep(fpath,'\','/');

    %% compile if needed - mexw64 file should work as is on most Windows machines
    % requires compiler; recommend VS2019
    % uncomment below two lines and run
    % clc
    % mex ptu_data.cpp

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


    %% Prep for sFCS/lineFCS Carpet

    % sFCS or lineFCS

    % create delta marker time
    sizemrkr = length(marktime);
    mktimedelta = zeros(sizemrkr,1);
    for i=2:sizemrkr
        mktimedelta(i) = marktime(i) - marktime(i-1);
    end

    % get number of pixels along/around each cycle: Nc in C++ code
    % use PixX + 1 for sFCS and pixX for lineFCS
    PixX = meta{1,3};
    % get time between cycles
    % use line start, mrk#4, to line end, mrk#2, for lineFCS
    % use line start, to line start for sFCS, mrk#4->mrk#4; here the line end
    % to line start gap #2->#4 is 1 x pixel, hence above pixX + 1
    % assumes 1,4,2,4,2... markers, i.e. single frame marker at the start
    tarr = [];
    if (sfcs)
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

end

