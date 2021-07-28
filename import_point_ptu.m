function [macrot, microt, decay] = import_point_ptu(fpath)
% fpath: full path with standard Windows backslash '\' and using '', not ""

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
    [macrot, microt, decay, marktime, marktype, meta] = ptu_data(fpath, 0);

    %% Remove Points Prior to Frame Marker and Set First Point to Time Zero

    [row,] = find(macrot > marktime(1)); % assumes first marker is frame marker
    % pre measurement noise removed from macro and micro - must match perfectly
    % Note: do this only once
    macrot = macrot(row(1):end);
    microt = microt(row(1):end);

    macrot = macrot - marktime(1); % this uses marker 1 as the t zero
    marktime = marktime - marktime(1); % also remove from marker time

end

