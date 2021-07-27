function [tst_macrot, tst_mrks, tst_mrk_type, tst_dt] = generate_tst_data(tst_nc, tst_c, tst_cyc_dt, tst_fb_pixels)

tst_mrks = [];
tst_mrk_type = [];
tst_macrot = [];

% test data needs to have a deltat that's divisible by Nc without remainder 
if ( rem( tst_cyc_dt, tst_nc ) )
    disp("Delta T (tst_dt) must be exactly divisible by Nc (tst_nc).");
    return 
end

tst_dt = tst_cyc_dt / tst_nc; % one complete cycle, dt, made up of nc pixels
tst_fb_dt = tst_fb_pixels * tst_dt; % set 'flyback' time
tst_line_dt = tst_cyc_dt - tst_fb_dt; % set line time

% use line start, mrk#4, to line end, mrk#2, for lineFCS
% use line start, to line start for sFCS, mrk#4->mrk#4; here the line end
% to line start gap #2->#4 is 1 x pixel, hence above pixX + 1
% assumes 1,4,2,4,2... markers, i.e. single frame marker at the start
mrk_type = [4,2];

tst_mrk_cnt = 1;
tst_mrks(tst_mrk_cnt) = 0;
tst_mrk_type(tst_mrk_cnt) = 1;



for c = 2:(tst_c + 2)
    % marker #4 line start
    tst_mrk_cnt = tst_mrk_cnt + 1;
    tst_mrks(tst_mrk_cnt) = tst_mrks(tst_mrk_cnt - 1) + tst_fb_dt;
    tst_mrk_type(tst_mrk_cnt) = mrk_type( 1 + rem(tst_mrk_cnt,2) );
    
    % marker #2 line end; except for last
    if (c ~= (tst_c + 2))
        tst_mrk_cnt = tst_mrk_cnt + 1;
        tst_mrks(tst_mrk_cnt) = tst_mrks(tst_mrk_cnt - 1) + tst_line_dt;
        tst_mrk_type(tst_mrk_cnt) = mrk_type( 1 + rem(tst_mrk_cnt,2) );
    end
end

tst_macro_cnt = 1;
for c = 1:tst_c
    tmp_mrk_time = tst_mrks(c*2);
    for n = 1:tst_nc
        for m = 1:3
            tst_macrot(tst_macro_cnt) = tmp_mrk_time ...
                + ((n - 1) * tst_dt) + m;
            tst_macro_cnt = tst_macro_cnt + 1;
        end
        
    end
end

%[tst_pixCounts, tst_photonCarpet] = sfcs_carpet_(uint32(tst_dt), uint32(tst_nc), uint32(tst_c), uint16(tst_macrot), uint64(tst_macrot), uint64(tst_mrks), uint8(tst_mrks));

tst_macrot = tst_macrot.';
tst_mrks = tst_mrks.';
tst_mrk_type = tst_mrk_type.';


end

