% Get values of the interpolated timecourse profile 
% for the given timepoints. 
% Points for interpolation are given as struct tc.
% Uses linear interpolation with timepoints for queries in [h].
% Timecourses are periodic with T=24h.
%
% Matthias Koenig
% 2014-06-04

function [glc, lac] = f_timecourse(timepoints, tc)
    % handle periodicity
    timepoints = mod(timepoints, 24);
    glc = interp1(tc.time, tc.glc, timepoints, 'linear');
    lac = interp1(tc.time, tc.lac, timepoints, 'linear');
end