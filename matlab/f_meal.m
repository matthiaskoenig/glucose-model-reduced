% Get the interpolated timecourse profile for the given timepoints time
% from given timecourse information.
% Uses linear interpolation with timepoints for queries in [h].
% Timecourses are periodic with T=24h.
function [glc] = f_meal(timepoints, tc)
    % handle periodicity
    timepoints = mod(timepoints, 24);
    glc = interp1(tc.time, tc.glc, timepoints, 'linear');
end
