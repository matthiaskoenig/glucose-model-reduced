% functions to get the numerical input profiles for every time point
function [glc] = f_meal(time, tc)
    % time in [s], query in [h]
    time = mod(time, 24);
    glc = interp1(tc.time,tc.glc,time);
end
