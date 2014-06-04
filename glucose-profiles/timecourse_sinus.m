function [tc] = timecourse_sinus()
T = 24; % [h] period

tc.time = 0:0.5:24; % [h]
tc.glc = 5.5 - 2.0*cos(2*pi/T *tc.time); % [mM] [3.5-7.5]
tc.lac = 4.0 + 2.0*sin(2*pi/T *tc.time); % [mM] [2.0-6.0]
end
