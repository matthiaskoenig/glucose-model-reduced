% Plot the various timecourses ('stationary', '1meal', '3meals')
%
% Matthias Koenig
% 2014-06-04
clear all
addpath('../glucose-profiles/')

disp('*** Load data ***')
tcs = {timecourse_1meal(), timecourse_3meals(), timecourse_sinus()};
names = {'1 Meal', '3 Meals', 'Sinus'}

disp('*** Plot data ***')
time = linspace(0,48, 5000);
fig1 = figure('Name' , 'Timecourses', 'Color', [1 1 1]);

% glc
for k=1:numel(tcs)
    subplot(2,numel(tcs),k)
    tc = tcs{k};
    [glc, lac] = f_timecourse(time, tc);
    plot(time, glc, 'bo-'), hold on
    plot(tc.time, tc.glc, 'ko-'), hold off
    ylabel('glucose [mmol/L]')
    xlabel('time [h]')
    title(names{k})
    axis square
end
% lac
for k=1:numel(tcs)
    subplot(2,numel(tcs),numel(tcs)+k)
    tc = tcs{k};
    [glc, lac] = f_timecourse(time, tc);
    plot(time, lac, 'bo-'), hold on
    plot(tc.time, tc.lac, 'ko-'), hold off
    ylabel('lactate [mmol/L]')
    xlabel('time [h]')
    title(names{k})
    axis square
end

% save the figure
print(fig1, '../results/timecourse_profiles.png', '-dpng')