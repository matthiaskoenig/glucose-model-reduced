function [] = timecourse_sim()
% Plotting the timecourse information
clc, format compact, clear all, close all
addpath('../glucose-profiles/')

disp('*** Load data ***')
tc1 = timecourse_1meal()
tc3 = timecourse_3meals()

disp('*** Plot data ***')
ttest = linspace(0,48,200);
fig1 = figure('Name' , 'Timecourse glucose', 'Color', [1 1 1]);
subplot(1,2,1)
plot(ttest, f_meal(ttest, tc1), 'bo-'), hold on
plot(tc1.time, tc1.glc, 'ko-'), hold off
ylabel('glucose [mmol/L]')
xlabel('time [h]')
title('1 Meal')
axis square

subplot(1,2,2)
plot(ttest, f_meal(ttest, tc3), 'bo-'), hold on
plot(tc3.time, tc3.glc, 'ko-'), hold off
xlabel('time [h]')
ylabel('glucose [mmol/L]')
title('3 Meals')
axis square

% Perform the timecourse simulation for the 1 Meal and 3 Meal simulation
% TODO:
% perform simulations with the data as input data
% how is the data used as input?

end
