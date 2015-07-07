function [] = fit_glycogen()
% Fitting the glycogen timecourse for the 1meal with altered model
% structure
clear all
%close all

% load the glycogen data points for fitting
[h, m] = csvread_with_headers('MV12_Matlab_ODE-1meal.csv');
ind_time = find(ismember(h,'time'));
ind_gly = find(ismember(h,'gly'));
time = m(:, ind_time);
gly = m(:, ind_gly);

%% glc changes
glc = linspace(4,8,40)
figure()
plot(glc, glc, '-o')
close all
figure()

glc_new =  glc + 2.5*(glc - 6.2)/6.2
glc_new2 =  glc + 5.0*(glc - 6.2)/6.2
plot(glc, glc, 'k-'), hold on
plot(glc, glc_new, 'r-o'), hold on
plot(glc, glc_new2, 'b-o'), hold on
%plot(glc, glc.^2 - 5.5*glc+5.5, 'b-o'), hold on
grid on



% define the optimization function
% fit parameters:
% solve the system with the given fit parameters
[t_data, f_data, c_data] = mv_solve_flow(2);

% calculate the glycogen timecourse
gly_sim = c_data(5,:)';

% calculate value to minimize
F = sum( (gly-gly_sim).^2 )

fig1 = figure()
plot(time, gly, 'k-', time, gly_sim, 'b-')

pars0 = [0.38] % Starting guess
lb = [0.3]     % lower bounds
ub = [0.7]     % upper bounds
% [pars, resnorm] = lsqnonlin(@myfun, pars0, lb, ub); % Call optimizer


function F = myfun(pars)
    alpha = pars(1)
    % solve the system with the given fit parameters
    [t_data, f_data, c_data] = mv_solve_flow(2);

    % calculate the glycogen timecourse
    gly_sim = c_data(5,:)';

    % calculate value to minimize
    F = sum( (gly-gly_sim).^2 )
end
    
end