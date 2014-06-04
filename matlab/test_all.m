% Call various scripts for testing.
%
% Matthias Koenig
% 2014-06-04

% [1] fit kinetics
fig_fit_kinetics

% [2] timecourse information
fig_timecourses

% [4] timecourse integration
mv_solve_ode(1) % stationary
mv_solve_ode(2) % 1 meal
mv_solve_ode(3) % 3 meals
mv_solve_ode(4) % sinus