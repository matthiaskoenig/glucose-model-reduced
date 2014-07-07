function [t_data, f_data, c_data] = mv_solve_flow(modus_index)
% mv_solve_ode : Solves the reduced model with ODE solver.
% See mv_solve_euler for the solution by Euler method and desripton
% of model, variables and parameters.
%
% @author: Matthias Koenig (matthias.koenig[AT]charite.de)
% @date:   2014-07-07

VERSION = 14;
addpath('../glucose-profiles/')

global modus tc
modus_sel = {'stationary', '1meal', '3meals', 'sinus'};
if (nargin == 0)
    modus = modus_sel(1);
else
   modus = modus_sel(modus_index); 
end
if (strcmp(modus, '1meal'))
    tc = timecourse_1meal();
elseif (strcmp(modus, '3meals'))
    tc = timecourse_3meals();
elseif (strcmp(modus, 'sinus'))
    tc = timecourse_sinus();
end
disp(modus)

global f_liquid  f_solid
f_liquid = 0.2;                  % [0,1] liquid fraction of volume (sinusoids & Space of Disse)
f_solid  = 0.7;                  % [0,1] solid fraction of volume (hepatocytes)

% stationary initial conditions
C_glyc = 500;                    % [mmol/L] (90 mg/ml) glycogen storage density 

c_init = [
                3                   % glc_pp   [mmol/L]
                4                   % lac_pp   [mmol/L]
                3                   % glc_ext   [mmol/L]
                4                   % lac_ext   [mmol/L]
                0.5*C_glyc/f_solid  % glyc      [mmol/L]
];

% initial conditions for timecourse profiles
if (strcmp(modus, '1meal') || strcmp(modus, '3meals') || strcmp(modus, 'sinus'))
    [glc, lac] = f_timecourse(0/3600, tc);
    c_init(1,:) = glc;  % [mM] periportal
    c_init(2,:) = lac;  % [mM] periportal
    c_init(3,:) = glc;  % [mM] extern
    c_init(4,:) = lac;  % [mM] extern
end

%% Integration
DT = 60*3600;             % [s] endtime simulation (DT)
% DT = 10*3600;                  % [s] endtime simulation (DT)

                 
options = odeset('AbsTol', 1E-6, 'RelTol', 1E-6);
%[t_data, c_data] = ode15s(@mv_dxdt_flow, [0:60:DT], c_init, options);
[t_data, c_data] = ode15s(@mv_dxdt_flow, [0, DT], c_init, options);

c_data = c_data';
f_data = zeros(size(c_data)); 

% Set the concentrations from profile
if (strcmp(modus, '1meal') || strcmp(modus, '3meals') || strcmp(modus, 'sinus'))
    [glc, lac] = f_timecourse(t_data/3600, tc);
    c_data(1,:) = glc;    % [mM]
    c_data(2,:) = lac;    % [mM]
end

% Calculate the fluxes [mmol/l/s]
for k = 1:length(t_data)
    f_data(:,k) = mv_dxdt_flow(t_data(k), c_data(:,k));
end
% glycogen per volume
c_data(6,:) = c_data(5,:)*f_solid;  % [mM] glycogen concentration in reference volume [0, 500]

if (1)

%% store the results
headers = {'time', 'glc_pp', 'lac_pp', 'glc', 'lac', 'gly', 'gly_vol', ...
                 'f_glc_pp', 'f_lac_pp', 'f_glc', 'f_lac', 'f_gly'} 
res =[t_data, c_data', f_data'];
fname = strcat('MV', int2str(VERSION), '_Matlab_ODE-', modus, '.csv')
csvwrite_with_headers(fname{1}, res, headers);

%% plot the results
fig_name = strcat('MV', int2str(VERSION), '_Matlab_ODE-', modus)
fig_integration_flow(t_data, f_data, c_data, fig_name{1});

end
