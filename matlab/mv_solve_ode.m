function [t_data, f_data, c_data] = mv_solve_ode(modus_index)
% mv_solve_ode : Solves the reduced model with ODE solver.
% See mv_solve_euler for the solution by Euler method and desripton
% of model, variables and parameters.
%
% @author: Matthias Koenig (matthias.koenig[AT]charite.de)
% @date:   2014-06-04

VERSION = 12;

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

global f_liquid  f_solid V_sim
f_liquid = 0.25;                  % [0,1] liquid fraction of volume (sinusoids & Space of Disse)
f_solid  = 0.75;                  % [0,1] solid fraction of volume (hepatocytes)
V_sim = 0.05;                     % [L] simulation volume (reference volume of model V_ref=1.0L

% stationary initial conditions
C_glyc = 500;                    % [mmol/L] (90 mg/ml) glycogen storage density 
c_init = [
                9                   % glc_ext   C1 [mmol/L]
                0.5*C_glyc/f_solid  % glyc      C2 [mmol/L]
                4                   % lac_ext   C3 [mmol/L] 
];

% initial conditions for timecourse profiles
if (strcmp(modus, '1meal') || strcmp(modus, '3meals') || strcmp(modus, 'sinus'))
    [glc, lac] = f_timecourse(0/3600, tc);
    c_init(1,:) = glc;  % [mM]
    c_init(3,:) = lac;  % [mM]
end

%% Integration
%DT = 60*3600;                  % [s] endtime simulation (DT)
DT = 10;                  % [s] endtime simulation (DT)
                 
options = odeset('AbsTol', 1E-6, 'RelTol', 1E-6);
[t_data, c_data] = ode15s(@mv_dxdt, [0:10:DT], c_init, options);
% tic
% [t,c] = ode15s(@mv_dxdt, [0, 2000], c_init, odeset('RelTol', 1e-9, 'AbsTol', 1e-9));
% toc
c_data = c_data';
f_data = zeros(size(c_data)); 

% Set the concentrations from profile
if (strcmp(modus, '1meal') || strcmp(modus, '3meals') || strcmp(modus, 'sinus'))
    [glc, lac] = f_timecourse(t_data/3600, tc);
    c_data(1,:) = glc;    % [mM]
    c_data(3,:) = lac;    % [mM]
end

% Calculate the fluxes [mmol/l/s]
for k = 1:length(t_data)
    f_data(:,k) = mv_dxdt(t_data(k), c_data(:,k));
end
% glycogen per volume
c_data(4,:) = c_data(2,:)*f_solid;  % [mM] glycogen concentration in reference volume [0, 500]

%% store the results
headers = {'time', 'glc', 'gly', 'lac', 'gly_vol', 'f_glc', 'f_gly', 'f_lac'} 
res =[t_data, c_data', f_data'];
fname = strcat('../results/','MV', int2str(VERSION), '_Matlab_ODE-', modus, '.csv')
csvwrite_with_headers(fname{1}, res, headers);

%% plot the results
fig_name = strcat('MV', int2str(VERSION), '_Matlab_ODE-', modus)
fig_integration(t_data, f_data, c_data, fig_name{1});

end