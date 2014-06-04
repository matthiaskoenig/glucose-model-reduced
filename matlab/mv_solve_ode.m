function [t_data, f_data, c_data] = mv_solve_ode()
% mv_solve_ode : Solves the reduced model with ODE solver.
% See mv_solve_euler for the solution by Euler method and desripton
% of model, variables and parameters.
%
% @author: Matthias Koenig (matthias.koenig[AT]charite.de)
% @date:   2013-07-18

% time course information
global modus tc
% modus = '1meal'
modus = '3meals'
if (strcmp(modus, '1meal'))
    tc = timecourse_1meal();
elseif (strcmp(modus, '3meals'))
    tc = timecourse_3meals();
end

global f_liquid  
global f_solid  
f_liquid = 0.2;                  % [0,1] liquid fraction of volume (sinusoids & Space of Disse)
f_solid  = 0.7;                  % [0,1] solid fraction of volume (hepatocytes)

C_glyc = 500;                    % [mmol/L] (90 mg/ml) glycogen storage density 
C_glyc_solid = C_glyc/f_solid;   % [mmol/L] glycogen storage density in solid 

% initial conditions
c_init = [
                9                 % glc_ext   C1 [mmol/L]
                0.5*C_glyc_solid  % glyc      C2 [mmol/L]
                4                 % lac_ext   C3 [mmol/L] 
];

%% Integration
DT = 60*60*60;                  % [s] endtime simulation (DT)
[t_data, c_data] = ode15s(@mv_dxdt, [0 DT], c_init);
c_data = c_data';
f_data = zeros(size(c_data)); 

% Calculate the fluxes [mmol/l/s]
for k = 1:length(t_data)
    f_data(:,k) = mv_dxdt(t_data(k), c_data(:,k));
end

% Set the  concentrations from profile
if (strcmp(modus, '1meal') || strcmp(modus, '3meals'))
    c_data(1,:) = f_meal(t_data/3600, tc);    % [mM]
    c_data(3,:) = 3*ones(size(t_data));     % [mM]
end

%% Calculate additional data
V_sim = 1;      % [mL] simulation volume
c_data(4,:) = c_data(2,:)*f_solid;            % Glycogen concentration in simulation volume [0, 500]
c_data(5,:) = c_data(2,:)*f_solid*V_sim/1000; % [mmol] glycogen in simulation volume

%% plot the results
fig_mv(t_data, f_data, c_data, 'MV8 Matlab ODE');

return