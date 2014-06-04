function [t_data, f_data, c_data] = mv_solve_euler(modus_index)
%  mv_solve_euler : Solves the reduced model in Matlab with Euler method.
%
%  Model v8 : Reduced Model of Hepatic Glucose Metabolism 
%  @author: Matthias Koenig (matthias.koenig[AT]charite.de)
%  @date:   2013-07-18
%
%     Koenig M., Bulik S. and Holzhuetter HG.
%     Quantifying the Contribution of the Liver to the Homeostasis of Plasma 
%     Glucose: A Detailed Kinetic Model of Hepatic Glucose Metabolism 
%     PLoS Comput Biol. 2012 Jun;8(6):e1002577. Epub 2012 Jun 21.
%
% Model reduction based on Quasi-Steady-State approximation of full kinetic
% model with 
%     c1 -> glc_ext (blood glucose) [mmol/L]
%     c2 -> glyc    (hepatic glycogen) [mmol/L]
%     c3 -> lac_ext (blood lactate) [mmol/L]
%
%     v1 -> HGU : HGU > 0 (hepatic glucose utilization)
%                 HGU < 0 (hepatic glucose prodcution)
%     v2 -> GLY : GLY > 0 (glycogenolysis)
%                 GLY < 0 (gluconeogenesis)
%     v1 -> GS  : GS  > 0 (glycogen synthesis)
%                 GS  < 0 (glycogenolysis)
% 
% Glycogen values are always between [0, 500] mM in full model
% corresponding to [0, 500]*VL/VSIM in simulation volume V_sim.

% -------------------------------------------------------------------------
% Model scenarios for simulation
% -------------------------------------------------------------------------
% 1. Switch between hepatic glucose production/hepatic glucose utilization
%       and glycogen utilization/glycogen storage via given glucose
%       profile.
%   Phase 1 [t0, t1) -> glc_ext = 3.5mM;  lac_ext = 3mM; glycogen(t0) = 1/2max
%   Phase 2 [t1, t2) -> glc_ext = 8mM;
%   Phase 3 [t2, t3) -> glc_ext = 3.5mM;
% expected result: in Phase 1 with low glucose, glucose is produced by
%                   the liver via gluconeogenesis and from glycogen stores.
%                   Glycogen is degraded in this phase.
%           
% 2. Complete filling of glycogen stores 
%    set in [t0, t1) -> glc_ext = 9mM; lac_ext = 3mM; glycogen(t0) = 0  
%
% 3. Complete emptying of glycogen stores 
%    set in [t0, t1) -> glc_ext = 3.5mM; lac_ext = 3mM; glycogen(t0) = min  
%
% 4. Complete filling of glycogen stores 
%    set in [t0, t1) -> glc_ext = 9mM; lac_ext = 3mM; glycogen(t0) = max
% 
% 5. Dynamic daily glucose profiles (1 meal & 3 meal)
% results: analog to the hard switch in Simulation 1, but more dynamic
%           behavior
% -------------------------------------------------------------------------

disp('EULER NOT UPDATED - USE THE ODE IMPLENTATION')
return

% endtime and timesteps (DT and DTC)
% time steps change the accuracy of the solver, but should not change
% the solutions as long as DCT is small compared to the kinetics
format compact; clear all
DT = 2*3600;             % [s] endtime
DTC = 1;                  % [s] size of time steps

t_data = zeros(1, DT/DTC+1); 
for k = 2:length(t_data)
    t_data(1,k) = t_data(1,k-1) + DTC;
end

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

C_glyc = 500;                    % [mmol/L] (90 mg/ml) glycogen storage density 
C_glyc_solid = C_glyc/f_solid;   % [mmol/L] glycogen storage density in solid 

% stationary initial conditions
c_init = [
                9                 % glc_ext   C1 [mmol/L]
                0.5*C_glyc_solid  % glyc      C2 [mmol/L]
                4                 % lac_ext   C3 [mmol/L] 
];
% Set initial concentrations from profile
if (strcmp(modus, '1meal') || strcmp(modus, '3meals') || strcmp(modus, 'sinus'))
    [glc, lac] = f_timecourse(0/3600, tc);
    c_init(1,:) = glc;  % [mM]
    c_init(3,:) = lac;  % [mM]
end


% initial conditions of simulation volume
NC = 3;  % number of concentrations
c_data = zeros(NC, length(t_data));
f_data = zeros(NC, length(t_data));
t_data(1,1) = 0;
c_data(:,1) = c_init;

%% Integration Euler
display('-----------------')
display('Time Matlab Euler')

tic
for k = 2:length(t_data)
    % get fluxes and concentrations for next time point
    f_data(:,k) = mv_dxdt(t_data(k-1), c_data(:, k-1));   % [mmol/s/l]

    % glucose in fluid
    c_data(1,k) = c_data(1, k-1) + DTC *f_data(1,k);       % [mmol/l]  
    % glycogen in solid
    c_data(2,k) = c_data(2, k-1) + DTC *f_data(2,k);       % [mmol/l]
    % lactate in fluid
    c_data(3,k) = c_data(3, k-1) + DTC *f_data(3,k);       % [mmol/l]
end
toc

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

% plot the results
fig_name = strcat('MV9_Matlab_Euler-', modus)
fig_integration(t_data, f_data, c_data, fig_name{1});

return
