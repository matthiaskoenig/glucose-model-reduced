function [t_data, f_data, c_data] = mv_solve_euler()
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

% endtime and timesteps (DT and DTC)
% time steps change the accuracy of the solver, but should not change
% the solutions as long as DCT is small compared to the kinetics
format compact; clear all
DT = 0.2*60*60;           % [s] endtime
DTC = 1;                  % [s] size of time steps

t_data = zeros(1, DT/DTC+1); 
for k = 2:length(t_data)
    t_data(1,k) = t_data(1,k-1) + DTC;
end

% Simulation Volume in [ml]
% V_sim = 1500; % [ml] total human liver volume
% V_sim = 0.9;  % [ml] total mouse liver volume
V_sim = 1;      % [mL] simulation volume

% Every simulation volume consists of different fractions [f_i]
global f_liquid  f_solid  
f_liquid = 0.2;                  % [0,1] liquid fraction of volume (sinusoids & Space of Disse)
f_solid  = 0.7;                  % [0,1] solid fraction of volume (hepatocytes)
f_rest   = 1-f_liquid-f_solid;   % [0,1] volume fraction not in solid and liquid 

% Every simulation volume has a maximum glycogen density (capacity)
% calculated from human model and in line with rat/mouse: 500mmmol/L -> 90mg/gLW -> 90mg/[ml_Vsim]
% The storage capacity is fully in the solid phase;
C_glyc = 500;                    % [mmol/L] (90 mg/ml) glycogen storage density 
C_glyc_solid = C_glyc/f_solid;   % [mmol/L] glycogen storage density in solid 

% storage of timepoints, concentrations and fluxes (F)
NC = 3;  % number of concentrations

% initial conditions of simulation volume
c_data = zeros(NC, length(t_data));
f_data = zeros(NC, length(t_data));
t_data(1,1) = 0;
c_data(:,1) = [
                9                  % [mmol/L]  glc_ext   C1
                0.5*C_glyc_solid   % [mmol/L]  glyc      C2  (half-filled glycogen stores solid phase)
                4                  % [mmol/L]  lac_ext   C3       
];

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

%% Calculate additional data
c_data(4,:) = c_data(2,:)*f_solid;            % Glycogen concentration in simulation volume [0, 500]
c_data(5,:) = c_data(2,:)*f_solid*V_sim/1000; % [mmol] glycogen in simulation volume

% plot the results
fig_mv(t_data, f_data, c_data, 'MV8 Matlab Euler');


return
