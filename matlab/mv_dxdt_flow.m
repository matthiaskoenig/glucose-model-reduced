function [dxdt] = mv_dxdt_flow(t, x)
%% mv_dxdt_flow : ODE rate laws for model
%  Definition of the rate laws and calculation of the fluxes for the given
%  parameters and concentrations. This function is called by the ODE
%  integration routine.
%  
%   dxdt                vector of fluxes (changes in ODE)
%   c                   vector of concentrations

% Units
%   time: seconds     [sec] ,  
%   concentration     [mM] = [mmol/l] 
%   volume:           [litre]
%   dxdt:             [mmol/s/litre]
%
%   author: Matthias Koenig 
%           Charite Berlin
%           Computational Systems Biochemistry Berlin
%           matthias.koenig@charite.de
%   date:   2014-07-07


% [A] hard coded time profiles in the periportal compartment, i.e.
% set the concentrations from profile
global modus tc
if (strcmp(modus, '1meal') || strcmp(modus, '3meals') || strcmp(modus, 'sinus'))
    [glc, lac] = f_timecourse(t/3600, tc);
    x(1) = glc;    % [mM] glucose periportal (V_pp)
    x(2) = lac;    % [mM] lactate periportal (V_pp)
end

% [B] Calculate changes
dxdt = zeros(size(x));

% [B2] metabolic changes
global f_liquid f_solid

% Changes are relative to the respective volumes
% hgu, gs, gly are [mmol/s/L] simulation Volume
% the changes in the liquid and solid volumes are calculated with the 
% respective fractions

% glc_ext (x3), glyc (x5), lac_ext (x4)
f = fit_kinetics_flow();
hgu = f.hgu(x(3), x(5), x(4));
gly = f.gly(x(3), x(5), x(4));
gs = f.gs(x(3), x(5), x(4));
dxdt(3) =  -1 * hgu/f_liquid; % [mmol/s/L]  glucose V_liquid 
dxdt(4) =  +2 * gly/f_liquid; % [mmol/s/L]  lactate V_liquid
dxdt(5) =  +1 * gs/f_solid;   % [mmol/s/L]  glycogen V_solid                

% [B1] flow changes (glucose and lactate)
L_sin = 550E-6;           % [m] sinusoid length
v_flow = 180E-6;          % [m/s] sinusoidal flow
t_transit = L_sin/v_flow; % [s] transit time  
k_flow = 1/t_transit;     % [1/s] change rate due to flow (alpha)
% k_flow = 1/13.577511022376964;

dxdt(3) = dxdt(3) + k_flow*(x(1) - x(3));  % [mmol/s/L] glucose V_liquid  
dxdt(4) = dxdt(4) + k_flow*(x(2) - x(4));  % [mmol/s/L] lactate V_liquid

% periportal concentrations assumed constant (given by profile)
dxdt(1) = 0;       % [mmol/s/L]  glucose V_periportal 
dxdt(2) = 0;       % [mmol/s/L]  lactate V_periportal 


