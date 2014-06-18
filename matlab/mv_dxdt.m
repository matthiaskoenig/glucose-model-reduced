function [dxdt] = mv_dxdt(t, x)
%% mv_dxdt : ODE rate laws for model
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
%   date:   2014-06-04

global f_liquid f_solid V_sim
f = fit_kinetics();

% hard coded time profiles
global modus tc
% Set the concentrations from profile
if (strcmp(modus, '1meal') || strcmp(modus, '3meals') || strcmp(modus, 'sinus'))
    [glc, lac] = f_timecourse(t/3600, tc);
    x(1) = glc;    % [mM]
    x(3) = lac;    % [mM]
end

% glc_ext (x1), glyc (x2), lac_ext (x3)
hgu = f.hgu(x(1), x(2), x(3));
gly = f.gly(x(1), x(2), x(3));
gs = f.gs(x(1), x(2), x(3));

% Changes are relative to the respective volumes
% hgu, gs, gly are [mmol/s/L] simulation Volume
% the changes in the liquid and solid volumes are calculated with the 
% respective fractions
V_ref = 1.0;                % [l] reference volume 
dxdt = zeros(size(x));
dxdt(1) = -1 * hgu/f_liquid * V_sim/V_ref;   % [mmol/s/L] 
dxdt(2) = +1 * gs/f_solid   * V_sim/V_ref;   % [mmol/s/L]                  
dxdt(3) = +2 * gly/f_liquid * V_sim/V_ref;   % [mmol/s/L]


% test with constant glucose and lactate concentrations (liquid)
%dxdt(1) = 0;                
%dxdt(3) = 0;


