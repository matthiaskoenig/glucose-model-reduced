% Calculates the changes in glucose, lactate and glycogen for a
% reference volume of 1L liver; 

% All fluxes are at this point in [mmol/s] per 1L reference volume, 
% i.e. [mmol/L/s]

function [f] = fit_kinetics()
    global f_solid  
    
    f.hgu = @hgu;
    f.gly = @gly;
    f.gs  = @gs;
    
    % lac_ext becomes limiting for HGP, GS and GNG at low concentrations
    k_lac = 0.05;   % [mM]
    k_glc = 0.05;   % [mM]
    
    % x (glc_ext), y (glycogen), z (lac_ext)
    % HGU (>0) and HGP (<0)
    function [res] = hgu(x, y, z)   % [mmole/s/L]
        % The glycogen storage capacity depends on Vf and the glycogen
        % variable has to be transformed accordingly -> [0,500]
        y = y*f_solid; 
        C = [   0.002037960420379  -0.000367490977632  -0.069301032419012  ... 
               -0.000002823120484   0.011282864074433   3.740276159346358  ...
               -0.000000181515364   0.000157328485157   -0.100050917438436 ...
               -18.414978834613287];
           
        % kinetics depends mainly on glc_ext and glycogen_ext
        res = C(1)*x^3 + C(2)*x^2*y + C(3)*x^2 + C(4)*x*y^2 + C(5)*x*y ...
            +C(6)*x + C(7)*y^3 + C(8)*y^2 + C(9)*y + C(10);
        % HGP limited by lactate
        if (res < 0)
            res = res * z/(z + k_lac); 
        end
        % HGU limited by glc
        if (res > 0)
            res = res * x/(x + k_glc); 
        end
    end

    %GNG and GLY
    function [res] = gly(x, y, z)  % [mmole/s/L]
        % The glycogen storage capacity depends on Vf and the glycogen
        % variable has to be transformed accordingly;
        y = y*f_solid; 
        C = [  -0.013260401508103  -0.000078240970095   0.478235644004833  ...
                0.000002861605817   0.000932752106971  -2.492569641130055  ...
                0.000000166945924  -0.000125285017396   0.015354944655784  ...
               -4.975026288067225];
        res = C(1)*x^3 + C(2)*x^2*y + C(3)*x^2 + C(4)*x*y^2 + C(5)*x*y ...
            +C(6)*x + C(7)*y^3 + C(8)*y^2 + C(9)*y + C(10);
        
        % GNG limited by lactate
        if (res < 0)
            res = res * z/(z + k_lac); 
        end
        % GLY limited by glc
        if (res > 0)
            res = res * x/(x + k_glc); 
        end
    end

    %GS
    function [res] = gs(x, y, z)
        res = hgu(x,y,z)-gly(x, y, z);
    end
end