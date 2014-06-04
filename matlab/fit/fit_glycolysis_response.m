% Fit the glycolysis response curves with polynoms to reuse in other models

%% Load the data for analysis (dependendy glucose and glycogen)
out_folder = '/home/mkoenig/Desktop/kinetic_model_results/simulations/110314_data/glucose_glycogen_dependency/T2DM'
conditions = {'normal', 'diabetes', 'insulin_restored', 'glucagon_restored'}

condition = 'normal' 
in_file = strcat(out_folder,  '/glucose_glycogen_dependency_', condition);
load(in_file);
    
glc_ext  = mpars.glc_ext;
glycogen = mpars.glycogen;
t_span   = mpars.t_span;

% indeces of boundaries
    t_eval = 100;           % time point of evaluation
    glc_min = glc_ext(1); glc_max = glc_ext(end);
    glycogen_min = glycogen(1); glycogen_max = glycogen(end);
    
    tmp = find(glc_ext>=glc_min); glc_min_ind = tmp(1);
    tmp = find(glc_ext<=glc_max); glc_max_ind = tmp(end);
    tmp = find(glycogen>=glycogen_min); glycogen_min_ind = tmp(1);
    tmp = find(glycogen<=glycogen_max); glycogen_max_ind = tmp(end);
    tmp = find(glycogen>=250); glycogen_250_ind = tmp(1);
    tmp = find(t_span>=t_eval); t_eval_ind = tmp(1);

    x = glc_ext(glc_min_ind:glc_max_ind);
    y = glycogen(glycogen_min_ind:glycogen_max_ind);
    
    % panel 1 - HGP and HGU
    z1 = v_full(glc_min_ind:glc_max_ind, glycogen_min_ind:glycogen_max_ind, t_eval_ind, 1);

    % panel 2 - glycolysis and gluconeogenesis
    z2 = v_full(glc_min_ind:glc_max_ind, glycogen_min_ind:glycogen_max_ind,  t_eval_ind, 4);
    
    % panel 3 - glycogenolyis and glycogen synthesis
    z3 = -v_full(glc_min_ind:glc_max_ind, glycogen_min_ind:glycogen_max_ind, t_eval_ind, 5);
    

 
%% Fit the test data sets and plot the results
fit_z1 = zeros(length(x), length(y));
fit_z2 = zeros(length(x), length(y));
fit_z3 = zeros(length(x), length(y));
indepvar  = zeros(length(x)*length(y), 2);
depvar_z1 = zeros(length(x)*length(y), 1);
depvar_z2 = zeros(length(x)*length(y), 1);
depvar_z3 = zeros(length(x)*length(y), 1);
counter = 1;
for i = 1:length(x)
    for k = 1:length(y)
        indepvar(counter,1) = x(i);
        indepvar(counter,2) = y(k);
        depvar_z1(counter,1) = z1(i,k);
        depvar_z2(counter,1) = z2(i,k);
        depvar_z3(counter,1) = z3(i,k);
        counter = counter + 1;
    end
end

% Fit polynom and calculation of the actual values
 modelterms = 3;
 polymodel_z1 = polyfitn(indepvar,depvar_z1,modelterms)
 polymodel_z1.Coefficients
 polymodel_z2 = polyfitn(indepvar,depvar_z2,modelterms)
 polymodel_z2.Coefficients
 polymodel_z3 = polyfitn(indepvar,depvar_z3,modelterms)
 polymodel_z3.Coefficients
 polymodel_z3.ModelTerms
 ypred_z1 = polyvaln(polymodel_z1, indepvar);
 ypred_z2 = polyvaln(polymodel_z2, indepvar);
 ypred_z3 = polyvaln(polymodel_z3, indepvar);
 counter = 1
 
 % Compare with calculated fit curves
 fs = fit_curves_glucose_metabolism();
 
 for i = 1:length(x)
     for k = 1:length(y)
         fit_z1(i,k) = ypred_z1(counter);
         fit_z1(i,k) = fs.z1(x(i), y(k));
         fit_z2(i,k) = ypred_z2(counter);
         fit_z2(i,k) = fs.z2(x(i), y(k));
         fit_z3(i,k) = ypred_z3(counter);
         fit_z3(i,k) = fs.z3(x(i), y(k));
         counter = counter + 1;
     end
 end

 
 %% Plot the test data and the resulting fit
 edge_alpha = 0.1;
 surf_axis  = [0 12 0 15 0 1]
 figure(1)
 figure('Name', '3D polynomal fit', 'Color', [1 1 1], 'OuterPosition', [0 0 800 1600])
 for k=1:12
    subplot(3,4,k)
    switch k
        case 1
            p1 = surf(x,y,z1');
            title('z1', 'FontWeight', 'bold')
        case 2
            p1 = surf(x,y,z2');
            title('z2', 'FontWeight', 'bold')
        case 3
            p1 = surf(x,y,z3');
            title('z3', 'FontWeight', 'bold')       
        case 4
            p1 = surf(x,y,z2' + z3');
            title('z2 + z3', 'FontWeight', 'bold')       
        
            
        case 5
            p1 = surf(x,y,fit_z1');
            title('fit z1', 'FontWeight', 'bold')
        case 6
            p1 = surf(x,y,fit_z2');
            title('fit z2', 'FontWeight', 'bold')
        case 7
            p1 = surf(x,y,fit_z3');
            title('fit z3', 'FontWeight', 'bold')       
        case 8
            p1 = surf(x,y,fit_z2' + fit_z3');
            title('fit z3', 'FontWeight', 'bold')       
            
        case 9
            p1 = surf(x,y,z1'-fit_z1');
            title('diff', 'FontWeight', 'bold')
        case 10
            p1 = surf(x,y,z2'-fit_z2');
            title('diff', 'FontWeight', 'bold')
        case 11
            p1 = surf(x,y,z3'-fit_z3');
            title('diff', 'FontWeight', 'bold')       
                 
    end
    set(p1, 'EdgeAlpha', edge_alpha)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    %axis(surf_axis)
    axis square
 end
 
