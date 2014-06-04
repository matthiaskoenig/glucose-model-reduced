% Plot the fit kinetics of the reduced glucose model;
%
% Matthias Koenig (matthias.koenig[AT]charite.de)
% 2013-06-04


close all
clear all
fprintf('Reduced Glucose Model Response\n')
Vf = 1;  % Volume factor

f = fit_kinetics(Vf);
x1 = 0:0.5:12;   % glc_ext
x2 = 0:10:500;   % glycogen
x3 = 0:0.1:4;    % lac_ext

% Precalculate all the values for the plot
hgu = zeros(length(x1), length(x2), length(x3));
gly = zeros(length(x1), length(x2), length(x3));
gs = zeros(length(x1), length(x2), length(x3));
for k=1:length(x1)
    for p=1:length(x2)
        for q=1:length(x3)
            hgu(k,p,q) = f.hgu(x1(k),x2(p),x3(q));
            gly(k,p,q) = f.gly(x1(k),x2(p),x3(q));
            gs(k,p,q) = f.gs(x1(k),x2(p),x3(q));
        end
    end
end

% Generate the different subplots
figure = figure('Name', 'Polynomal Fits', 'Color',[1 1 1]', ...
    'Position', [100 100 1000 800]);
for k=1:9
   subplot(3,3,k)
   switch k
       case 1
           data = hgu(:,:,3);
           tit = 'HGU';
       case 2
           data = gly(:,:,3);
           tit = 'GLY';
       case 3
           data = gs(:,:,3);
           tit = 'GS';
       case 4
           data = hgu(:,:,2);
           tit = 'HGU';
       case 5
           data = gly(:,:,2);
           tit = 'GLY';
       case 6
           data = gs(:,:,2);
           tit = 'GS';
       case 7
           data = hgu(:,:,1);
           tit = 'HGU';
       case 8
           data = gly(:,:,1);
           tit = 'GLY';
       case 9
           data = gs(:,:,1);
           tit = 'GS';
   end
   data = squeeze(data);
   p = surf(data)
   p = pcolor(x1, x2, data');
   set(p, 'EdgeAlpha', 0.1);
   title(tit, 'FontWeight', 'bold')
   xlabel('glc_{ext}', 'FontWeight', 'bold')
   ylabel('glyc', 'FontWeight', 'bold')
   %colorbar
   axis square
end
