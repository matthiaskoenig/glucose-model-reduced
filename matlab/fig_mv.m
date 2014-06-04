function [y] = fig_mv7(t, f, c, ftext)
% plot details about the solution
figure('Name', ftext, 'Color', [1 1 1], 'OuterPosition', [0 0 1200 800])
t_unit = 'min'
switch t_unit
    case 's'
        t = t;
    case 'min'
        t = t/60;
    case 'h'
        t = t/3600;
end

Nc = size(f,1);

% plot the concentrations
for k=1:Nc
    subplot(2,Nc,k)
    data = c(k,:);    
    plot(t, data); hold on;
    switch k
        case 1
            title('Glucose extern', 'FontWeight', 'bold')
        case 2
            title('Glycogen', 'FontWeight', 'bold')
        case 3
            title('Lactate extern', 'FontWeight', 'bold')
    end
    axis square
    xlabel(strcat('time [', t_unit, ']'))
    ylabel('C [mM]')
end
% plot the additional concentrations in the plot
subplot(2,3,2)
plot(t, c(4,:), 'k-'), hold on; 
plot(t, c(5,:), 'r-'), hold on;


% plot the fdata
for k=1:Nc
    subplot(2,Nc,Nc+k)
    data = f(k,:);    
    plot(t, data);
    switch k
        case 1
            title('fdata glc_{ext}', 'FontWeight', 'bold')
        case 2
            title('fdata glyc', 'FontWeight', 'bold')
        case 3
            title('fdata lac_{ext}', 'FontWeight', 'bold')
    end
    axis square
    xlabel(strcat('time [', t_unit, ']'))
    ylabel('v [mM/s]')
end

