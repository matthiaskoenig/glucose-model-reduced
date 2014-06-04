% Test the fitting of 3D polynoms with test data

%% Generate the test data set and plot the resulting data
N = 70;
 x = linspace(0,12,N+1);
 y = linspace(0,15,2*N+1);
 z      = zeros(length(x), length(y));
 z_ns   = zeros(length(x), length(y));
 z_fit  = zeros(length(x), length(y));
 z_fit_ns  = zeros(length(x), length(y));
 
 mx = 3;
 my = 5;
 sigmax = 3;
 sigmay = 0.2;
 ns = 0.4;
 for i = 1:length(x)
     for k = 1:length(y)
        z(i,k) = exp(-(mx-x(i))^2/sigmax^2).*exp(-(mx-y(k))^2/sigmax^2);
        z_ns(i,k) = z(i,k) + ns * z(i,k) * (rand(1)-0.5);
     end
 end

 
 
%% Fit the test data set and plot the results
indepvar  = zeros(length(x)*length(y), 2);
depvar    = zeros(length(x)*length(y), 1);
depvar_ns = zeros(length(x)*length(y), 1);
counter = 1;
for i = 1:length(x)
    for k = 1:length(y)
        indepvar(counter,1) = x(i);
        indepvar(counter,2) = y(k);
        depvar(counter,1) = z(i,k);
        depvar_ns(counter,1) = z_ns(i,k);
        counter = counter + 1;
    end
end

%{
figure(1)
subplot(1,2,1)
plot(indepvar(:,1), depvar(:), 'o')
subplot(1,2,2)
plot(indepvar(:,2), depvar(:), 'o')

figure(2)
subplot(1,3,1)
plot(indepvar(:,1), 'o')
subplot(1,3,2)
plot(indepvar(:,2), 'o')
subplot(1,3,3)
plot(depvar(:), 'o')
%}

% Fit polynom and calculation of the actual values
 modelterms = 4;
 polymodel = polyfitn(indepvar,depvar,modelterms);
 polymodel.ModelTerms
 polymodel.Coefficients
 
 ypred = polyvaln(polymodel, indepvar);
 counter = 1
 for i = 1:length(x)
     for k = 1:length(y)
         z_fit(i,k) = ypred(counter);
         counter = counter + 1;
     end
 end
 
polymodel = polyfitn(indepvar,depvar_ns,modelterms);
 ypred_ns = polyvaln(polymodel, indepvar);
 counter = 1
 for i = 1:length(x)
     for k = 1:length(y)
         z_fit(i,k) = ypred(counter);
         z_fit_ns(i,k) = ypred_ns(counter);
         counter = counter + 1;
     end
 end
 
 %{
 % 2nd order
 c2 = [-0.0019 0.0063 -0.0535 8.2924e-04 -0.0792 0.7911];
 f2 = @(x,y) c2(1)*x^2 + c2(2)*x.*y + c2(3)*x + c2(4)* y^2 + c2(5)*y + c2(6)
 z_fit_f2  = zeros(length(x), length(y));
 for i = 1:length(x)
     for k = 1:length(y)
         z_fit_f2(i,k) = f2(x(i), y(k));
     end
 end
 z_fit = z_fit_f2;

 % 3rd order
 c3 = [0.0017    0.0004   -0.0363   -0.0002    0.0041    0.1263    0.0009   -0.0191    0.0402  0.4598];
 f3 = @(x,y) c3(1)*x^3 + c3(2)*x^2.*y + c3(3)*x^2 + c3(4)*x*y^2 + c3(5)*x*y + c3(6)*x ...
           + c3(7)*y^3 + c3(8)*y^2 + c3(9)*y + c3(10);
 z_fit_f3  = zeros(length(x), length(y));
 for i = 1:length(x)
     for k = 1:length(y)
         z_fit_f3(i,k) = f3(x(i), y(k));
     end
 end
 z_fit = z_fit_f3;
 
 %4rd order
%} 
c4 = [-0.000179592381288  -0.000383112240982   0.008928543667041  -0.000011722582207 ...
       0.007481166159607  -0.121596115247509  -0.000203421457280   0.004537379595064 ...
      -0.058254162191878   0.498464510771109  -0.000195079825952   0.008008634466379 ...
        -0.103182874252539   0.426333758596291  -0.173671700212682];
 f4 = @(x,y) c4(1)*x^4 + c4(2)*x^3*y + c4(3)*x^3 + c4(4)*x^2*y^2 + c4(5)*x^2*y + c4(6)*x^2 ...
           + c4(7)*x*y^3 + c4(8)*x*y^2 + c4(9)*x*y + c4(10)*x + c4(11)*y^4 ...
           + c4(12)*y^3 + c4(13)*y^2 + c4(14)*y + c4(15);
 z_fit_f4  = zeros(length(x), length(y));
 for i = 1:length(x)
     for k = 1:length(y)
         z_fit_f4(i,k) = f4(x(i), y(k));
     end
 end
 z_fit = z_fit_f4;
 
 
 %% Plot the test data and the resulting fit
 edge_alpha = 0.1;
 surf_axis  = [0 12 0 15 0 1]
 figure(1)
 figure('Name', '3D polynomal fit', 'Color', [1 1 1], 'OuterPosition', [0 0 800 1200])
 for k=1:6
    subplot(3,2,k)
    switch k
        case 1
            p1 = surf(x,y,z');
            title('data', 'FontWeight', 'bold')
        case 2
            p1 = surf(x,y,z_ns');
            title('data noise', 'FontWeight', 'bold')
        case 3
            p1 = surf(x,y,z_fit');
            title('data fit', 'FontWeight', 'bold')       
        case 4
            p1 = surf(x,y,z_fit_ns');
            title('data fit', 'FontWeight', 'bold')   
        case 5
            p1 = surf(x,y,z'-z_fit');
            title('data fit', 'FontWeight', 'bold')       
        case 6
            p1 = surf(x,y,z_ns'-z_fit_ns');
            title('data fit', 'FontWeight', 'bold')       
    end
    set(p1, 'EdgeAlpha', edge_alpha)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis(surf_axis)
    axis square
 end
 
