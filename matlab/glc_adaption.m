glc = linspace(4,8,40)
figure()
plot(glc, glc, '-o')
close all
figure()

glc_new =  glc + 2.5*(glc - 6.2)/6.2
glc_new2 =  glc + 5.0*(glc - 6.2)/6.2
plot(glc, glc, 'k-'), hold on
plot(glc, glc_new, 'r-o'), hold on
plot(glc, glc_new2, 'b-o'), hold on
%plot(glc, glc.^2 - 5.5*glc+5.5, 'b-o'), hold on
grid on
