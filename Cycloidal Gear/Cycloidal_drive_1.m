
clear;
clc;



N = 13 ;%Number of rollers
Rr = 5 ;%Radius of the roller
R = 40 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E =2 ;% Eccentricity - offset from input shaft to a cycloidal disk



for i = 1 : 1 : 361

    t(i) = i / 180 * pi;

    phi(i) = atan(sin((1-N)*t(i))/((R/(E*N))-cos((1-N)*t(i))));
    
    Cf(i, :, :) = [R-Rr*cos(phi(i));Rr*sin(phi(i));1];
    
    M1(i, :, :) = [cos(-N*t(i)), -sin(-N*t(i)), -E*cos(-N*t(i));
                   sin(-N*t(i)), cos(-N*t(i)), -E*sin(-N*t(i));
                   0, 0, 1];
    M2(i, :, :) = [cos(1-N*phi(i)), sin(1-N*phi(i)), 0;
                   -sin(1-N*phi(i)), cos(1-N*phi(i)), 0;
                   0, 0, 1];
    
    C(i, :, :) = M1(i, :, :).*M2(i, :, :).*Cf(i, :, :);

    
    

end

x = C(1,:);
y = C(2,:);

figure;
plot(x,y,'LineWidth',2);
xlabel('X','fontname','Times New Roman','fontsize',20');
ylabel('Y','fontname','Times New Roman','fontsize',20');
axis equal;
xlim([-50,50]);
ylim([-50,50]);
xticks(-50:10:50);
yticks(-50:10:50);


box on;
grid on;


% fid = fopen('Cycloidal_Drive.scr','w');
% fprintf(fid,'spline ');
% fprintf(fid,'%f,%f\n', [x; y]);
% fclose(fid);