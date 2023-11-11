
clear;
clc;



N = 24 ;%Number of rollers
Rr = 9 ;%Radius of the roller
R = 150 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E =3 ;% Eccentricity - offset from input shaft to a cycloidal disk

for i = 1 : 1 : 361

    t(i) = i / 180 * pi;
    


    x1(i) = (R*cos(t(i)))-(Rr*cos(t(i)+atan(sin((1-N)*t(i))/((R/(E*N))-cos((1-N)*t(i))))))-(E*cos(N*t(i)));
    y1(i) = -(R*sin(t(i)))+(Rr*sin(t(i)+atan(sin((1-N)*t(i))/((R/(E*N))-cos((1-N)*t(i))))))+(E*sin(N*t(i)));
    
    x2(i) = (R*cos(t(i)))+Rr*cos(t(i)+atan(sin((1+N)*t(i))/((R/(E*N))-cos((1+N)*t(i)))))+(E*cos(N*t(i)));
    y2(i) = (R*sin(t(i)))+Rr*sin(t(i)+atan(sin((1+N)*t(i))/((R/(E*N))-cos((1+N)*t(i)))))-(E*sin(N*t(i)));

end

figure;
plot(x1,y1,'LineWidth',2);
xlabel('X','fontname','Times New Roman','fontsize',20');
ylabel('Y','fontname','Times New Roman','fontsize',20');
axis equal;
xlim([-(1.2*R+10),(1.2*R+10)]);
ylim([-(1.2*R+10),(1.2*R+10)]);
xticks(-(1.2*R+10):10:(1.2*R+10));
yticks(-(1.2*R+10):10:(1.2*R+10));


box on;
grid on;

figure;
hold on
plot(x2,y2,'LineWidth',2);
xlabel('X','fontname','Times New Roman','fontsize',20');
ylabel('Y','fontname','Times New Roman','fontsize',20');
axis equal;
xlim([-(1.2*R+10),(1.2*R+10)]);
ylim([-(1.2*R+10),(1.2*R+10)]);
xticks(-(1.2*R+10):10:(1.2*R+10));
yticks(-(1.2*R+10):10:(1.2*R+10));


hold off
box on;
grid on;


% fid = fopen('Cycloidal_Drive.scr','w');
% fprintf(fid,'spline ');
% fprintf(fid,'%f,%f\n', [x; y]);
% fclose(fid);