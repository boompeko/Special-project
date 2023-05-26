
clear;
clc;



N = 10 ;%Number of rollers
Rr = 5 ;%Radius of the roller
R = 40 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E =2 ;% Eccentricity - offset from input shaft to a cycloidal disk

for i = 1 : 1 : 361

    t(i) = i / 180 * pi;
    


    x(i) = (R*cos(t(i)))-(Rr*cos(t(i)+atan(sin((1-N)*t(i))/((R/(E*N))-cos((1-N)*t(i))))))-(E*cos(N*t(i)));
    y(i) = (-R*sin(t(i)))+(Rr*sin(t(i)+atan(sin((1-N)*t(i))/((R/(E*N))-cos((1-N)*t(i))))))+(E*sin(N*t(i)));

end

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