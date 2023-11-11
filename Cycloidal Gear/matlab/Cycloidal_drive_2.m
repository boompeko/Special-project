
clear;
clc;



N = 13 ;%Number of rollers
Rr = 5 ;%Radius of the roller
R = 40 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E =2 ;% Eccentricity - offset from input shaft to a cycloidal disk



for i = 1 : 1 : 361

    t(i) = i / 180 * pi;

    phi(i) = atan(sin((1+N)*t(i))/((R/(E*N))-cos((1+N)*t(i))));
    
    Cf = [R+Rr*cos(phi(i));Rr*sin(phi(i));1];
    
    [Rc(i)] = Rrcal_hypocycloid(N,R,E,Rr,t(i));

  
    Ccf = [R+(Rr+Rc(i))*cos(phi(i));(Rr+Rc(i))*sin(phi(i));1];

    M1 = [cos(N*t(i)), sin(N*t(i)), 0;
          -sin(N*t(i)), cos(N*t(i)), 0;
          0, 0, 1];
    M2 = [cos((1+N)*t(i)), -sin((1+N)*t(i)), E;
          sin((1+N)*t(i)), cos((1+N)*t(i)), 0;
          0, 0, 1];
    
    C = M1*M2*Cf;
    Cc = M2*Ccf;
    
    x(i) = C(1,:);
    y(i) = C(2,:);

    xc(i) = Cc(1,:);
    yc(i) = Cc(2,:);
    

end


a = figure('Visible', 'on');
plot(x,y,'LineWidth',2);
xlabel('X','fontname','Times New Roman','fontsize',20');
ylabel('Y','fontname','Times New Roman','fontsize',20');
title('擺線輪輪廓');
axis equal;
xlim([-(R+2*Rr),(R+2*Rr)]);
ylim([-(R+2*Rr),(R+2*Rr)]);
xticks(-2*R:10:2*R);
yticks(-2*R:10:2*R);

box on;
grid on;
axis square;
