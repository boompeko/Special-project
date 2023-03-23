
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              parameter input                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rb = 40; 
h = 20; 
beta = 100 / 180 * pi;
beta2 = 110 / 180 * pi;
theta_D = (190 / 180) * pi;
rf = 10;
offset = 12;


for i = 1 : 1 : 3600
    theta(i)=i*2*pi/3600;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Rise¡G0~120                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1000
    S(i)=h*theta(i)/beta-h/(2*pi)*sin(2*pi*theta(i)/beta);
    L(i)=((rf+rb)^2-offset^2)^0.5+S(i);
    q(i)=h/beta-h/beta*cos(2*pi*theta(i)/beta);
    a(i)=2*pi*h/(beta^2)*sin(2*pi*theta(i)/beta);
    phi(i) = atan((q(i)-offset)/L(i));
    X(i)=L(i)*cos(theta(i))-offset*sin(theta(i))-rf*cos(theta(i)-phi(i));
    Y(i)=L(i)*sin(theta(i))+offset*cos(theta(i))-rf*sin(theta(i)-phi(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Fall¡G190~310                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1901:3000
    S(i)=h-(h*(theta(i)-theta_D)/(beta2)-h/(2*pi)*sin(2*pi*(theta(i)-theta_D)/(beta2)));
    L(i)=((rf+rb)^2-offset^2)^0.5+S(i);
    q(i)=-((h/(beta2)-h/(beta2)*cos(2*pi*(theta(i)-theta_D)/(beta2))));
    a(i)=-2*pi*h/((beta2)^2)*sin(2*pi*(theta(i)-theta_D)/(beta2));
    phi(i) = atan((q(i)-offset)/L(i));
    X(i)=L(i)*cos(theta(i))-offset*sin(theta(i))-rf*cos(theta(i)-phi(i));
    Y(i)=L(i)*sin(theta(i))+offset*cos(theta(i))-rf*sin(theta(i)-phi(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           High Dwell¡G120~190                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1001:1900
    S(i)=h;
    L(i)=((rf+rb)^2-offset^2)^0.5+S(i);
    q(i)=0;
    a(i)=0;
    phi(i) = atan((q(i)-offset)/L(i));
    X(i)=L(i)*cos(theta(i))-offset*sin(theta(i))-rf*cos(theta(i)-phi(i));
    Y(i)=L(i)*sin(theta(i))+offset*cos(theta(i))-rf*sin(theta(i)-phi(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Low Dwell¡G 300~360                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=3001:3600
    S(i)=0;
    L(i)=((rf+rb)^2-offset^2)^0.5+S(i);
    q(i)=0;
    a(i)=0;
    phi(i) = atan((q(i)-offset)/L(i));
    X(i)=L(i)*cos(theta(i))-offset*sin(theta(i))-rf*cos(theta(i)-phi(i));
    Y(i)=L(i)*sin(theta(i))+offset*cos(theta(i))-rf*sin(theta(i)-phi(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Generate contour figure                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
axis equal;
plot(theta/(2*pi)*360,S,'LineWidth',2);
xlabel('\theta  (deg)','fontname','Times New Roman','fontsize',20');
ylabel('S  (\theta)(mm)','fontname','Times New Roman','fontsize',20');
xticks(0:60:360);
ylim([0,25]);
xlim([0,360]);
box on;
grid on;


figure;
plot(theta/(2*pi)*360,q,'LineWidth',2);
xlabel('\theta  (deg)','fontname','Times New Roman','fontsize',20');
ylabel('V  (mm/sec)','fontname','Times New Roman','fontsize',20');
xticks(0:60:360);
xlim([0,360]);
box on;
grid on;


figure;
plot(theta/(2*pi)*360,a,'LineWidth',2);
xlabel('\theta  (deg)','fontname','Times New Roman','fontsize',20');
ylabel('a  (m/sec^2)','fontname','Times New Roman','fontsize',20');
xticks(0:60:360);
xlim([0,360]);
box on;
grid on;



figure;
hold on 
axis ([-80 60 -60 80]);
plot([4,4],[-8,0],'LineWidth',2,'Color','k')
plot([-4,-4],[-8,0],'LineWidth',2,'Color','k')
plot([-8,8],[-8,-8],'LineWidth',2,'Color','k')
plot([-8,-4],[-11,-8],'LineWidth',2,'Color','k')
plot([-4,0],[-11,-8],'LineWidth',2,'Color','k')
plot([0,4],[-11,-8],'LineWidth',2,'Color','k')
plot([4,8],[-11,-8],'LineWidth',2,'Color','k')
plot(4.*cos(0 : 0.01 : 2*pi) , 4.*sin( 0: 0.01 : 2*pi),'LineWidth',2,'Color','k');
plot(2.5.*cos(0 : 0.01 : 2*pi) , 2.5.*sin( 0: 0.01 : 2*pi),'LineWidth',2,'Color','k');
plot(40.*cos(0 : 0.01 : 2*pi) , 40.*sin( 0: 0.01 : 2*pi),'LineWidth',2,'LineStyle','--');

plot(X,Y,'LineWidth',2);axis square
xlabel('X','fontname','Times New Roman','fontsize',20');
ylabel('Y','fontname','Times New Roman','fontsize',20');
title('cam profiles','fontname','Times New Roman','fontsize',20');
box on;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Generate contour.scr                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = X';
Y = Y';

fid = fopen('CamProfile_Trans_Flat_Face_Cycloid.scr','w');
fprintf(fid,'spline ');
fprintf(fid,'%f,%f ', [X; Y]);
fclose(fid);
