
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              parameter input                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rb = 40; 
phi = 25 / 180 * pi; 
beta = 120 / 180 * pi;
beta2 = 120 / 180 * pi;
theta_D = 160 / 180 * pi;
rf = 8;
l=52;
f=80;



for i = 1 : 1 : 360
    theta(i)=i / 180 * pi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Rise：0~120                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:120
    S(i)=phi*(theta(i)/beta-1/(2*pi)*sin(2*pi*theta(i)/beta));
    xi(i)=acos((l^2+f^2-(rb+rf)^2)/(2*l*f))+S(i);
    V(i)=phi/beta-phi/beta*cos(2*pi*theta(i)/beta);
    a(i)=2*pi*phi/(beta^2)*sin(2*pi*theta(i)/beta);
    q(i)=f*abs(V(i))/(1-abs(V(i)));
    QC(i)=(l^2+(f+q(i))^2-2*l*(f+q(i))*cos(xi(i)))^0.5;
    alpha(i)=asin((l*sin(xi(i))/QC(i)));
    X(i)=q(i)*cos(theta(i)+pi)+(QC(i)-rf)*cos(theta(i)+alpha(i));
    Y(i)=q(i)*sin(theta(i)+pi)+(QC(i)-rf)*sin(theta(i)+alpha(i));
    P(i)=((X(i)+rf*cos(theta(i)+alpha(i)))^2+(Y(i)+rf*sin(theta(i)+alpha(i)))^2)^0.5;
    P_(i)=(X(i)^2+Y(i)^2)^0.5;
    X_(i)=P_(i)*cos(theta(i));
    Y_(i)=P_(i)*sin(theta(i));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Fall：160~280                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=161:280
    S(i)=phi*(1-(theta(i)-theta_D)/(beta2)+1/(2*pi)*sin(2*pi*(theta(i)-theta_D)/(beta2)));
    xi(i)=acos((l^2+f^2-(rb+rf)^2)/(2*l*f))+S(i);
    V(i)=-((phi/(beta2)-phi/(beta2)*cos(2*pi*(theta(i)-theta_D)/(beta2))));
    q(i)=f*abs(V(i))/(1-abs(V(i)));
    a(i)=-2*pi*phi/((beta2)^2)*sin(2*pi*(theta(i)-theta_D)/(beta2));
    QC(i)=(l^2+(f+q(i))^2-2*l*(f+q(i))*cos(xi(i)))^0.5;
    alpha(i)=asin((l*sin(xi(i))/QC(i)));
    X(i)=q(i)*cos(theta(i)+pi)+(QC(i)-rf)*cos(theta(i)+alpha(i));
    Y(i)=q(i)*sin(theta(i)+pi)+(QC(i)-rf)*sin(theta(i)+alpha(i));
    P(i)=((X(i)+rf*cos(theta(i)+alpha(i)))^2+(Y(i)+rf*sin(theta(i)+alpha(i)))^2)^0.5;
    P_(i)=(X(i)^2+Y(i)^2)^0.5;
    X_(i)=P_(i)*cos(theta(i));
    Y_(i)=P_(i)*sin(theta(i));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           High Dwell：120~160                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=121:160
    S(i)=phi;
    xi(i)=acos((l^2+f^2-(rb+rf)^2)/(2*l*f))+S(i);
    V(i)=0;
    q(i)=0;
    a(i)=0;
    QC(i)=(l^2+(f+q(i))^2-2*l*(f+q(i))*cos(xi(i)))^0.5;
    alpha(i)=asin((l*sin(xi(i))/QC(i)));
    X(i)=q(i)*cos(theta(i)+pi)+(QC(i)-rf)*cos(theta(i)+alpha(i));
    Y(i)=q(i)*sin(theta(i)+pi)+(QC(i)-rf)*sin(theta(i)+alpha(i));
    P(i)=((X(i)+rf*cos(theta(i)+alpha(i)))^2+(Y(i)+rf*sin(theta(i)+alpha(i)))^2)^0.5;
    P_(i)=(X(i)^2+Y(i)^2)^0.5;
    X_(i)=P_(i)*cos(theta(i));
    Y_(i)=P_(i)*sin(theta(i));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Low Dwell： 280~360                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=281:360
    S(i)=0;
    xi(i)=acos((l^2+f^2-(rb+rf)^2)/(2*l*f))+S(i);
    V(i)=0;
    q(i)=0;
    a(i)=0;
    QC(i)=(l^2+(f+q(i))^2-2*l*(f+q(i))*cos(xi(i)))^0.5;
    alpha(i)=asin((l*sin(xi(i))/QC(i)));
    X(i)=q(i)*cos(theta(i)+pi)+(QC(i)-rf)*cos(theta(i)+alpha(i));
    Y(i)=q(i)*sin(theta(i)+pi)+(QC(i)-rf)*sin(theta(i)+alpha(i));
    P(i)=((X(i)+rf*cos(theta(i)+alpha(i)))^2+(Y(i)+rf*sin(theta(i)+alpha(i)))^2)^0.5;
    P_(i)=(X(i)^2+Y(i)^2)^0.5;
    X_(i)=P_(i)*cos(theta(i));
    Y_(i)=P_(i)*sin(theta(i));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            close 361                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


X(361) = X(1);
Y(361) = Y(1);
X_(361) = X_(1);
Y_(361) = Y_(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Generate contour figure                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
axis equal;
plot(theta/(2*pi)*360,P_,'LineWidth',2);
xlabel('\theta  (deg)','fontname','Times New Roman','fontsize',20');
ylabel('P  (mm)','fontname','Times New Roman','fontsize',20');
xticks(0:60:360);
ylim([35,60]);
xlim([0,360]);
box on;
grid on;


figure;
axis equal;
plot(theta/(2*pi)*360,P,'LineWidth',2);
xlabel('\theta  (deg)','fontname','Times New Roman','fontsize',20');
ylabel('P  (mm)','fontname','Times New Roman','fontsize',20');
xticks(0:60:360);
ylim([35,60]);
xlim([0,360]);
box on;
grid on;


figure;
axis equal;
plot(theta/(2*pi)*360,S,'LineWidth',2);
xlabel('\theta  (deg)','fontname','Times New Roman','fontsize',20');
ylabel('S  (\phi)(rad)','fontname','Times New Roman','fontsize',20');
xticks(0:60:360);
ylim([-0.2,1]);
xlim([0,360]);
box on;
grid on;


figure;
plot(theta/(2*pi)*360,V,'LineWidth',2);
xlabel('\theta  (deg)','fontname','Times New Roman','fontsize',20');
ylabel('V  (rad/sec)','fontname','Times New Roman','fontsize',20');
xticks(0:60:360);
xlim([0,360]);
ylim([-0.8,0.9]);
box on;
grid on;


figure;
plot(theta/(2*pi)*360,a,'LineWidth',2);
xlabel('\theta  (deg)','fontname','Times New Roman','fontsize',20');
ylabel('a  (rad/sec^2)','fontname','Times New Roman','fontsize',20');
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
xlim([-100,50]);
ylim([-80,80]);
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

plot(X_,Y_,'LineWidth',2);axis square
xlabel('X_','fontname','Times New Roman','fontsize',20');
ylabel('Y_','fontname','Times New Roman','fontsize',20');
title('cam profiles','fontname','Times New Roman','fontsize',20');
xlim([-100,50]);
ylim([-80,80]);
box on;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Generate contour.scr                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fid = fopen('CamProfile_rocker_roller_Cycloid.scr','w');
fprintf(fid,'spline ');

for i = 1: 361
    
    fprintf(fid,'%f,%f ', [X(i),Y(i)]);

end
fclose(fid);

fid = fopen('CamProfile_rocker_roller_Cycloid_t.scr','w');
fprintf(fid,'spline ');

for i = 1: 361
    
    fprintf(fid,'%f,%f ', [X_(i),Y_(i)]);

end
fclose(fid);



