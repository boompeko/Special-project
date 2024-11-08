
clear;
clc;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%參數
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%       !!!!!!!!!!!!!!!!!!!  (R/N) > E   !!!!!!!!!!!!!!!!!!!!!!

N = 15 ;%Number of rollers
Rr = 18 ;%Radius of the roller
R = 120 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E =5 ;% Eccentricity - offset from input shaft to a cycloidal disk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作圖
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CUT = 10;%切割倍數
tick = 40;
SHOW = sprintf("N = %d, Rr = %d, R = %d, E = %d", N, Rr, R, E);
file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\output'; %家裡電腦
%file_path = 'C:\Users\Johnny Jou\Documents\GitHub\Special-project\Cycloidal Gear\output';  %筆記電腦


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%計算曲率中心和接觸點
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : 1 : ((360*CUT*N)+15)

    phi(i) = i / (180*CUT*N) * pi;

    psi(i) = atan(-sin((1+N)*phi(i))/((R/(E*N))-cos((1+N)*phi(i))));
   
    
    L(i) = (R^2+(E*N)^2-2*R*(E*N)*cos((1+N)*phi(i)))^0.5 + Rr;

    [newRc(i)] = newRrcal_hypocycloid(N,R,E,Rr,phi(i)); %曲率半徑

    X(i) = cos(N*phi(i))*E*(N+1) + cos(phi(i)-psi(i))*(L(i)); %接觸點C
    Y(i) = -sin(N*phi(i))*E*(N+1) + sin(phi(i)-psi(i))*(L(i));

    Kx(i)= cos(N*phi(i))*E*(N+1) + cos(phi(i)-psi(i))*(L(i)+newRc(i)); % 點K
    Ky(i)= -sin(N*phi(i))*E*(N+1) + sin(phi(i)-psi(i))*(L(i)+newRc(i));

    Orx(i) = cos(N*phi(i))*E*(N+1) + cos(phi(i)-psi(i))*(L(i)-Rr);  % 點Or
    Ory(i) = -sin(N*phi(i))*E*(N+1) + sin(phi(i)-psi(i))*(L(i)-Rr);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%繪製擺線輪輪廓
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a = figure('Visible', 'on');
plot(X,Y,'LineWidth',2,'Color','k');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%xlabel('X (mm)','fontname','Times New Roman','fontsize',20');
%ylabel('Y (mm)','fontname','Times New Roman','fontsize',20');
%title('擺線輪輪廓','fontsize',20');
axis equal;
xlim([-160,160]);
ylim([-160,160]);
xticks(-2*R:tick:2*R);
yticks(-2*R:tick:2*R);

box on;
grid on;
axis square;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%繪製曲率半徑
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z(1) = 0;
bRc(1) = newRc(360*CUT*(N-1));

for i = 1 : 1 : ceil(360*CUT)

    z(i+1) = i/CUT;
    bRc(i+1) = newRc(i);

end


b = figure('Visible', 'off');
plot(z,bRc,'LineWidth',2, 'Color','b');
%xlabel('input angle (deg)','fontname','Times New Roman','fontsize',18');
%ylabel('radius of curvature (mm)','fontname','Times New Roman','fontsize',18');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
title('曲率半徑','FontSize',16);
xlim([0,360]);
ylim([-R,R]);
xticks(0:60:(360));
yticks(-R:tick:R);

box on;
grid on;
axis square;


box on;
grid on;
axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%誤差分析
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:1:360*CUT*15/16
    
    OM(i) = atan(Ky(i)/Kx(i));
    
    phiin(i) = -i / (180*CUT) * pi ;
    phi2(i) =  -i / (180*CUT) * pi -OM(i) ;
    
    Ocx(i) = (E)*cos(phiin(i)); %點Oc
    Ocy(i) = (E)*sin(phiin(i));

    b1(i) = (Kx(i)^2 + Ky(i)^2)^0.5;
    b2(i) = newRc(i)+Rr;
    b3(i) = R;
    b4(i) = E;
  
    error = 0.015; %長度誤差
    input_error = 0.1*pi/180; %角度誤差
    
    H(i) = -2*b3(i)*b4(i)*sin(phi2(i));
    
    I(i) = -2*b3(i)*b4(i)*cos(phi2(i)) + 2*b1(i)*b3(i);

    J(i) = b1(i)^2 - b2(i)^2 + b3(i)^2 + b4(i)^2 - 2*b1(i)*b4(i)*cos(phi2(i));

    
    theta(i) = atan2((Ory(i) - Ocy(i)),(Orx(i) - Ocx(i)))-OM(i); % 輸出角 依照點與點之間位置算
    theta1(i) = 2*atan((H(i)-(H(i)^2+I(i)^2-J(i)^2)^0.5)/(I(i)+J(i))); % 輸出角 依據講義直接由HIJ 推
    theta2(i) = 2*atan((H(i)+(H(i)^2+I(i)^2-J(i)^2)^0.5)/(I(i)+J(i)));
    theta3(i) = - phiin(i)/N - OM(i);

    error_1 = 0.022;

    E1(i) = ((-b3(i)*cos(theta(i))+b1(i)-b4(i)*cos(phi2(i))))/(H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*error_1*180/pi;
    
   

    error_r = 0.025;
        
    E2_r(i) = cos(asin(E*(N+1)*sin(phi2(i)-phi(i)+psi(i))/(X(i)^2+Y(i)^2)^0.5))*((-b2(i)))/(H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*error_r*180/pi;
    
    if E2_r(i)> 0.5
    
        E2_r(i) = NaN;

    end

    error_Rr = 0.011;
        
    E2_Rr(i) = ((-b2(i)))/(H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*error_Rr*180/pi;
    
    if E2_Rr(i)> 0.5
    
        E2_Rr(i) = NaN;

    end

    error_3 = 0.022; 

        
    E3(i) = (b4(i)*cos(phi2(i)-theta(i))+b3(i)-b1(i)*cos(theta(i)))/(H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*error*180/pi;
    
    if E3(i)> 0.5
    
        E3(i) = NaN;

    end

    error_4 = 0.008; 

    E4(i) = (b3(i)*cos(phi2(i)-theta(i))+b4(i)-b1(i)*cos(phi2(i)))/(H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*error_4*180/pi;
    
    if E4(i)> 0.5
    
        E4(i) = NaN;

    end

    %Ephi(i) = (a3(i)*a4(i)*cos(phi2(i))*sin(theta(i))-a3(i)*a4(i)*sin(phi2(i))*cos(theta(i))+a1(i)*a4(i)*sin(phi2(i)))/ (H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*input_error*180/pi;
    Ephi(i) = (b3(i)*b4(i)*sin(theta(i)-phi2(i))+b1(i)*b4(i)*sin(phi2(i)))/ (H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*input_error*180/pi;

    max(i) = abs(E3(i))+abs(E2_r(i))+abs(E2_Rr(i))+abs(E4(i));
    rss(i) = ((E3(i))^2+(E2_r(i))^2+(E2_Rr(i))^2+(E4(i))^2)^0.5;



    f(i) = i/CUT ;
    


end



TO1 = figure('Visible', 'on');
hold on
plot(f,E2_r,'LineWidth',2,'Color','r','LineStyle','-');
plot(f,E2_Rr,'LineWidth',2,'Color','r','LineStyle','--');
plot(f,E3,'LineWidth',2,'Color','g','LineStyle','-');
plot(f,E4,'LineWidth',2,'Color','g','LineStyle','--');
plot(f,Ephi,'LineWidth',2,'Color','b','LineStyle','-');
legend("\epsilon_{r}","\epsilon_{Rr}","\epsilon_{3}","\epsilon_{4}","\epsilon_{\phi}",'Location','northeastoutside');
set(legend,'FontSize',20');
hold off
set(gca,'FontSize',14');
xlabel('Input angle (deg)','fontname','Times New Roman','fontsize',20');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',20');
%title('epicycloid reducer error','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.1,0.1]);
xticks(0:60:(360));

box on;
grid on;

TO2 = figure('Visible', 'on');
hold on
plot(f,rss,'LineWidth',2,'Color','r','LineStyle','-');
plot(f,max,'LineWidth',2,'Color','g','LineStyle','--');
legend("\epsilon_{rss}","\epsilon_{max}",'Location','northeastoutside');
set(legend,'FontSize',20');
% lgd = legend;
% lgd.Layout.Tile = 2;
hold off
set(gca,'FontSize',14');
xlabel('Input angle (deg)','fontname','Times New Roman','fontsize',20');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',20');
%title('epicycloid reducer error','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([0,0.5]);
xticks(0:60:(360));


box on;
grid on;


E1_f = figure('Visible', 'off');
hold on
plot(f,E1,'LineWidth',2,'Color','b');
%plot(f,bRc,'LineWidth',2, 'Color','r','LineStyle','--');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',20');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',20');
%title('\epsilon1','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.5,0.5]);
xticks(0:60:(360));
set(gca,'FontSize',16');
hold off
box on;
grid on;
axis square;

E2_r_f = figure('Visible', 'off');
hold on
plot(f,E2_r,'LineWidth',2,'Color','b');
%plot(f,bRc,'LineWidth',2, 'Color','r','LineStyle','--');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',20');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',20');
%title('\epsilon2','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.5,0.5]);
xticks(0:60:(360));
set(gca,'FontSize',16');
hold off
box on;
grid on;
axis square;

E2_Rr_f = figure('Visible', 'off');
hold on
plot(f,E2_Rr,'LineWidth',2,'Color','b');
%plot(f,bRc,'LineWidth',2, 'Color','r','LineStyle','--');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',20');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',20');
%title('\epsilon2','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.5,0.5]);
xticks(0:60:(360));
set(gca,'FontSize',16');
hold off
box on;
grid on;
axis square;

E3_f = figure('Visible', 'off');
hold on
plot(f,E3,'LineWidth',2,'Color','b');
%plot(f,bRc,'LineWidth',2, 'Color','r','LineStyle','--');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',20');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',20');
%title('\epsilon3','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.5,0.5]);
xticks(0:60:(360));
set(gca,'FontSize',16');
hold off
box on;
grid on;
axis square;

E4_f = figure('Visible', 'off');
hold on
plot(f,E4,'LineWidth',2,'Color','b');
%plot(f,bRc,'LineWidth',2, 'Color','r','LineStyle','--');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',20');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',20');
%title('\epsilon4','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.5,0.5]);
xticks(0:60:(360));
set(gca,'FontSize',16');
hold off
box on;
grid on;
axis square;

Ephi_f = figure('Visible', 'off');
hold on
plot(f,Ephi,'LineWidth',2,'Color','b');
set(gca,'FontSize',14');
%plot(f,bRc,'LineWidth',2, 'Color','r','LineStyle','--');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',20');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',20');
%title('\epsiloninput','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.02,0.02]);
xticks(0:60:(360));

hold off
box on;
grid on;
axis square;

j = figure('Visible', 'off');
hold on
%plot(f,abs(theta1),'LineWidth',2,'Color','b');
plot(f,(theta),'LineWidth',2,'Color','g');
%plot(f,abs(theta2),'LineWidth',2,'Color','r');
plot(f,(theta3),'LineWidth',2,'Color','k');

xlim([0,360]);
% % 
a = figure('Visible', 'off');
hold on
%plot(f,(theta1),'LineWidth',2,'Color','b');
plot(f,(theta),'LineWidth',2,'Color','g');
%plot(f,(theta2),'LineWidth',2,'Color','r');
%plot(f,(theta3),'LineWidth',2,'Color','k');
xlim([0,360]);

b= figure('Visible', 'off');
hold on
plot(f,(theta1),'LineWidth',2,'Color','b');
%plot(f,(theta),'LineWidth',2,'Color','g');
%plot(f,(theta2),'LineWidth',2,'Color','r');
%plot(f,(theta3),'LineWidth',2,'Color','k');
xlim([0,360]);

c= figure('Visible', 'off');
hold on
%plot(f,(theta1),'LineWidth',2,'Color','b');
%plot(f,(theta),'LineWidth',2,'Color','g');
plot(f,(theta2),'LineWidth',2,'Color','r');
%plot(f,(theta3),'LineWidth',2,'Color','k');

xlim([0,360]);

k = figure('Visible', 'off');
hold on
plot(Kx,Ky,'LineWidth',2,'Color','b');
xlim([0,360]);

box on;
grid on;
axis square;


% file_name = '誤差分析.png';
% 
% 
% full_file_path = fullfile(file_path, file_name);
% 
% saveas(gcf, full_file_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%生成scr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% scrfilename = fullfile(file_path, 'Cycloidal_Drive_hypocycloid.scr');
% fid = fopen(scrfilename,'w');
% fprintf(fid,'spline ');
% fprintf(fid,'%f,%f\n', [X; Y]);
% 
% fclose(fid);