
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


CUT = 50;%切割倍數
tick = 20;
SHOW = sprintf("N = %d, Rr = %d, R = %d, E = %d", N, Rr, R, E);
file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\output'; %家裡電腦
%file_path = 'C:\Users\Johnny Jou\Documents\GitHub\Special-project\Cycloidal Gear\output';  %筆記電腦


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%計算曲率中心和接觸點
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : 1 : ((360*CUT*N)+1)

    phi(i) = i / (180*CUT) * pi;

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
plot(X,Y,'LineWidth',2);
xlabel('X','fontname','Times New Roman','fontsize',20');
ylabel('Y','fontname','Times New Roman','fontsize',20');
title('擺線輪輪廓',SHOW);
axis equal;
xlim([-(R+2*Rr),(R+2*Rr)]);
ylim([-(R+2*Rr),(R+2*Rr)]);
xticks(-2*R:tick:2*R);
yticks(-2*R:tick:2*R);

box on;
grid on;
axis square;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%誤差分析
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i=1:1:ceil(360*CUT/(N+1))
    
    OM(i) = atan(Ky(i)-Kx(i));
    
    phi(i) = i / (180*CUT) * pi - OM(i) ;
    
    
    b1(i) = (Kx(i))^2 + (Ky(i))^2;
    b2(i) = newRc(i)+Rr;
    b3(i) = R;
    b4(i) = E;
  
    
    error = 0.015;
    
    H(i) = -2*b3(i)*b4(i)*sin(phi(i));
    
    I(i) = -2*b3(i)*b4(i)*cos(phi(i)) + 2*b1(i)*b3(i);

    J(i) = b1(i)^2 - b2(i)^2 + b3(i)^2 + b4(i)^2 - 2*b1(i)*b4(i)*cos(phi(i));
    
    zeta(i) = atan((Ory(i) - y(i))/(Orx(i) - x(i))) - OM(i) ;
    
%     zeta1(i) = t(i)/N -OM(i) ;
%     zeta2(i) = 2*atan((H(i)-(H(i)^2+I(i)^2-J(i)^2)^0.5)/(I(i)+J(i))) -OM(i);


    E1(i) = ((-b3(i)*cos(zeta(i))+b1(i)-b4(i)*cos(phi(i))))/(H(i)*cos(zeta(i))-I(i)*sin(zeta(i)))*2*error*180/pi;

    E2(i) = ((b2(i)))/(H(i)*cos(zeta(i))-I(i)*sin(zeta(i)))*2*error*180/pi;
    
    E3(i) = ((b4(i)*cos(phi(i)-zeta(i))+b3(i)-b1(i)*cos(phi(i))))/(H(i)*cos(zeta(i))-I(i)*sin(zeta(i)))*2*error*180/pi;

    E4(i) = ((b3(i)*cos(phi(i)-zeta(i))+b4(i)-b1(i)*cos(phi(i)))/(H(i)*cos(zeta(i))-I(i)*sin(zeta(i)))*2*error)*180/pi;
    
    max(i) = abs(E1(i))+abs(E3(i))+abs(E2(i))+abs(E4(i));
    rms(i) = ((E1(i))^2+(E2(i))^2+(E3(i))^2+(E4(i))^2)^0.5;

    f(i) = i/CUT ;

end


g = figure('Visible', 'on');
hold on
plot(f,E1,'LineWidth',2,'Color','b');
plot(f,E2,'LineWidth',2,'Color','r');
plot(f,E3,'LineWidth',2,'Color','g');
plot(f,E4,'LineWidth',2,'Color','m');
plot(f,rms,'LineWidth',2,'Color','c');
plot(f,max,'LineWidth',2,'Color','k');
legend("\epsilon1","\epsilon2","\epsilon3","\epsilon4","\epsilonrms","\epsilonmax")
hold off
xlabel('cycloidal disk rotation angle (θ)','fontname','Times New Roman','fontsize',18');
ylabel('Errors (θ)','fontname','Times New Roman','fontsize',18');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
title('hypocycloid reducer error','fontname','標楷體','FontSize',16);
xlim([0,ceil(360/(N+1))-1]);
%ylim([-0.5/1000,0.5/1000]);
% xticks(0:5:(360/(N-1)));
%yticks(-0.5/1000:10:0.5/1000);

box on;
grid on;
axis square;
% 

% j = figure('Visible', 'on');
% hold on
% plot(f,zeta1,'LineWidth',2,'Color','b');
% plot(f,zeta,'LineWidth',2,'Color','g');
% plot(f,zeta2,'LineWidth',2,'Color','k');
% hold off

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