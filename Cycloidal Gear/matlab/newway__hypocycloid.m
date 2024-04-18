
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
tick = 20;
SHOW = sprintf("N = %d, Rr = %d, R = %d, E = %d", N, Rr, R, E);
file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\output'; %家裡電腦
%file_path = 'C:\Users\Johnny Jou\Documents\GitHub\Special-project\Cycloidal Gear\output';  %筆記電腦


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%計算曲率中心和接觸點
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : 1 : ((360*CUT)+1)

    t(i) = i / (180*CUT) * pi;

    phi(i) = atan(-sin((1+N)*t(i))/((R/(E*N))-cos((1+N)*t(i))));
   
    
    L(i) = (R^2+(E*N)^2-2*R*(E*N)*cos((1+N)*t(i)))^0.5 + Rr;

    [newRc(i)] = newRrcal_hypocycloid(N,R,E,Rr,t(i));

    X(i) = cos(N*t(i))*E*(N+1) + cos(t(i)-phi(i))*(L(i));
    Y(i) = -sin(N*t(i))*E*(N+1) + sin(t(i)-phi(i))*(L(i));

    Kx(i)= E*cos(N*t(i))*(N + 1)+(L(i)+ Rr)*cos(t(i) - phi(i));
    Ky(i)= -E*sin(N*t(i))*(N + 1)+(L(i)+ Rr)*sin(t(i) - phi(i));

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



for i=1:1:ceil(360*CUT/(N-1))
    
    t(i) = i / (180*CUT) * pi ;
    x = (E)*cos((1-N)*t(i));
    y = (E)*sin((1-N)*t(i));

    zeta = atan((Xc1(i) - x)/(Yc1(i) - y));
    
    
    a1 = (Kx(i))^2 + (Ky(i))^2;
    a2 = newRc(i)+Rr;
    a3 = R;
    a4 = E;

    error = 0.01;
    
    H(i) = -2*a4*a3*cos((1-N)*t(i));
    
    I(i) = -2*a4*a3*cos((1-N)*t(i)) + 2*a1*a3;

    J(i) = a1^2 - a2^2 + a3^2 + a4^2 +2*a1*a4*cos((1-N)*t(i));


    E1(i) = ((-a3*cos(zeta)+a1-a4*cos((1-N)*t(i)))*2*error)/(H(i)*cos(zeta)-I(i)*sin(zeta));

    E2(i) = ((a2)*2*error)/(H(i)*cos(zeta)-I(i)*sin(zeta));
    
    E3(i) = ((-a4*cos((1-N)*t(i)-zeta)+a3)*2*error)/(H(i)*cos(zeta)-I(i)*sin(zeta));

    E4(i) = ((-a3*cos((1-N)*t(i)-zeta)+a4-a1*cos((1-N)*t(i)))*2*error)/(H(i)*cos(zeta)-I(i)*sin(zeta));
    

    f(i) = i/CUT ;

end


g = figure('Visible', 'on');
hold on
plot(f,E1,'LineWidth',2,'Color','b');
plot(f,E2,'LineWidth',2,'Color','r');
plot(f,E3,'LineWidth',2,'Color','g');
plot(f,E4,'LineWidth',2,'Color','m');
hold off
xlabel('cycloidal disk rotation angle (θ)','fontname','Times New Roman','fontsize',18');
ylabel('Errors (θ)','fontname','Times New Roman','fontsize',18');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%title(SHOW,'Error due only to \Deltaa1 ','fontname','標楷體','FontSize',16);
xlim([0,ceil(360/(N-1))]);
%ylim([-0.5/1000,0.5/1000]);
% xticks(0:5:(360/(N-1)));
%yticks(-0.5/1000:10:0.5/1000);




box on;
grid on;
axis square;
% 
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