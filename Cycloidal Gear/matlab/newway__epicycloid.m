
clear;
clc;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%參數
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       !!!!!!!!!!!!!!!!!!!  (R/(N*E)) > 1  !!!!!!!!!!!!!!!!!!!!!!

N = 15 ;%Number of rollers
Rr = 18 ;%Radius of the roller
R = 120 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E =5 ;% Eccentricity - offset from input shaft to a cycloidal disk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作圖
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CUT = 100;%切割倍數
tick = 40;
SHOW = sprintf("N = %d, Rr = %d, R = %d, E = %d", N, Rr, R, E);
file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\output'; %家裡電腦
%file_path = 'C:\Users\Johnny Jou\Documents\GitHub\Special-project\Cycloidal Gear\output';  %筆記電腦


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%計算曲率中心和接觸點
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : 1 : ((360*CUT*(N-1))+14)

    phi(i) = i / (180*CUT*(N-1)) * pi;  % phi = phi3

    psi(i) = atan(sin((1-N)*phi(i))/((R/(E*N))-cos((1-N)*phi(i))));
   
    
    L(i) = (R^2+(E*N)^2-2*R*(E*N)*cos((1-N)*phi(i)))^0.5 - Rr;

    X(i) = cos(N*phi(i))*E*(N-1) + cos(phi(i)+psi(i))*(L(i)); %輪廓
    Y(i) = -sin(N*phi(i))*E*(N-1) - sin(phi(i)+psi(i))*(L(i));

    [newRc(i)] = newRrcal_epicycloid(N,R,E,Rr,phi(i));  %曲率半徑

    Kx(i) = R-(Rr-newRc(i))*cos(psi(i));  %曲率中心
    Ky(i) = (Rr-newRc(i))*sin(psi(i));
    k_test(i) = Rr-newRc(i);
    
    er(i) = (X(i)^2+Y(i)^2)^0.5;

    tphi(i) = (asin((E*(N-1))/er(i)*sin(pi-psi(i) - phi(i)*(1-N))))*180/pi;

%     Xc2(i) = (R)*cos((N-1)*t(i))-(Rr+newRc(i))*cos(phi(i)+(N-1)*t(i));
%     Yc2(i) = (R)*sin((N-1)*t(i))-(Rr+newRc(i))*sin(phi(i)+(N-1)*t(i));


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%繪製擺線輪輪廓
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a = figure('Visible', 'on');
plot(X,Y,'LineWidth',2,'Color','b');
xlabel('X (mm)','fontname','Times New Roman','fontsize',18');
ylabel('Y (mm)','fontname','Times New Roman','fontsize',18');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
%title('擺線輪輪廓',SHOW);
title('擺線輪輪廓','fontsize',20');
axis equal;
xlim([-(R),(R)]);
ylim([-(R),(R)]);
xticks(-2*R:tick:2*R);
yticks(-2*R:tick:2*R);

box on;
grid on;
axis square;

% file_name = '外擺線輪.png';


% full_file_path = fullfile(file_path, file_name);
% 
% saveas(gcf, full_file_path);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%繪製曲率半徑
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z(1) = 0;
bRc(1) = newRc(360*CUT*(N-1));
cRc(1) = 1/newRc(360*CUT*(N-1));

for i = 1 : 1 : ceil(360*CUT)

    z(i+1) = i/CUT;
    bRc(i+1) = newRc(i);
    cRc(i+1) = 1/newRc(i);
end


b = figure('Visible', 'on');
plot(z,bRc,'LineWidth',2, 'Color','b');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',18');
ylabel('radius of curvature (mm)','fontname','Times New Roman','fontsize',18');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
title('曲率半徑','FontSize',16);
xlim([0,360]);
ylim([-R,R]);
xticks(0:60:(360));
yticks(-R:tick:R);

box on;
grid on;
axis square;

c = figure('Visible', 'on');
plot(z,cRc,'LineWidth',2, 'Color','b');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',18');
ylabel('radius of curvature (mm)','fontname','Times New Roman','fontsize',18');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
title('曲率半徑','FontSize',16);
xlim([0,360]);
ylim([-0.1,0.1]);
xticks(0:60:(360));
yticks(-R:tick:R);

box on;
grid on;
axis square;


box on;
grid on;
axis square;

% file_name = 'radius_of_curvature_epicycloid.png';
% 
% 
% full_file_path = fullfile(file_path, file_name);
% 
% saveas(gcf, full_file_path);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%誤差分析
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for i=1:1:(360*CUT+1)
    
    phi2(i) =  -i / (180*CUT) * pi   ; % 輸入角phi2
    
 
    
    Ocx(i) = (E)*cos(phi2(i)); %點Oc
    Ocy(i) = (E)*sin(phi2(i));


    a1(i) = R;
    a2(i) = -newRc(i)+Rr;
    a3(i) = ((Kx(i) - Ocx(i))^2+(Ky(i) - Ocy(i))^2)^0.5;
    a4(i) = E;
    
    
    error = 0.015; %長度誤差
    input_error = 0.1*pi/180; %角度誤差
  
    H(i) = -2*a3(i)*a4(i)*sin(phi2(i));
    
    I(i) = -2*a3(i)*a4(i)*cos(phi2(i)) + 2*a1(i)*a3(i);

    J(i) = a1(i)^2 - a2(i)^2 + a3(i)^2 + a4(i)^2 - 2*a1(i)*a4(i)*cos(phi2(i));

    theta1(i) = 2*atan((H(i)-(H(i)^2+I(i)^2-J(i)^2)^0.5)/(I(i)+J(i))); % 輸出角 依據講義直接由HIJ 推
    theta(i) = atan2((Ky(i) - Ocy(i)),(Kx(i) - Ocx(i))); % 輸出角 依照點與點之間位置算
    theta2(i) = 2*atan((H(i)+(H(i)^2+I(i)^2-J(i)^2)^0.5)/(I(i)+J(i)));

    E1(i) = ((-a3(i)*cos(theta(i))+a1(i)-a4(i)*cos(phi2(i))))/(H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*error*180/pi;

    E2(i) = ((-a2(i)))/(H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*error*180/pi;
    
    E3(i) = (a4(i)*cos(phi2(i)-theta(i))+a3(i)-a1(i)*cos(theta(i)))/(H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*error*180/pi;

    E4(i) = (a3(i)*cos(phi2(i)-theta(i))+a4(i)-a1(i)*cos(phi2(i)))/(H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*error*180/pi;
    
    %Ephi(i) = (a3(i)*a4(i)*cos(phi2(i))*sin(theta(i))-a3(i)*a4(i)*sin(phi2(i))*cos(theta(i))+a1(i)*a4(i)*sin(phi2(i)))/ (H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*input_error*180/pi;
    Ephi(i) = (a3(i)*a4(i)*sin(theta(i)-phi2(i))+a1(i)*a4(i)*sin(phi2(i)))/ (H(i)*cos(theta(i))-I(i)*sin(theta(i)))*2*input_error*180/pi;

    max(i) = abs(E1(i))+abs(E2(i))+abs(E3(i))+abs(E4(i));
    rss(i) = ((E1(i))^2+(E2(i))^2+(E3(i))^2+(E4(i))^2)^0.5;


    f(i) = i/CUT ;

end


g = figure('Visible', 'on');
hold on
plot(f,E1,'LineWidth',2,'Color','b');
plot(f,E2,'LineWidth',2,'Color','r');
plot(f,E3,'LineWidth',2,'Color','g');
plot(f,E4,'LineWidth',2,'Color','m');
plot(f,Ephi,'LineWidth',2,'Color','y');
plot(f,rss,'LineWidth',2,'Color','c');
plot(f,max,'LineWidth',2,'Color','k');
legend("\epsilon1","\epsilon2","\epsilon3","\epsilon4","\epsilon\phi2","\epsilonrms","\epsilonmax")
hold off
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',18');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',18');
set(gca, 'Fontname', 'Times New Roman','FontSize',14);
title('epicycloid reducer error','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.1,0.1]);
xticks(0:60:(360));
box on;
grid on;
axis square;

E1_f = figure('Visible', 'on');
hold on
plot(f,E1,'LineWidth',2,'Color','b');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',18');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',18');
title('\epsilon1','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.5,0.5]);
xticks(0:60:(360));
hold off
box on;
grid on;
axis square;

E2_f = figure('Visible', 'on');
hold on
plot(f,E2,'LineWidth',2,'Color','b');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',18');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',18');
title('\epsilon2','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.5,0.5]);
xticks(0:60:(360));
hold off
box on;
grid on;
axis square;

E3_f = figure('Visible', 'on');
hold on
plot(f,E3,'LineWidth',2,'Color','b');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',18');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',18');
title('\epsilon3','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.5,0.5]);
xticks(0:60:(360));
hold off
box on;
grid on;
axis square;

E4_f = figure('Visible', 'on');
hold on
plot(f,E4,'LineWidth',2,'Color','b');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',18');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',18');
title('\epsilon4','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.5,0.5]);
xticks(0:60:(360));
hold off
box on;
grid on;
axis square;

Ephi_f = figure('Visible', 'on');
hold on
plot(f,Ephi,'LineWidth',2,'Color','b');
xlabel('input angle (deg)','fontname','Times New Roman','fontsize',18');
ylabel('Errors (deg)','fontname','Times New Roman','fontsize',18');
title('\epsiloninput','fontname','標楷體','FontSize',16);
xlim([0,360]);
ylim([-0.02,0.02]);
xticks(0:60:(360));
hold off
box on;
grid on;
axis square;

j = figure('Visible', 'off');
hold on
plot(f,theta1,'LineWidth',2,'Color','b');
plot(f,theta,'LineWidth',2,'Color','g');
plot(f,theta2,'LineWidth',2,'Color','r');
xlim([0,360]);
% % 


k = figure('Visible', 'off');
hold on
plot(Kx,Ky,'LineWidth',2,'Color','b');
xlim([0,360]);

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
%     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%迴圈生成圖形幀並保存為GIF         %  Stationary ring gear type epicycloid reducer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% for i=1:1:ceil(360*CUT/(N-1))
%     
%     t(i) = i / (180*CUT) * pi ;
%     x = (E)*cos((1-N)*t(i));
%     y = (E)*sin((1-N)*t(i));
% 
%     c = figure('Visible', 'off');
%     hold on
%     plot([Xc1(i),R],[Yc1(i),0],'LineWidth',2);
%     plot([0,x],[0,y],'LineWidth',2);
%     plot([x,Xc1(i)],[y,Yc1(i)],'LineWidth',2);
% 
%     %覆蓋多於線段
% 
%     fill(1.5.*cos(0 : 0.01 : 2*pi)+x , 1.5.*sin( 0: 0.01 : 2*pi)+y,'w-');
%     fill(1.5.*cos(0 : 0.01 : 2*pi)+Xc1(i) , 1.5.*sin( 0: 0.01 : 2*pi)+Yc1(i),'w-');
%     fill(2.5.*cos(0 : 0.01 : 2*pi) , 2.5.*sin( 0: 0.01 : 2*pi),'w-');
%     fill(2.5.*cos(0 : 0.01 : 2*pi)+R , 2.5.*sin( 0: 0.01 : 2*pi),'w-');
% 
% 
%     %非地桿接點
%     plot(1.5.*cos(0 : 0.01 : 2*pi)+x , 1.5.*sin( 0: 0.01 : 2*pi)+y,'LineWidth',2,'Color','k');
%     plot(1.5.*cos(0 : 0.01 : 2*pi)+Xc1(i) , 1.5.*sin( 0: 0.01 : 2*pi)+Yc1(i),'LineWidth',2,'Color','k');
%     
% 
% 
%     
%     %原點地桿
% 
%     offsetX = 0;
%     offsetY = 0;
%     
%     
%     plot([4+offsetX,4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
%     plot([-4+offsetX,-4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
%     plot([-8+offsetX,8+offsetX],[-8+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([-8+offsetX,-4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([-4+offsetX,0+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([0+offsetX,4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([4+offsetX,8+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot(4.*cos(0 : 0.01 : 2*pi)+offsetX , 4.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
%     plot(2.5.*cos(0 : 0.01 : 2*pi)+offsetX , 2.5.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
% 
% 
%     %針輪地桿
% 
%     offsetX = R;
%     offsetY = 0;
%     
% 
%     plot([4+offsetX,4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
%     plot([-4+offsetX,-4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
%     plot([-8+offsetX,8+offsetX],[-8+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([-8+offsetX,-4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([-4+offsetX,0+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([0+offsetX,4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([4+offsetX,8+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot(4.*cos(0 : 0.01 : 2*pi)+offsetX , 4.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
%     plot(2.5.*cos(0 : 0.01 : 2*pi)+offsetX , 2.5.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
%     
%     box on;
%     grid on;
%     axis equal;
%     xlabel('X (mm)','fontname','Times New Roman','fontsize',20');
%     ylabel('Y (mm)','fontname','Times New Roman','fontsize',20');
%     title(SHOW,'等效連桿');
%     
%     xlim([-2*E,1.5*R]);
%     ylim([-(0.8*R)/2,(0.8*R)/2]);
%     xticks(-2*R:tick:2*R);
%     yticks(-2*R:tick:2*R);
%     
%     hold off
%     
%     % 擷取圖形幀
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
% 
%     giffilename = fullfile(file_path, 'Newway_Stationary_ring_move_pic_epicycloid.gif');
%     if i==1
%         imwrite(I,map,giffilename,'gif','writeMode','overwrite','LoopCount',inf,'delaytime',0.0001/CUT,'loopcount',inf);
%     else
%         imwrite(I,map,giffilename,'gif','writeMode','append','delaytime',0.0001/CUT);
%     end
% 
%     fprintf("Stationary ring epicycloid countdown %d \n",(ceil(360*CUT/(N-1))-i))
%     
%     
% 
%     close;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%迴圈生成圖形幀並保存為GIF         %  Rotating ring gear type epicycloid reducer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% for i=1:1:ceil(360*CUT/(N-1))
%     
%     t(i) = i / (180*CUT) * pi ;
%     x = (R)*cos((N-1)*t(i));
%     y = (R)*sin((N-1)*t(i));
% 
%     c = figure('Visible', 'off');
%     hold on
%     plot([Xc2(i),E],[Yc2(i),0],'LineWidth',2);
%     plot([0,x],[0,y],'LineWidth',2);
%     plot([x,Xc2(i)],[y,Yc2(i)],'LineWidth',2);
% 
%     %覆蓋多於線段
% 
%     fill(1.5.*cos(0 : 0.01 : 2*pi)+x , 1.5.*sin( 0: 0.01 : 2*pi)+y,'w-');
%     fill(1.5.*cos(0 : 0.01 : 2*pi)+Xc2(i) , 1.5.*sin( 0: 0.01 : 2*pi)+Yc2(i),'w-');
%     fill(2.5.*cos(0 : 0.01 : 2*pi) , 2.5.*sin( 0: 0.01 : 2*pi),'w-');
%     fill(2.5.*cos(0 : 0.01 : 2*pi)+E , 2.5.*sin( 0: 0.01 : 2*pi),'w-');
% 
% 
%     %非地桿接點
%     plot(1.5.*cos(0 : 0.01 : 2*pi)+x , 1.5.*sin( 0: 0.01 : 2*pi)+y,'LineWidth',2,'Color','k');
%     plot(1.5.*cos(0 : 0.01 : 2*pi)+Xc2(i) , 1.5.*sin( 0: 0.01 : 2*pi)+Yc2(i),'LineWidth',2,'Color','k');
%     
% 
% 
%     
%     %原點地桿
% 
%     offsetX = 0;
%     offsetY = 0;
%     
%     
%     plot([4+offsetX,4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
%     plot([-4+offsetX,-4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
%     plot([-8+offsetX,8+offsetX],[-8+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([-8+offsetX,-4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([-4+offsetX,0+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([0+offsetX,4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([4+offsetX,8+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot(4.*cos(0 : 0.01 : 2*pi)+offsetX , 4.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
%     plot(2.5.*cos(0 : 0.01 : 2*pi)+offsetX , 2.5.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
% 
% 
%     %針輪地桿
% 
%     offsetX = E;
%     offsetY = 0;
%     
% 
%     plot([4+offsetX,4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
%     plot([-4+offsetX,-4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
%     plot([-8+offsetX,8+offsetX],[-8+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([-8+offsetX,-4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([-4+offsetX,0+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([0+offsetX,4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot([4+offsetX,8+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
%     plot(4.*cos(0 : 0.01 : 2*pi)+offsetX , 4.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
%     plot(2.5.*cos(0 : 0.01 : 2*pi)+offsetX , 2.5.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
%     
%     box on;
%     grid on;
%     axis equal;
%     xlabel('X (mm)','fontname','Times New Roman','fontsize',20');
%     ylabel('Y (mm)','fontname','Times New Roman','fontsize',20');
%     title(SHOW,'等效連桿');
%     
%     xlim([-2*E,1.5*R]);
%     ylim([-(0.8*R)/2,(0.8*R)/2]);
%     xticks(-2*R:tick:2*R);
%     yticks(-2*R:tick:2*R);
%     
%     hold off
%     
%     % 擷取圖形幀
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
% 
%     giffilename = fullfile(file_path, 'Newway_Stationary_ring_move_pic_epicycloid.gif');
%     if i==1
%         imwrite(I,map,giffilename,'gif','writeMode','overwrite','LoopCount',inf,'delaytime',0.0001/CUT,'loopcount',inf);
%     else
%         imwrite(I,map,giffilename,'gif','writeMode','append','delaytime',0.0001/CUT);
%     end
% 
%     fprintf("Rotating ring epicycloid countdown %d \n",(ceil(360*CUT/(N-1))-i))
%     
%     
% 
%     close;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%生成scr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% scrfilename = fullfile(file_path, 'Cycloidal_Drive_epicycloid.scr');
% fid = fopen(scrfilename,'w');
% fprintf(fid,'spline ');
% fprintf(fid,'%f,%f\n', [X; Y]);
% fclose(fid);

[ans] = newRrcal_epicycloid(15,120,5,18,0.057560592920962710041210399673131);

1/ans
