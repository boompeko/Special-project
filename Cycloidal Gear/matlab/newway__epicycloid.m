
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

    phi(i) = atan(sin((1-N)*t(i))/((R/(E*N))-cos((1-N)*t(i))));
   
    
    L(i) = (R^2+(E*N)^2-2*R*(E*N)*cos((1-N)*t(i)))^0.5 - Rr;

    X(i) = cos(N*t(i))*E*(N-1) + cos(t(i)+phi(i))*(L(i));
    Y(i) = -sin(N*t(i))*E*(N-1) - sin(t(i)+phi(i))*(L(i));

    [newRc(i)] = newRrcal_epicycloid(N,R,E,Rr,t(i));

    Xc1(i) = R-(Rr-newRc(i))*cos(phi(i));
    Yc1(i) = (Rr-newRc(i))*sin(phi(i));

    Xc2(i) = (R)*cos((N-1)*t(i))-(Rr+newRc(i))*cos(phi(i)+(N-1)*t(i));
    Yc2(i) = (R)*sin((N-1)*t(i))-(Rr+newRc(i))*sin(phi(i)+(N-1)*t(i));


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
scrfilename = fullfile(file_path, 'Cycloidal_Drive_epicycloid.scr');
fid = fopen(scrfilename,'w');
fprintf(fid,'spline ');
fprintf(fid,'%f,%f\n', [X; Y]);
fclose(fid);