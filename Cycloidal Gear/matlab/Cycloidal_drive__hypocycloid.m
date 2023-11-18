
clear;
clc;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%參數
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       !!!!!!!!!!!!!!!!!!!  (R/N) > E   !!!!!!!!!!!!!!!!!!!!!!

N = 24 ;%Number of rollers
Rr = 9 ;%Radius of the roller
R = 150 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E =3 ;% Eccentricity - offset from input shaft to a cycloidal disk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作圖
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CUT = 10;%切割倍數
tick = 20;
SHOW = sprintf("N = %d, Rr = %d, R = %d, E = %d", N, Rr, R, E);
% file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\movie\pic'; 
file_path = 'C:\Users\Johnny Jou\Documents\GitHub\Special-project\Cycloidal Gear\output';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%計算曲率中心和接觸點
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : 1 : ((360*CUT)+1)

    t(i) = i / (180*CUT) * pi;

    phi(i) = atan(sin((1+N)*t(i))/((R/(E*N))-cos((1+N)*t(i))));
    
    Cf = [R+Rr*cos(phi(i));Rr*sin(phi(i));1];
    
    [Rc(i)] = Rrcal_hypocycloid(N,R,E,Rr,t(i));

  
    Ccf = [R+(Rr-Rc(i))*cos(phi(i));(Rr-Rc(i))*sin(phi(i));1];

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%繪製擺線輪輪廓
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a = figure('Visible', 'on');
plot(x,y,'LineWidth',2);
xlabel('X','fontname','Times New Roman','fontsize',20');
ylabel('Y','fontname','Times New Roman','fontsize',20');
title('擺線輪輪廓',SHOW);
axis equal;
xlim([-(1.3*R),(1.3*R)]);
ylim([-(1.3*R),(1.3*R)]);
xticks(-2*R:tick:2*R);
yticks(-2*R:tick:2*R);

box on;
grid on;
axis square;

file_name = 'contour_hypocycloid.png';


full_file_path = fullfile(file_path, file_name);

saveas(gcf, full_file_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%繪製曲率半徑
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z(1) = 0;
bRc(1) = Rc(360*CUT);

for i = 1 : 1 : ceil(360*CUT/(N+1))

    z(i+1) = i/CUT ;
    bRc(i+1) = Rc(i);

end


b = figure('Visible', 'on');
plot(z,bRc,'LineWidth',2);
xlabel('cycloidal disk rotation angle(θ)','fontname','Times New Roman','fontsize',20');
ylabel('radius of curvature','fontname','Times New Roman','fontsize',20');
title('曲率半徑',SHOW);
xlim([0,ceil(360/(N+1))-1]);
ylim([-1.3*R,1.3*R]);
xticks(0:5:(360/(N-1)));
yticks(-2*R:tick:2*R);

box on;
grid on;
axis square;

file_name = 'radius_of_curvature_hypocycloid.png';


full_file_path = fullfile(file_path, file_name);

saveas(gcf, full_file_path);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%迴圈生成圖形幀並保存為GIF         %  Rotating ring gear type hypocycloid reducer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for i=1:1:ceil(360*CUT/(N+1))
%     
%     t(i) = i / (180*CUT) * pi ;
%     X = (R)*cos((N+1)*t(i))+E;
%     Y = (R)*sin((N+1)*t(i));
% 
%     c = figure('Visible', 'off');
%     hold on
%     plot([xc(i),0],[yc(i),0],'LineWidth',2);
%     plot([E,X],[0,Y],'LineWidth',2);
%     plot([X,xc(i)],[Y,yc(i)],'LineWidth',2);
% 
%     %覆蓋多於線段
% 
%     fill(1.5.*cos(0 : 0.01 : 2*pi)+X , 1.5.*sin( 0: 0.01 : 2*pi)+Y,'w-');
%     fill(1.5.*cos(0 : 0.01 : 2*pi)+xc(i) , 1.5.*sin( 0: 0.01 : 2*pi)+yc(i),'w-');
%     fill(2.5.*cos(0 : 0.01 : 2*pi) , 2.5.*sin( 0: 0.01 : 2*pi),'w-');
%     fill(2.5.*cos(0 : 0.01 : 2*pi)+E , 2.5.*sin( 0: 0.01 : 2*pi),'w-');
% 
% 
%     %非地桿接點
%     plot(1.5.*cos(0 : 0.01 : 2*pi)+X , 1.5.*sin( 0: 0.01 : 2*pi)+Y,'LineWidth',2,'Color','k');
%     plot(1.5.*cos(0 : 0.01 : 2*pi)+xc(i) , 1.5.*sin( 0: 0.01 : 2*pi)+yc(i),'LineWidth',2,'Color','k');
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
%     xlabel('X','fontname','Times New Roman','fontsize',20');
%     ylabel('Y','fontname','Times New Roman','fontsize',20');
%     title('等效連桿',SHOW);
%     
%     xlim([-(R+R/8),R+R/8]);
%     ylim([-(R+R/8),(R+R/8)]);
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
%     giffilename = fullfile(file_path, 'move_pic_hypocycloid.gif');
%     if i==1
%         imwrite(I,map,giffilename,'gif','writeMode','overwrite','LoopCount',inf,'delaytime',0.0001/CUT,'loopcount',inf);
%     else
%         imwrite(I,map,giffilename,'gif','writeMode','append','delaytime',0.0001/CUT);
%     end
% 
%     fprintf("hypocycloid_countdown %d \n",((360*CUT/(N+1))-i))
%     
%     
% 
%     close;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%生成scr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


scrfilename = fullfile(file_path, 'Cycloidal_Drive_hypocycloid.scr');
fid = fopen(scrfilename,'w');
fprintf(fid,'spline ');
fprintf(fid,'%f,%f\n', [x; y]);
fclose(fid);