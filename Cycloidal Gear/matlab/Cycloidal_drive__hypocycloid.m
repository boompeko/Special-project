
clear;
clc;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%參數
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       !!!!!!!!!!!!!!!!!!!  (R/N) > E   !!!!!!!!!!!!!!!!!!!!!!

% N = 15 ;%Number of rollers
% Rr = 18 ;%Radius of the roller
% R = 120 ;%Radius of the rollers PCD (Pitch Circle Diamater)
% E = 5 ;% Eccentricity - offset from input shaft to a cycloidal disk

N = 9 ;%Number of rollers
Rr = 15 ;%Radius of the roller
R = 90 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E = 9 ;% Eccentricity - offset from input shaft to a cycloidal disk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作圖
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CUT = 5;%切割倍數
tick = 40;
SHOW = sprintf("N = %d, Rr = %d, R = %d, E = %d", N, Rr, R, E);
file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\output'; %家裡電腦
%file_path = 'C:\Users\Johnny Jou\Documents\GitHub\Special-project\Cycloidal Gear\output';  %筆記電腦

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
hold on
plot(x,y,'LineWidth',2,'Color','b');
xlabel('X (mm)','fontname','Times New Roman','fontsize',20');
ylabel('Y (mm)','fontname','Times New Roman','fontsize',20');
title(SHOW,'內擺線輪輪廓','fontsize',16);
axis equal;
xlim([-(1.3*R),(1.3*R)]);
ylim([-(1.3*R),(1.3*R)]);
xticks(-2*R:tick:2*R);
yticks(-2*R:tick:2*R);

box on;
grid on;
axis square;
hold off
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
xlabel('cycloidal disk rotation angle (θ)','fontname','Times New Roman','fontsize',20');
ylabel('radius of curvature (mm)','fontname','Times New Roman','fontsize',20');
title(SHOW,'曲率半徑','fontsize',16);
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
%     xlabel('X (mm)','fontname','Times New Roman','fontsize',20');
%     ylabel('Y (mm)','fontname','Times New Roman','fontsize',20');
%     title(SHOW,'等效連桿');
%     
%     xlim([-(R+R/5+2*E),R+R/5+2*E]);
%     ylim([-(R+R/8+2*E),(R+R/8)+2*E]);
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
%     fprintf("hypocycloid_countdown %d \n",(ceil((360*CUT/(N+1))-i)))
%     
%     
% 
%     close;
% end
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%生成scr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


scrfilename = fullfile(file_path, 'Cycloidal_Drive_hypocycloid.scr');
fid = fopen(scrfilename,'w');
fprintf(fid,'spline ');
fprintf(fid,'%f,%f\n', [x; y]);
fclose(fid);



for i=1:1:ceil(360*CUT/(N+1))
    
    t(i) = i / (180*CUT) * pi ;
    
    x = (E)*cos(t(i));
    y = (E)*sin(t(i));
    
    zeta1(i) = -t(i)/N ;
    zeta(i) = -atan((Ory - y)/(Orx - x)) ;

    A(i) = zeta1(i)/zeta(i);
    
    OM(i) = atan(Ky(i)-Kx(i));
    
    a1(i) = (Kx(i))^2 + (Ky(i))^2;
    a2(i) = newRc(i)+Rr;
    a3(i) = R;
    a4(i) = E;

    error = 0.015;
    
    H(i) = -2*a4(i)*a3(i)*cos(t(i));
    
    I(i) = -2*a4(i)*a3(i)*cos(t(i)) + 2*a1(i)*a3(i);

    J(i) = a1(i)^2 - a2(i)^2 + a3(i)^2 + a4(i)^2 - 2*a1(i)*a4(i)*cos(t(i));

    

    E1(i) = ((-a3(i)*cos(zeta(i))+a1(i)-a4(i)*cos(t(i)))*2*error)/(H(i)*cos(zeta(i))-I(i)*sin(zeta(i)));

    E2(i) = ((a2(i))*2*error)/(H(i)*cos(zeta(i))-I(i)*sin(zeta(i)));
    
    E3(i) = ((a4(i)*cos(t(i)-zeta(i))+a3(i)-a1(i)*cos(t(i)))*2*error)/(H(i)*cos(zeta(i))-I(i)*sin(zeta(i)));

    E4(i) = ((a3(i)*cos(t(i)-zeta(i))+a4(i)-a1(i)*cos(t(i)))*2*error)/(H(i)*cos(zeta(i))-I(i)*sin(zeta(i)));
    

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


j = figure('Visible', 'on');
hold on
plot(f,A,'LineWidth',2,'Color','b');
