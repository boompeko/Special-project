
clear;
clc;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%參數
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       !!!!!!!!!!!!!!!!!!!  (R/(N*E)) > 1   !!!!!!!!!!!!!!!!!!!!!!

N = 7 ;%Number of rollers
Rr = 10 ;%Radius of the roller
R = 80 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E =10 ;% Eccentricity - offset from input shaft to a cycloidal disk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作圖
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CUT = 1;%切割倍數
tick = 10;
SHOW = sprintf("N = %d, Rr = %d, R = %d, E = %d", N, Rr, R, E);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%計算曲率中心和接觸點
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1 : 1 : ((360*CUT)+1)

    t(i) = i / (180*CUT) * pi;
    phi(i) = atan(sin((1-N)*t(i))/((R/(E*N))-cos((1-N)*t(i))));
    
    Cf = [R-Rr*cos(phi(i));Rr*sin(phi(i));1];
    
    [Rc(i)] = Rrcal_epicycloid(N,R,E,Rr,t(i));

  
    Ccf = [R-(Rr+Rc(i))*cos(phi(i));(Rr+Rc(i))*sin(phi(i));1];

    M1 = [cos(-N*t(i)), -sin(-N*t(i)), -E*cos(-N*t(i));
          sin(-N*t(i)), cos(-N*t(i)), -E*sin(-N*t(i));
          0, 0, 1];
    M2 = [cos((1-N)*t(i)), sin((1-N)*t(i)), 0;
          -sin((1-N)*t(i)), cos((1-N)*t(i)), 0;
          0, 0, 1];
    
    C = M1*M2*Cf;
    Cc = M1*M2*Ccf;
%     
    x(i) = C(1,:);
    y(i) = C(2,:);

    xc(i) = Ccf(1,:);
    yc(i) = Ccf(2,:);
%     

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
xlim([-(R+R/10),(R+R/10)]);
ylim([-(R+R/10),(R+R/10)]);
xticks(-2*R:tick:2*R);
yticks(-2*R:tick:2*R);

box on;
grid on;
axis square;

file_name = 'contour.png';
file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\movie\pic'; 

full_file_path = fullfile(file_path, file_name);

saveas(gcf, file_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%繪製曲率半徑
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z(1) = 0;
bRc(1) = Rc(360*CUT/(N-1));

for i = 1 : 1 : (360*CUT/(N-1))

    z(i+1) = i/CUT ;
    bRc(i+1) = Rc(i);

end


b = figure('Visible', 'on');
plot(z,bRc,'LineWidth',2);
xlabel('cycloidal disk rotation angle(θ)','fontname','Times New Roman','fontsize',20');
ylabel('radius of curvature','fontname','Times New Roman','fontsize',20');
title('曲率半徑',SHOW);
xlim([0,(360/(N-1))]);
ylim([0,R]);
xticks(0:5:(360/(N-1)));
yticks(0:tick:R);

box on;
grid on;
axis square;

file_name = 'radius_of_curvature.png';
file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\movie\pic'; 

full_file_path = fullfile(file_path, file_name);

saveas(gcf, file_name);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%循环生成图形帧并保存为GIF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:1:(360*CUT/(N-1))
    
    t(i) = i / (180*CUT) * pi ;
    X = (E)*cos((1-N)*t(i));
    Y = (E)*sin((1-N)*t(i));

    c = figure('Visible', 'off');
    hold on
    plot([xc(i),R],[yc(i),0],'LineWidth',2);
    plot([0,X],[0,Y],'LineWidth',2);
    plot([X,xc(i)],[Y,yc(i)],'LineWidth',2);

    %覆蓋多於線段

    fill(1.5.*cos(0 : 0.01 : 2*pi)+X , 1.5.*sin( 0: 0.01 : 2*pi)+Y,'w-');
    fill(1.5.*cos(0 : 0.01 : 2*pi)+xc(i) , 1.5.*sin( 0: 0.01 : 2*pi)+yc(i),'w-');
    fill(2.5.*cos(0 : 0.01 : 2*pi) , 2.5.*sin( 0: 0.01 : 2*pi),'w-');
    fill(2.5.*cos(0 : 0.01 : 2*pi)+R , 2.5.*sin( 0: 0.01 : 2*pi),'w-');


    %非地桿接點
    plot(1.5.*cos(0 : 0.01 : 2*pi)+X , 1.5.*sin( 0: 0.01 : 2*pi)+Y,'LineWidth',2,'Color','k');
    plot(1.5.*cos(0 : 0.01 : 2*pi)+xc(i) , 1.5.*sin( 0: 0.01 : 2*pi)+yc(i),'LineWidth',2,'Color','k');
    


    
    %原點地桿

    offsetX = 0;
    offsetY = 0;
    
    
    plot([4+offsetX,4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
    plot([-4+offsetX,-4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
    plot([-8+offsetX,8+offsetX],[-8+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot([-8+offsetX,-4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot([-4+offsetX,0+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot([0+offsetX,4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot([4+offsetX,8+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot(4.*cos(0 : 0.01 : 2*pi)+offsetX , 4.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
    plot(2.5.*cos(0 : 0.01 : 2*pi)+offsetX , 2.5.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');


    %針輪地桿

    offsetX = R;
    offsetY = 0;
    

    plot([4+offsetX,4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
    plot([-4+offsetX,-4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
    plot([-8+offsetX,8+offsetX],[-8+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot([-8+offsetX,-4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot([-4+offsetX,0+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot([0+offsetX,4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot([4+offsetX,8+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
    plot(4.*cos(0 : 0.01 : 2*pi)+offsetX , 4.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
    plot(2.5.*cos(0 : 0.01 : 2*pi)+offsetX , 2.5.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
    
    box on;
    grid on;
    axis equal;
    xlabel('X','fontname','Times New Roman','fontsize',20');
    ylabel('Y','fontname','Times New Roman','fontsize',20');
    title('等效連桿',SHOW);
    
    xlim([-2*E,R+R/8]);
    ylim([-(R+R/8)/2,(R+R/8)/2]);
    xticks(-2*R:tick:2*R);
    yticks(-2*R:tick:2*R);
    
    hold off
    
    % 捕獲圖形幀
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    output_folder = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\movie';
    giffilename = fullfile(output_folder, 'move_pic.gif');
    if i==1
        imwrite(I,map,'move_pic.gif','gif','writeMode','overwrite','LoopCount',inf,'delaytime',0.0001,'loopcount',inf);
    else
        imwrite(I,map,'move_pic.gif','gif','writeMode','append','delaytime',0.0001);
    end

    fprintf("countdown %d \n",((360*CUT/(N-1))-i))
    
    

    close;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%生成scr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_folder = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\CAD';
scrfilename = fullfile(output_folder, 'Cycloidal_Drive_ epicycloid.scr');
fid = fopen(scrfilename,'w');
fprintf(fid,'spline ');
fprintf(fid,'%f,%f\n', [x; y]);
fclose(fid);