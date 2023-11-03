
clear;
clc;



N = 13 ;%Number of rollers
Rr = 5 ;%Radius of the roller
R = 40 ;%Radius of the rollers PCD (Pitch Circle Diamater)
E =2 ;% Eccentricity - offset from input shaft to a cycloidal disk



for i = 1 : 1 : 361

    t(i) = i / 180 * pi;

    phi(i) = atan(sin((1+N)*t(i))/((R/(E*N))-cos((1+N)*t(i))));
    
    Cf = [R+Rr*cos(phi(i));Rr*sin(phi(i));1];
%     
%     [Rc(i)] = Rrcal(N,R,E,Rr,t(i));
% 
%   
%     Ccf = [R-(Rr+Rc(i))*cos(phi(i));(Rr+Rc(i))*sin(phi(i));1];

    M1 = [cos(t(i)), sin(t(i)), 0;
          -sin(t(i)), cos(t(i)), 0;
          0, 0, 1];
    M2 = [cos((1+N)*t(i)), -sin((1+N)*t(i)), E;
          sin((1+N)*t(i)), cos((1+N)*t(i)), 0;
          0, 0, 1];
    
    C = M1*M2*Cf;
%     Cc = M1*M2*Ccf;
    
    x(i) = C(1,:);
    y(i) = C(2,:);

%     xc(i) = Ccf(1,:);
%     yc(i) = Ccf(2,:);
    

end



figure;
plot(x,y,'LineWidth',2);
xlabel('X','fontname','Times New Roman','fontsize',20');
ylabel('Y','fontname','Times New Roman','fontsize',20');
axis equal;
xlim([-50,50]);
ylim([-50,50]);
xticks(-50:10:50);
yticks(-50:10:50);

% for i = 1 : 1 : 30000
% 
%     t(i) = i /1000;
% end
% 
% figure;
% plot(t,Rc,'LineWidth',2);
% xlabel('X','fontname','Times New Roman','fontsize',20');
% ylabel('Y','fontname','Times New Roman','fontsize',20');
% axis equal;
% xlim([0,30]);
% ylim([0,10000]);
% xticks(0:1:30);
% yticks(0:100:10000);

% for i=1:1:361
%     
%     t(i) = i / 180 * pi / 12;
%     X = E*cos(t(i));
%     Y = E*sin(t(i));
% 
%     figure;
%     hold on
%     plot([xc(i),R],[yc(i),0],'LineWidth',2);
%     plot([0,X],[0,Y],'LineWidth',2);
%     plot([X,xc(i)],[Y,yc(i)],'LineWidth',2);
%     
%     xlabel('X','fontname','Times New Roman','fontsize',20');
%     ylabel('Y','fontname','Times New Roman','fontsize',20');
%     axis equal;
%     xlim([-50,50]);
%     ylim([-50,50]);
%     xticks(-50:10:50);
%     yticks(-50:10:50);
%     box on;
%     grid on;
% 
%     save_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\gif';  % 更改成你想要的路徑
%     file_name =  sprintf('plot_%d.png', i);
%     
%     % 儲存圖像到指定路徑
%     full_file_path = fullfile(save_path, file_name);
%     saveas(gcf, full_file_path);
% 
%     close all;
% 
% end








% fid = fopen('Cycloidal_Drive.scr','w');
% fprintf(fid,'spline ');
% fprintf(fid,'%f,%f\n', [x; y]);
% fclose(fid);