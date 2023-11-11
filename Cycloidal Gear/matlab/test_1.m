clear;
clc;
close all;

offsetX = 0;
offsetY = 0;

figure;
hold on 
axis ([-80 60 -60 80]);

% 绘制完整的线段
x = 0 :30;
y = x;
plot(x, y, 'b-', 'LineWidth', 2);


% 擦除部分线段，例如从x=4到x=7之间
x_erase = 20:25;
% y_erase = x_erase;
y_erase = NaN(size(x_erase));
plot(x_erase, y_erase, 'w-', 'LineWidth', 2);  % 使用白色线段擦除

plot([4+offsetX,4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
plot([-4+offsetX,-4+offsetX],[-8+offsetY,0+offsetY],'LineWidth',2,'Color','k')
plot([-8+offsetX,8+offsetX],[-8+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
plot([-8+offsetX,-4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
plot([-4+offsetX,0+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
plot([0+offsetX,4+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
plot([4+offsetX,8+offsetX],[-11+offsetY,-8+offsetY],'LineWidth',2,'Color','k')
plot(4.*cos(0 : 0.01 : 2*pi)+offsetX , 4.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
plot(2.5.*cos(0 : 0.01 : 2*pi)+offsetX , 2.5.*sin( 0: 0.01 : 2*pi)+offsetY,'LineWidth',2,'Color','k');
fill(2.5.*cos(0 : 0.01 : 2*pi)+offsetX , 2.5.*sin( 0: 0.01 : 2*pi)+offsetY,'w-')



hold off;

% 添加标题和标签
title('擦除部分线段的效果');
xlabel('X轴');
ylabel('Y轴');


box on;
grid on;
axis equal;
