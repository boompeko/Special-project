clc;
clear;

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


CUT = 1;%切割倍數
tick = 40;
SHOW = sprintf("N = %d, Rr = %d, R = %d, E = %d", N, Rr, R, E);
file_path = 'C:\Users\JOU\Desktop\git\Special-project\Cycloidal Gear\output'; %家裡電腦
%file_path = 'C:\Users\Johnny Jou\Documents\GitHub\Special-project\Cycloidal Gear\output';  %筆記電腦


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%計算曲率中心和接觸點
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : 1 : 1

    phi(i) = 0.057560592920962710041210399673131;  % phi = phi3

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
%誤差分析
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for i=1:1:1
    
    phi2(i) =  -0.057560592920962710041210399673131*(N-1)   ; % 輸入角phi2
    
 
    
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

