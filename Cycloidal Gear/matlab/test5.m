clc;
clear;

syms e1 e_r e_Rr e4;



phi2 =  -1 / (180) * pi   ; % 輸入角phi2




a1 = 120;
a2 = -2.019998129220696;
a3 = 1.170199971467818e+02;
a4 = 5;



H = -2*a3*a4*sin(phi2);

I = -2*a3*a4*cos(phi2) + 2*a1*a3;

J = a1^2 - a2^2 + a3^2 + a4^2 - 2*a1*a4*cos(phi2);


theta = 0.001247468139883; % 輸出角 依照點與點之間位置算

error = 0.022; %長度誤差

error_1 = 0.022;



E1_eqn = (H)*sin(theta+e1)+(I+2*error_1*a3)*cos(theta+e1)-(J+2*error_1*a1-2*error_1*a4*cos(phi2));

T1 = vpasolve(E1_eqn == 0, e1);
E1_sol = T1 *180/pi;

error_r = 0.022;
    


E2_r_eqn = (H)*sin(theta+e_r)+(I)*cos(theta+e_r)-(J-2*error_r*a2);

Tr = vpasolve(E2_r_eqn == 0, e_r);

if isempty(Tr)
    E2_r_sol = NaN;
else
    E2_r_sol =  Tr*180*pi;

end

error_Rr = 0.011;

E2_Rr_eqn = (H)*sin(theta+e_Rr)+(I)*cos(theta+e_Rr)-(J-2*error_Rr*a2);


TRr = vpasolve(E2_Rr_eqn == 0, e_Rr);
if isempty(TRr)
    E2_Rr_sol = NaN;
else
    E2_Rr_sol =  TRr*180*pi;

end

error_4 = 0.008; 
E4_eqn = (H-2*error_4*a3*sin(phi2))*sin(theta+e4)+(I-2*error_4*a3*cos(phi2)*cos(theta+e1))-(J+2*error_4*a4-2*error_4*a1*cos(phi2));

T4 = vpasolve(E4_eqn == 0, e4);
E4_sol = T4*180/pi;

% sol = fsolve(eqns == 0,0);
% disp('vpasolve 解的结果:');
% disp(sol*180/3.1415926);


