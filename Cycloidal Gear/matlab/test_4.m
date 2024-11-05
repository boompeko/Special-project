clc;
clear;


syms N R E Rr x


% % X = cos(N*t)*E*(N-1) + cos(t+atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t))))*((R^2+(E*N)^2-2*R*(E*N)*cos((1-N)*t))^0.5 - Rr);
% % Y = -sin(N*t)*E*(N-1) - sin(t+atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t))))*((R^2+(E*N)^2-2*R*(E*N)*cos((1-N)*t))^0.5 - Rr);
% 
% 
% X = cos(15*x)*5*(16) + cos(x-atan(-sin((16)*x)/((120/(5*15))-cos((16)*x))))*((120^2+(5*15)^2-2*120*(5*15)*cos((16)*x))^0.5 + 18);
% Y = -sin(15*x)*5*(16) + sin(x-atan(-sin((16)*x)/((120/(5*15))-cos((16)*x))))*((120^2+(5*15)^2-2*120*(5*15)*cos((16)*x))^0.5 + 18);
% 
% 
% dx = diff(X,x);
% d2x = diff(dx,x);
% 
% dy = diff(Y,x);
% d2y = diff(dy,x);
% 
% 
% Rc = ((dx)^2+(dy)^2)^1.5/(dx*d2y-dy*d2x);

Rc = ((5*5*cos(5*x)*(5 - 1) - cos(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*(18 - (120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5)))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) + 1) + (5*5*120*sin(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*sin(x*(5 - 1))*(5 - 1))/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))^2 + (sin(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*(18 - (120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5)))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) + 1) - 5*5*sin(5*x)*(5 - 1) + (5*5*120*cos(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*sin(x*(5 - 1))*(5 - 1))/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))^2)^(3/2)/((5*5*cos(5*x)*(5 - 1) - cos(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*(18 - (120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5)))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) + 1) + (5*5*120*sin(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*sin(x*(5 - 1))*(5 - 1))/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(sin(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*(18 - (120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(((2*sin(x*(5 - 1))^3*(5 - 1)^2)/(cos(x*(5 - 1)) - 120/(5*5))^3 - (sin(x*(5 - 1))*(5 - 1)^2)/(cos(x*(5 - 1)) - 120/(5*5)) + (3*cos(x*(5 - 1))*sin(x*(5 - 1))*(5 - 1)^2)/(cos(x*(5 - 1)) - 120/(5*5))^2)/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) - (((2*sin(x*(5 - 1))^3*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^3 + (2*cos(x*(5 - 1))*sin(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2)*((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1)^2) + cos(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*(18 - (120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5)))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) + 1)^2 - 5*5^2*cos(5*x)*(5 - 1) + (5*5*120*cos(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*cos(x*(5 - 1))*(5 - 1)^2)/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2) - (5^2*5^2*120^2*cos(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*sin(x*(5 - 1))^2*(5 - 1)^2)/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(3/2) - (2*5*5*120*sin(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*sin(x*(5 - 1))*(5 - 1)*(((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5)))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) + 1))/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2)) + (sin(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*(18 - (120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5)))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) + 1) - 5*5*sin(5*x)*(5 - 1) + (5*5*120*cos(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*sin(x*(5 - 1))*(5 - 1))/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(cos(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*(18 - (120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(((2*sin(x*(5 - 1))^3*(5 - 1)^2)/(cos(x*(5 - 1)) - 120/(5*5))^3 - (sin(x*(5 - 1))*(5 - 1)^2)/(cos(x*(5 - 1)) - 120/(5*5)) + (3*cos(x*(5 - 1))*sin(x*(5 - 1))*(5 - 1)^2)/(cos(x*(5 - 1)) - 120/(5*5))^2)/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) - (((2*sin(x*(5 - 1))^3*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^3 + (2*cos(x*(5 - 1))*sin(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2)*((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1)^2) - sin(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*(18 - (120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2))*(((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5)))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) + 1)^2 + 5*5^2*sin(5*x)*(5 - 1) - (5*5*120*sin(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*cos(x*(5 - 1))*(5 - 1)^2)/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2) + (5^2*5^2*120^2*sin(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*sin(x*(5 - 1))^2*(5 - 1)^2)/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(3/2) - (2*5*5*120*cos(x + atan(sin(x*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))))*sin(x*(5 - 1))*(5 - 1)*(((sin(x*(5 - 1))^2*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5))^2 + (cos(x*(5 - 1))*(5 - 1))/(cos(x*(5 - 1)) - 120/(5*5)))/(sin(x*(5 - 1))^2/(cos(x*(5 - 1)) - 120/(5*5))^2 + 1) + 1))/(120^2 + 5^2*5^2 - 2*5*5*120*cos(x*(5 - 1)))^(1/2)));
 


k = 1/Rc == 0;



ans = vpasolve(k,x,0.05);

ans

a = 0.057560592920962710041210399673131*180*14/pi;



disp(a);


% disp(dx);
% disp(dy);
% disp(d2x);
% disp(d2y);


%final_ = simplify(Rc);


%disp(final_);