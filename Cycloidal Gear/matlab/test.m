clc;
clear;

syms N x R E Rr t






% X = (R*cos(t))-(Rr*cos(t+atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t)))))-(E*cos(N*t));
% Y = (-R*sin(t))+(Rr*sin(t+atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t)))))+(E*sin(N*t));

X = (R*cos(t))+Rr*cos(t+atan(sin((1+N)*t)/((R/(E*N))-cos((1+N)*t))))+(E*cos(N*t));
Y = (R*sin(t))+Rr*sin(t+atan(sin((1+N)*t)/((R/(E*N))-cos((1+N)*t))))-(E*sin(N*t));



% phi = atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t)));
% 
% Cf = [R-Rr*cos(phi);Rr*sin(phi);1];
% 
% X = R-Rr*cos(atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t))));
% Y = Rr*sin(atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t))));

dx = diff(X,t);
d2x = diff(dx,t);

dy = diff(Y,t);
d2y = diff(dy,t);


Rc = ((dx)^2+(dy)^2)^1.5/(dx*d2y-dy*d2x);

% rt = ((dx)^2+(dy)^2)^1.5;
% 
% L = limit(rt, t, Inf);
% 
% disp(L);

% disp(((dx)^2+(dy)^2)^1.5);
% 
% disp(dx);
% disp(dy);
% disp(d2x);
% disp(d2y);
disp(Rc);
% 
% a = ((Rr*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) - R*sin(t) + E*N*sin(N*t))^2 + (Rr*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) - R*cos(t) + E*N*cos(N*t))^2)^(3/2);
%  