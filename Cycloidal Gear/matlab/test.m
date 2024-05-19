clc;
clear;

syms N x R E Rr t



% X = (R*cos(t))-(Rr*cos(t+atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t)))))-(E*cos(N*t));
% Y = (-R*sin(t))+(Rr*sin(t+atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t)))))+(E*sin(N*t));
% 
% X = (R*cos(t))+Rr*cos(t+atan(sin((1+N)*t)/((R/(E*N))-cos((1+N)*t))))+(E*cos(N*t));
% Y = (R*sin(t))+Rr*sin(t+atan(sin((1+N)*t)/((R/(E*N))-cos((1+N)*t))))-(E*sin(N*t));

%newway
syms N x R E Rr t


% X = cos(N*t)*E*(N-1) + cos(t+atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t))))*((R^2+(E*N)^2-2*R*(E*N)*cos((1-N)*t))^0.5 - Rr);
% Y = -sin(N*t)*E*(N-1) - sin(t+atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t))))*((R^2+(E*N)^2-2*R*(E*N)*cos((1-N)*t))^0.5 - Rr);


X = cos(N*t)*E*(N+1) + cos(t-atan(-sin((1+N)*t)/((R/(E*N))-cos((1+N)*t))))*((R^2+(E*N)^2-2*R*(E*N)*cos((1+N)*t))^0.5 + Rr);
Y = -sin(N*t)*E*(N+1) + sin(t-atan(-sin((1+N)*t)/((R/(E*N))-cos((1+N)*t))))*((R^2+(E*N)^2-2*R*(E*N)*cos((1+N)*t))^0.5 + Rr);


dx = diff(X,t);
d2x = diff(dx,t);

dy = diff(Y,t);
d2y = diff(dy,t);


Rc = ((dx)^2+(dy)^2)^1.5/(dx*d2y-dy*d2x);

k = 1/Rc == 0;

ans = vpasolve(k,t);

ans
% disp(dx);
% disp(dy);
% disp(d2x);
% disp(d2y);


%final_ = simplify(Rc);


%disp(final_);