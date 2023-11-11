function [Rc] = Rrcal_epicycloid(N, R, E, Rr, t)


Rc = abs(-((Rr*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) - R*sin(t) + E*N*sin(N*t))^2 + (Rr*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) - R*cos(t) + E*N*cos(N*t))^2)^(3/2)/((Rr*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) - R*cos(t) + E*N*cos(N*t))*(Rr*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((2*sin(t*(N - 1))^3*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N))^3 - (sin(t*(N - 1))*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N)) + (3*cos(t*(N - 1))*sin(t*(N - 1))*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N))^2)/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) - (((2*sin(t*(N - 1))^3*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^3 + (2*cos(t*(N - 1))*sin(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2)*((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(N*E))^2 + 1)^2) - R*cos(t) + Rr*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1)^2 + E*N^2*cos(N*t)) - (Rr*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) - R*sin(t) + E*N*sin(N*t))*(R*sin(t) + Rr*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((2*sin(t*(N - 1))^3*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N))^3 - (sin(t*(N - 1))*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N)) + (3*cos(t*(N - 1))*sin(t*(N - 1))*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N))^2)/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) - (((2*sin(t*(N - 1))^3*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^3 + (2*cos(t*(N - 1))*sin(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2)*((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(N*E))^2 + 1)^2) - Rr*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1)^2 - E*N^2*sin(N*t))));
 

% syms N R E Rr t
% 
% phi = atan(sin((1-N)*t)/((R/(E*N))-cos((1-N)*t)));
% 
% Cf = [R-Rr*cos(phi);Rr*sin(phi);1];
% 
% dx = diff(Cf(1),t);
% d2x = diff(dx,t);
% 
% dy = diff(Cf(2),t);
% d2y = diff(dy,t);
% 
% N =a;
% R = b;
% E = c;
% Rr = d;
% t = e;
% 
% Rc = ((dx)^2+(dy)^2)^1.5/(dx*d2y-dy*d2x);