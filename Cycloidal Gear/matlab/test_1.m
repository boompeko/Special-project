clear;
clc;
close all;

syms theta zeta a1 a2 a3 ;


a4 = (( a1*cos(theta)+ a3*cos(zeta) - a2)^2+(a1*sin(theta) + a3*sin(zeta))^2);
a4 = expand(a4);
a4 = collect(a4);
% a4 = simplify(a4);

a4

