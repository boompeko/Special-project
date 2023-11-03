clc;
clear;



syms x
f = 1/x;
L = limit(f, x, Inf);

disp(L);
