clc;
clear;

syms a b c d  theta zeta phi

A = (d*cos(phi)+c*cos(theta)-a*cos(zeta))^2 + (d*sin(phi)+c*sin(theta)-a*sin(zeta))^2;

B = (d*cos(phi-zeta)+c*cos(theta-zeta)-a)^2 + (d*sin(phi-zeta)+c*sin(theta-zeta))^2;

C = A - B;

disp(simplify(collect(expand(A))));
disp(simplify(collect(expand(B))));
disp(simplify(collect(expand(C))));

 
 