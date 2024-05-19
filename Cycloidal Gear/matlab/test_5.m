% 定义要解决的方程
syms x y;
eqns = [x^2 + y^2 == 1, x == y^2];

% 使用 vpasolve 求解方程
sol = vpasolve(eqns, [x, y]);
disp('vpasolve 解的结果:');
disp(sol);
