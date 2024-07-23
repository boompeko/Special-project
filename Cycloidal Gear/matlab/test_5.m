clc;
clear;

a=[1,2,3,4,5];
b=[1,2,3,4,5];

figure;
plot(a,b,LineStyle=":");
figure;
plot(a,b,LineStyle="-");
figure;
plot(a,b,LineStyle="none");
figure;
plot(a,b,LineStyle="-.");
figure;
plot(a,b,LineStyle="--");