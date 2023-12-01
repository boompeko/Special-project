function [Rc] = newRrcal_epicycloid(N, R, E, Rr, t)


Rc = ((E*N*cos(N*t)*(N - 1) - cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(Rr - (R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) + (E*N*R*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*sin(t*(N - 1))*(N - 1))/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))^2 + (sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(Rr - (R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) - E*N*sin(N*t)*(N - 1) + (E*N*R*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*sin(t*(N - 1))*(N - 1))/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))^2)^(3/2)/((E*N*cos(N*t)*(N - 1) - cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(Rr - (R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) + (E*N*R*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*sin(t*(N - 1))*(N - 1))/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(Rr - (R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(((2*sin(t*(N - 1))^3*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N))^3 - (sin(t*(N - 1))*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N)) + (3*cos(t*(N - 1))*sin(t*(N - 1))*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N))^2)/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) - (((2*sin(t*(N - 1))^3*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^3 + (2*cos(t*(N - 1))*sin(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2)*((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(N*E))^2 + 1)^2) + cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(Rr - (R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1)^2 - E*N^2*cos(N*t)*(N - 1) + (E*N*R*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*cos(t*(N - 1))*(N - 1)^2)/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2) - (E^2*N^2*R^2*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*sin(t*(N - 1))^2*(N - 1)^2)/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(3/2) - (2*E*N*R*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*sin(t*(N - 1))*(N - 1)*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1))/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2)) + (sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(Rr - (R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1) - E*N*sin(N*t)*(N - 1) + (E*N*R*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*sin(t*(N - 1))*(N - 1))/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(Rr - (R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(((2*sin(t*(N - 1))^3*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N))^3 - (sin(t*(N - 1))*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N)) + (3*cos(t*(N - 1))*sin(t*(N - 1))*(N - 1)^2)/(cos(t*(N - 1)) - R/(E*N))^2)/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) - (((2*sin(t*(N - 1))^3*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^3 + (2*cos(t*(N - 1))*sin(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2)*((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(N*E))^2 + 1)^2) - sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*(Rr - (R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2))*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1)^2 + E*N^2*sin(N*t)*(N - 1) - (E*N*R*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*cos(t*(N - 1))*(N - 1)^2)/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2) + (E^2*N^2*R^2*sin(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*sin(t*(N - 1))^2*(N - 1)^2)/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(3/2) - (2*E*N*R*cos(t + atan(sin(t*(N - 1))/(cos(t*(N - 1)) - R/(E*N))))*sin(t*(N - 1))*(N - 1)*(((sin(t*(N - 1))^2*(N - 1))/(cos(t*(N - 1)) - R/(E*N))^2 + (cos(t*(N - 1))*(N - 1))/(cos(t*(N - 1)) - R/(E*N)))/(sin(t*(N - 1))^2/(cos(t*(N - 1)) - R/(E*N))^2 + 1) + 1))/(R^2 + E^2*N^2 - 2*E*N*R*cos(t*(N - 1)))^(1/2)));
 
