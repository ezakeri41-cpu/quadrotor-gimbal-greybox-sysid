function [A,B,C,D,K,x0] = mimoFull_ODE(a1r,a0r,a1p,a0p,a1y,a0y, ...
                                       bxx,bxy,bxz,byx,byy,byz,bzx,bzy,bzz,Ts) %#ok<INUSD>
%MIMOFULL_ODE  3-input / 6-output grey-box, full 3x3 torque B matrix, no actuator lag.
%
%  x = [phi p  theta q  psi r],   u = [tau_x tau_y tau_z],   y = [angles ; rates]
%
%  phi_dd   = -a0r*phi   - a1r*p + bxx*u_x + bxy*u_y + bxz*u_z
%  theta_dd = -a0p*theta - a1p*q + byx*u_x + byy*u_y + byz*u_z
%  psi_dd   = -a0y*psi   - a1y*r + bzx*u_x + bzy*u_y + bzz*u_z
%
%  Naming: b<axis><input>, e.g. bxy = effect of tau_y on roll.

A = zeros(6);
A(1,2) = 1;  A(2,1) = -a0r;  A(2,2) = -a1r;
A(3,4) = 1;  A(4,3) = -a0p;  A(4,4) = -a1p;
A(5,6) = 1;  A(6,5) = -a0y;  A(6,6) = -a1y;

B = zeros(6,3);
B(2,:) = [bxx bxy bxz];
B(4,:) = [byx byy byz];
B(6,:) = [bzx bzy bzz];

C = zeros(6,6);
C(1,1) = 1;  C(2,3) = 1;  C(3,5) = 1;     % roll pitch yaw
C(4,2) = 1;  C(5,4) = 1;  C(6,6) = 1;     % rollrate pitchrate yawrate

D  = zeros(6,3);
K  = zeros(6,6);
x0 = zeros(6,1);
end
