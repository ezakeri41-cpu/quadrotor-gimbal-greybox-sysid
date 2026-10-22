function [A,B,C,D,K,x0] = pitchODE(a1,a0,b,Tm,Ts,outSel,useLag)
%PITCHODE  Grey-box structure of the gimbal pitch axis.
%
%   States (useLag = false):  x = [theta ; theta_dot]
%       theta_ddot = -a1*theta_dot - a0*theta + b*u
%
%   States (useLag = true):   x = [theta ; theta_dot ; tau]
%       theta_ddot = -a1*theta_dot - a0*theta + tau
%       tau_dot    = (-tau + b*u)/Tm          (ESC + rotor spin-up lag)
%
%   u   : NORMALIZED pitch torque command (tau_y_cmd, dimensionless)
%   a1  = c/Iy     [1/s]     viscous + aero damping of the rig
%   a0  = k/Iy     [1/s^2]   restoring stiffness of the rig (CG offset from pivot)
%   b   = ku/Iy    [rad/s^2] control effectiveness, ku = physical N*m per unit command
%   Tm  [s]        actuator (ESC+rotor) time constant
%
%   outSel : 'pitch' | 'rate' | 'both'
%
%   Note: with a normalized input, only the RATIOS a1, a0, b are identifiable.
%   Iy and ku cannot be separated from this data alone (see main script).
%
%   The pitch trim offset (theta ~ -4..-5 deg at rest) is NOT modelled here.
%   It is removed by detrending the data in the main script, which is exactly
%   equivalent to working with the deviation variable theta' = theta - theta0.
%   If detrending proves insufficient, add a constant bias term to the second
%   state equation and estimate it as a fifth parameter.

if useLag
    A  = [ 0     1     0
          -a0   -a1    1
           0     0   -1/Tm ];
    B  = [0; 0; b/Tm];
else
    A  = [ 0     1
          -a0   -a1 ];
    B  = [0; b];
end

n  = size(A,1);
x0 = zeros(n,1);

switch lower(outSel)
    case 'pitch'
        C = [1 zeros(1,n-1)];
    case 'rate'
        C = [0 1 zeros(1,n-2)];
    case 'both'
        C = [eye(2) zeros(2,n-2)];
    otherwise
        error('outSel must be ''pitch'', ''rate'' or ''both''.');
end

D = zeros(size(C,1),1);
K = zeros(n,size(C,1));     % output-error structure (no noise model)
end
