function [A,B,C,D,K,x0] = rollODE(a1,a0,b,Tm,Ts,outSel,useLag)
%ROLLODE  Grey-box structure of the gimbal roll axis.
%
%   States (useLag = false):  x = [phi ; phi_dot]
%       phi_ddot = -a1*phi_dot - a0*phi + b*u
%
%   States (useLag = true):   x = [phi ; phi_dot ; tau]
%       phi_ddot = -a1*phi_dot - a0*phi + tau
%       tau_dot  = (-tau + b*u)/Tm          (ESC + rotor spin-up lag)
%
%   u   : NORMALIZED roll torque command (tau_x_cmd, dimensionless)
%   a1  = c/Ix     [1/s]     viscous + aero damping of the rig
%   a0  = k/Ix     [1/s^2]   restoring stiffness of the rig (CG offset from pivot)
%   b   = ku/Ix    [rad/s^2] control effectiveness, ku = physical N*m per unit command
%   Tm  [s]        actuator (ESC+rotor) time constant
%
%   outSel : 'roll' | 'rate' | 'both'
%
%   Note: with a normalized input, only the RATIOS a1, a0, b are identifiable.
%   Ix and ku cannot be separated from this data alone (see main script).

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
    case 'roll'
        C = [1 zeros(1,n-1)];
    case 'rate'
        C = [0 1 zeros(1,n-2)];
    case 'both'
        C = [eye(2) zeros(2,n-2)];
    otherwise
        error('outSel must be ''roll'', ''rate'' or ''both''.');
end

D = zeros(size(C,1),1);
K = zeros(n,size(C,1));     % output-error structure (no noise model)
end
