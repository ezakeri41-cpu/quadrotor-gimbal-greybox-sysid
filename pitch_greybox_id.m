%% PITCH_GREYBOX_ID  Grey-box system identification of a 3-DOF gimbal pitch axis
%  Estimation : standalone_chirp_pitch_log(0.1amp0.2nbt).csv
%  Validation : standalone_chirp_pitch_log(0.1amp0.2nbt2).csv
%  Requires   : System Identification Toolbox + pitchODE.m on the path
%
%  Log columns: time, tau_y_cmd, roll_rad, pitch_rad, yaw_rad,
%               pitchspeed_rad_s, pwm1..pwm4
%
%  MODEL (see pitchODE.m):
%      theta_ddot = -a1*theta_dot - a0*theta + b*u   [+ first-order actuator lag]
%  As on the roll axis, a0 (stiffness) and a1 (damping) belong to the TEST RIG
%  (CG offset from the pivot, bearing friction, cable stiffness), not to the
%  free-flying quadrotor. They are identified anyway so they do not corrupt b.
%
%  PITCH-SPECIFIC NOTE: this axis sits at a trim of about -4..-5 deg at rest.
%  Detrending each record with its own mean is exactly equivalent to using the
%  deviation variable theta' = theta - theta0, so the same pipeline applies.
%  If the trim drifts within a run rather than staying constant, switch
%  DETREND_ORDER to 1 (linear detrend).

clear; clc; close all;

%% ---------------- User settings ----------------
estFile   = 'standalone_chirp_pitch_log(0.1amp0.2nbt).csv';
valFile   = 'standalone_chirp_pitch_log(0.1amp0.2nbt2).csv';

Ts        = 0.02;      % [s] uniform resampling period (logs are ~49.7 Hz)
USE_LAG   = true;      % true -> 3rd-order model with ESC/rotor lag Tm
USE_DELAY = true;      % true -> run delayest and put the delay in the model
DETREND_ORDER = 0;     % 0 = remove mean (constant trim), 1 = remove linear drift
LPFILTER  = false;     % optional zero-phase low-pass of the outputs
fcut      = 8;         % [Hz] cut-off if LPFILTER = true

% Chirp excitation band: roughly 0.3 ... 2 Hz.
% The identified model is only trustworthy inside this band.
fBand     = [0.3 2];   % [Hz]

%% ---------------- Load + resample onto a uniform grid ----------------
% iddata assumes a constant sample time; the logs have small timing jitter
% (dt = 20.11 +/- 0.17 ms), so interpolate onto an exact grid first.
[uE,yE,ypE,tE] = loadLog(estFile,Ts);
[uV,yV,ypV,tV] = loadLog(valFile,Ts);

% Quick look at the trim: is it constant, or drifting within the run?
fprintf('Estimation run: mean pitch %.2f deg, last 5 s mean %.2f deg\n', ...
        rad2deg(mean(yE)), rad2deg(mean(yE(tE > tE(end)-5))));
fprintf('Validation run: mean pitch %.2f deg, last 5 s mean %.2f deg\n', ...
        rad2deg(mean(yV)), rad2deg(mean(yV(tV > tV(end)-5))));
fprintf('If those two numbers differ by more than ~1 deg, set DETREND_ORDER = 1.\n\n');

if LPFILTER
    [bf,af] = butter(4,fcut/(0.5/Ts));
    yE  = filtfilt(bf,af,yE);   ypE = filtfilt(bf,af,ypE);
    yV  = filtfilt(bf,af,yV);   ypV = filtfilt(bf,af,ypV);
end

%% ---------------- iddata objects + detrending ----------------
% Offsets are physical and run-specific (pitch trim, IMU mounting bias), so
% each record is detrended with its OWN mean.
zE_pitch = iddata(yE ,uE,Ts,'OutputName','pitch','OutputUnit','rad', ...
                            'InputName','tau_y_cmd','InputUnit','-');
zE_rate  = iddata(ypE,uE,Ts,'OutputName','pitchrate','OutputUnit','rad/s', ...
                            'InputName','tau_y_cmd');
zE_both  = iddata([yE ypE],uE,Ts,'OutputName',{'pitch','pitchrate'});

zV_pitch = iddata(yV ,uV,Ts,'OutputName','pitch','InputName','tau_y_cmd');
zV_rate  = iddata(ypV,uV,Ts,'OutputName','pitchrate','InputName','tau_y_cmd');
zV_both  = iddata([yV ypV],uV,Ts,'OutputName',{'pitch','pitchrate'});

zE_pitch = detrend(zE_pitch,DETREND_ORDER);  zV_pitch = detrend(zV_pitch,DETREND_ORDER);
zE_rate  = detrend(zE_rate ,DETREND_ORDER);  zV_rate  = detrend(zV_rate ,DETREND_ORDER);
zE_both  = detrend(zE_both ,DETREND_ORDER);  zV_both  = detrend(zV_both ,DETREND_ORDER);

%% ---------------- Input/output delay estimate ----------------
nk_pitch = delayest(zE_pitch,3,1,0,25);   % search 0..25 samples (0..0.5 s)
nk_rate  = delayest(zE_rate ,3,1,0,25);
fprintf('delayest: pitch %d samples (%.0f ms), rate %d samples (%.0f ms)\n', ...
        nk_pitch,1e3*nk_pitch*Ts,nk_rate,1e3*nk_rate*Ts);
fprintf('(These are usually resonance phase, not true dead time - see roll axis.)\n');

% With USE_LAG = true the ESC/rotor lag is already a model state, so a separate
% transport delay is unnecessary and the two compete for the same phase lag.
if USE_DELAY && ~USE_LAG
    delayPitch = nk_pitch*Ts;  delayRate = nk_rate*Ts;
else
    delayPitch = 0;            delayRate = 0;
end

%% ---------------- Grey-box model template ----------------
% Initial guesses from a quick equation-error / simulation-error fit of this
% exact data set: wn ~ 5.5 rad/s (a0 ~ 30), a1 ~ 2.2, b ~ 35, Tm ~ 0.15 s.
% idgrey wants an N-by-2 cell array {name, value}; bounds and Free flags
% are set afterwards on m.Structure.Parameters.
par0 = {'a1', 2.2 ; 'a0', 30 ; 'b', 35 ; 'Tm', 0.15};

parMin  = [0   0   -Inf 0.01];
parMax  = [Inf Inf  Inf 0.30];   % Tm > ~0.3 s is not a physical ESC/rotor lag
parFree = [true true true USE_LAG];   % Tm frozen if the lag state is unused

opt = greyestOptions('InitialState','estimate', ...
                     'DisturbanceModel','none', ...
                     'SearchMethod','lsqnonlin', ...
                     'Display','on');

mkModel = @(outSel,delay) setDelay( ...
    setPars(idgrey(@pitchODE,par0,'c',{outSel,USE_LAG},0, ...
                   'InputName','tau_y_cmd'), parMin,parMax,parFree), delay);

%% ---------------- Estimation ----------------
fprintf('\n--- pitch model (output = theta) ---\n');
mPitch = greyest(zE_pitch, mkModel('pitch',delayPitch), opt);

fprintf('\n--- pitch-rate model (output = theta_dot) ---\n');
mRate  = greyest(zE_rate, mkModel('rate',delayRate), opt);

fprintf('\n--- combined 2-output model ---\n');
% The two outputs have different units and very different magnitudes
% (rad vs rad/s), so weight them by their own variance - otherwise the
% rate channel dominates the cost and the pitch fit collapses.
optBoth = opt;
optBoth.OutputWeight = diag([1/var(zE_both.y(:,1)), 1/var(zE_both.y(:,2))]);
mBoth  = greyest(zE_both, mkModel('both',0), optBoth);

showPars(mPitch,'pitch model');
showPars(mRate ,'rate model');
showPars(mBoth ,'2-output model');

%% ---------------- Validation on the second run ----------------
figure('Name','Validation - pitch');    compare(zV_pitch, mPitch);  grid on
figure('Name','Validation - rate');     compare(zV_rate , mRate );  grid on
figure('Name','Validation - 2-output'); compare(zV_both , mBoth );  grid on

figure('Name','Residuals - pitch');     resid(zV_pitch, mPitch);
figure('Name','Residuals - rate');      resid(zV_rate , mRate );

[~,fitT] = compare(zV_pitch,mPitch);   [~,fitQ] = compare(zV_rate,mRate);
fprintf('\nValidation NRMSE fit: pitch %.1f %%, rate %.1f %%\n',fitT,fitQ);

% Self-check: a large gap between estimation and validation fit means you are
% fitting run-to-run variation rather than dynamics.
[~,fitTe] = compare(zE_pitch,mPitch);
fprintf('Estimation NRMSE fit (pitch): %.1f %%\n',fitTe);

%% ---------------- Frequency-domain cross-check ----------------
G = spa(zV_pitch);
figure('Name','Frequency response');
ssBoth = ss(mBoth);              % idgrey does not allow channel subreferencing
bode(ss(mPitch),'b',ssBoth(1,1),'g--',G,'r'); grid on
legend('grey-box (pitch)','grey-box (2-output)','spa (data)');
xline(2*pi*fBand(1),'k:'); xline(2*pi*fBand(2),'k:');
title('Trust the model only between the dotted lines');

%% ---------------- Physical interpretation ----------------
p  = getpvec(mPitch);
a1 = p(1); a0 = p(2); b = p(3);
wn = sqrt(max(a0,0));
fprintf('\nRig natural frequency  wn = %.2f rad/s (%.2f Hz)\n',wn,wn/2/pi);
fprintf('Rig damping ratio    zeta = %.3f\n',a1/(2*wn));
fprintf('Control effectiveness   b = %.2f rad/s^2 per unit command\n',b);

% Implied trim angle: if the rig behaves like a pendulum with a constant
% disturbance torque, theta0 = -(that torque)/a0. Compare the measured trim
% (-4..-5 deg) against what a0 would need to produce it - a useful sanity
% check on whether the offset is a spring/bias effect or something else.
theta0_meas = mean(yE);
fprintf('Measured trim = %.2f deg -> implied bias torque = %.3f rad/s^2\n', ...
        rad2deg(theta0_meas), -a0*theta0_meas);

% -------- Recovering Iy (this data alone is NOT enough) --------
%   (a) Thrust-stand route: log thrust T vs PWM for one motor, fit
%       T = c2*pwm^2 + c1*pwm + c0. The mixer maps tau_y_cmd -> pwm with the
%       front/rear pairs differenced (check the scaling in your own logs), then
%           ku = l * ( T(pwm_high) - T(pwm_low) )   per unit u,  and  Iy = ku/b.
%   (b) Inertia route: get Iy from CAD or a bifilar-pendulum swing test of the
%       whole gimbal assembly, then  ku = b*Iy.
%
% Sanity check on a0: with motors OFF, deflect the rig in pitch and let it
% swing. It should oscillate at wn/2pi Hz about the -4..-5 deg trim.
Iy = NaN;   % [kg*m^2] from route (b), if you have it
if ~isnan(Iy)
    fprintf('Implied ku = %.4f N*m per unit command\n', b*Iy);
end

%% ---------------- Model for the FREE-FLYING aircraft ----------------
% On the gimbal the plant is  b/(s^2 + a1*s + a0) [times the actuator lag].
% In free flight the pendulum stiffness a0 disappears and the rig friction a1
% is replaced by (much smaller) aerodynamic damping, so what carries over to
% controller design is b and Tm, not the whole transfer function.
s = tf('s');
if USE_LAG
    Tm = p(4);
    G_flight = b/(s^2) * 1/(Tm*s+1);   % theta(s)/u(s)
else
    G_flight = b/(s^2);
end
G_flight  %#ok<NOPTS>
G_flight_rate = s*G_flight;            % theta_dot(s)/u(s)
save('pitch_id_results.mat','mPitch','mRate','mBoth','G_flight','G_flight_rate');

%% ================= local functions =================
function [u,y,yp,t] = loadLog(fname,Ts)
    T  = readtable(fname);
    tr = T.time - T.time(1);
    t  = (0:Ts:floor(tr(end)/Ts)*Ts).';
    u  = interp1(tr,T.tau_y_cmd       ,t,'linear');
    y  = interp1(tr,T.pitch_rad       ,t,'linear');
    yp = interp1(tr,T.pitchspeed_rad_s,t,'linear');
end

function m = setDelay(m,d)
    if d > 0, m.InputDelay = d; end
end

function m = setPars(m,pmin,pmax,pfree)
    for k = 1:numel(m.Structure.Parameters)
        m.Structure.Parameters(k).Minimum = pmin(k);
        m.Structure.Parameters(k).Maximum = pmax(k);
        m.Structure.Parameters(k).Free    = pfree(k);
    end
end

function showPars(m,name)
    p  = getpvec(m);
    nm = {'a1','a0','b','Tm'};
    fprintf('\n%s:\n',name);
    for k = 1:numel(p)
        fprintf('   %-3s = %8.4f\n',nm{k},p(k));
    end
    fprintf('   free-run fit on estimation data: %s %%\n', ...
            num2str(m.Report.Fit.FitPercent(:).','%.1f  '));
end
