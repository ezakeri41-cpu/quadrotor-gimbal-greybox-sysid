%% ROLL_GREYBOX_ID  Grey-box system identification of a 3-DOF gimbal roll axis
%  Estimation : standalone_chirp_roll_log.csv
%  Validation : standalone_chirp_roll_log_0_1amp0_2nbt2_.csv
%  Requires   : System Identification Toolbox + rollODE.m on the path
%
%  Log columns: time, tau_x_cmd, roll_rad, pitch_rad, yaw_rad,
%               rollspeed_rad_s, pwm1..pwm4
%
%  MODEL (see rollODE.m):
%      phi_ddot = -a1*phi_dot - a0*phi + b*u        [+ first-order actuator lag]
%  The a0 (stiffness) and a1 (damping) terms are NOT in the textbook
%  quadrotor model phi_ddot = tau_x/Ix. They belong to the TEST RIG
%  (CG offset below the pivot -> pendulum spring, bearing friction,
%  cable stiffness). They must be identified anyway, otherwise they are
%  absorbed into b and corrupt the control-effectiveness estimate.

clear; clc; close all;

%% ---------------- User settings ----------------
estFile   = 'standalone_chirp_roll_log.csv';
valFile   = 'standalone_chirp_roll_log(0.1amp0.2nbt3).csv';

Ts        = 0.02;      % [s] uniform resampling period (logs are ~49.8 Hz)
USE_LAG   = true;      % true -> 3rd-order model with ESC/rotor lag Tm
USE_DELAY = true;      % true -> run delayest and put the delay in the model
LPFILTER  = false;     % optional zero-phase low-pass of the outputs
fcut      = 8;         % [Hz] cut-off if LPFILTER = true

% Chirp excitation band (from the logged input): roughly 0.4 ... 2 Hz.
% The identified model is only trustworthy inside this band.
fBand     = [0.4 2];   % [Hz]

%% ---------------- Load + resample onto a uniform grid ----------------
% iddata assumes a constant sample time; the logs have small jitter
% (dt = 20.09 +/- 0.17 ms), so interpolate onto an exact grid first.
[uE,yE,ypE,tE] = loadLog(estFile,Ts);
[uV,yV,ypV,tV] = loadLog(valFile,Ts);

if LPFILTER
    [bf,af] = butter(4,fcut/(0.5/Ts));
    yE  = filtfilt(bf,af,yE);   ypE = filtfilt(bf,af,ypE);
    yV  = filtfilt(bf,af,yV);   ypV = filtfilt(bf,af,ypV);
end

%% ---------------- iddata objects + detrending ----------------
% Offsets are physical and run-specific (IMU mounting bias, gimbal not
% perfectly level, small trim in the command), so each record is detrended
% with its OWN means. Use detrend(...,1) instead of 0 if you see a drift.
zE_roll = iddata(yE ,uE,Ts,'OutputName','roll','OutputUnit','rad', ...
                           'InputName','tau_x_cmd','InputUnit','-');
zE_rate = iddata(ypE,uE,Ts,'OutputName','rollrate','OutputUnit','rad/s', ...
                           'InputName','tau_x_cmd');
zE_both = iddata([yE ypE],uE,Ts,'OutputName',{'roll','rollrate'});

zV_roll = iddata(yV ,uV,Ts,'OutputName','roll','InputName','tau_x_cmd');
zV_rate = iddata(ypV,uV,Ts,'OutputName','rollrate','InputName','tau_x_cmd');
zV_both = iddata([yV ypV],uV,Ts,'OutputName',{'roll','rollrate'});

zE_roll = detrend(zE_roll,0);   zV_roll = detrend(zV_roll,0);
zE_rate = detrend(zE_rate,0);   zV_rate = detrend(zV_rate,0);
zE_both = detrend(zE_both,0);   zV_both = detrend(zV_both,0);

%% ---------------- Input/output delay estimate ----------------
nk_roll = delayest(zE_roll,3,1,0,25);   % search 0..25 samples (0..0.5 s)
nk_rate = delayest(zE_rate,3,1,0,25);
fprintf('delayest: roll %d samples (%.0f ms), rate %d samples (%.0f ms)\n', ...
        nk_roll,1e3*nk_roll*Ts,nk_rate,1e3*nk_rate*Ts);

% With USE_LAG = true the ESC/rotor lag is already a model state, so a
% separate transport delay is usually unnecessary and the two fight each
% other. Only apply the delay if it is clearly non-zero.
if USE_DELAY && ~USE_LAG
    delayRoll = nk_roll*Ts;  delayRate = nk_rate*Ts;
else
    delayRoll = 0;           delayRate = 0;
end

%% ---------------- Grey-box model template ----------------
% Initial guesses: wn ~ 4.3 rad/s (a0 ~ 18), zeta ~ 0.3 (a1 ~ 2.6),
% b ~ 16 rad/s^2 per unit command, Tm ~ 0.15 s. These came from a
% quick equation-error fit of this exact data set.
% idgrey wants an N-by-2 cell array {name, value}; bounds and Free flags
% are set afterwards on m.Structure.Parameters.
par0 = {'a1', 2.6 ; 'a0', 18 ; 'b', 16 ; 'Tm', 0.15};

parMin  = [0   0   -Inf 0.01];
parMax  = [Inf Inf  Inf 0.30];   % Tm > ~0.3 s is not a physical ESC/rotor lag
parFree = [true true true USE_LAG];   % Tm frozen if the lag state is unused

opt = greyestOptions('InitialState','estimate', ...
                     'DisturbanceModel','none', ...
                     'SearchMethod','lsqnonlin', ...
                     'Display','on');

mkModel = @(outSel,delay) setDelay( ...
    setPars(idgrey(@rollODE,par0,'c',{outSel,USE_LAG},0, ...
                   'InputName','tau_x_cmd'), parMin,parMax,parFree), delay);

%% ---------------- Estimation ----------------
fprintf('\n--- roll model (output = phi) ---\n');
mRoll = greyest(zE_roll, mkModel('roll',delayRoll), opt);

fprintf('\n--- roll-rate model (output = phi_dot) ---\n');
mRate = greyest(zE_rate, mkModel('rate',delayRate), opt);

fprintf('\n--- combined 2-output model ---\n');
% The two outputs have different units and very different magnitudes
% (rad vs rad/s), so weight them by their own variance - otherwise the
% rate channel dominates the cost and the roll fit collapses.
optBoth = opt;
optBoth.OutputWeight = diag([1/var(zE_both.y(:,1)), 1/var(zE_both.y(:,2))]);
mBoth = greyest(zE_both, mkModel('both',0), optBoth);

showPars(mRoll,'roll model');
showPars(mRate,'rate model');
showPars(mBoth,'2-output model');

%% ---------------- Validation on the second run ----------------
% mBoth has 2 outputs, so it must be compared against the 2-output record.
figure('Name','Validation - roll');  compare(zV_roll, mRoll);  grid on
figure('Name','Validation - rate');  compare(zV_rate, mRate);  grid on
figure('Name','Validation - 2-output'); compare(zV_both, mBoth); grid on

figure('Name','Residuals - roll');   resid(zV_roll, mRoll);
figure('Name','Residuals - rate');   resid(zV_rate, mRate);

[~,fitR] = compare(zV_roll,mRoll);   [~,fitP] = compare(zV_rate,mRate);
fprintf('\nValidation NRMSE fit: roll %.1f %%, rate %.1f %%\n',fitR,fitP);

% Self-check: is the fit on the estimation set much better than on the
% validation set? If yes you are over-fitting run-to-run variation.
[~,fitRe] = compare(zE_roll,mRoll);
fprintf('Estimation NRMSE fit (roll): %.1f %%\n',fitRe);

%% ---------------- Frequency-domain cross-check ----------------
% Non-parametric estimate from the same data, compared with the model.
G = spa(zV_roll);
figure('Name','Frequency response');
% idgrey does not allow channel subreferencing -> convert to ss first.
ssBoth = ss(mBoth);
bode(ss(mRoll),'b',ssBoth(1,1),'g--',G,'r'); grid on
legend('grey-box (roll)','grey-box (2-output)','spa (data)');
xline(2*pi*fBand(1),'k:'); xline(2*pi*fBand(2),'k:');
title('Trust the model only between the dotted lines');

%% ---------------- Physical interpretation ----------------
p  = getpvec(mRoll);
a1 = p(1); a0 = p(2); b = p(3);
wn = sqrt(max(a0,0));
fprintf('\nRig natural frequency  wn = %.2f rad/s (%.2f Hz)\n',wn,wn/2/pi);
fprintf('Rig damping ratio    zeta = %.3f\n',a1/(2*wn));
fprintf('Control effectiveness   b = %.2f rad/s^2 per unit command\n',b);

% -------- Recovering Ix (this data alone is NOT enough) --------
% Identified: b = ku/Ix, where ku [N*m] is the physical torque produced by
% one unit of tau_x_cmd. Two unknowns, one equation. Break the tie with ONE
% independent measurement:
%
%   (a) Thrust-stand route: log thrust T vs PWM for one motor, fit
%       T = c2*pwm^2 + c1*pwm + c0. In your logs the mixer maps
%       tau_x_cmd -> pwm as pwm1,4 = 1200 - 707*u and pwm2,3 = 1200 + 707*u
%       (check the 707 against your own logs). Then
%           ku = l * ( T(1200+707) - T(1200-707) )      per unit u,
%       with l the motor arm's roll moment arm, and  Ix = ku/b.
%
%   (b) Inertia route: get Ix from CAD or a bifilar-pendulum swing test of
%       the whole gimbal assembly, then  ku = b*Ix.
%
% Sanity check on a0: with motors OFF, deflect the rig and let it swing.
% It should oscillate at wn/2pi ~ 0.7 Hz. If it does, a0 is the rig
% pendulum, not a modelling artefact.
l  = NaN;   % [m] roll moment arm
Ix = NaN;   % [kg*m^2] from route (b), if you have it
if ~isnan(Ix)
    fprintf('Implied ku = %.4f N*m per unit command\n', b*Ix);
end

%% ---------------- Model for the FREE-FLYING aircraft ----------------
% On the gimbal the plant is  b/(s^2 + a1*s + a0)  [times the actuator lag].
% In free flight the pendulum stiffness a0 disappears and the rig friction
% a1 is replaced by (much smaller) aerodynamic damping. So the quantity you
% actually carry over to controller design is b (and Tm), not the whole TF:
s = tf('s');
if USE_LAG
    Tm = p(4);
    G_flight = b/(s^2) * 1/(Tm*s+1);
else
    G_flight = b/(s^2);
end
G_flight  %#ok<NOPTS>   % phi(s)/u(s) predicted for free flight
save('roll_id_results.mat','mRoll','mRate','mBoth','G_flight');

%% ================= local functions =================
function [u,y,yp,t] = loadLog(fname,Ts)
    T  = readtable(fname);
    tr = T.time - T.time(1);
    t  = (0:Ts:floor(tr(end)/Ts)*Ts).';
    u  = interp1(tr,T.tau_x_cmd     ,t,'linear');
    y  = interp1(tr,T.roll_rad      ,t,'linear');
    yp = interp1(tr,T.rollspeed_rad_s,t,'linear');
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
    p = getpvec(m);
    nm = {'a1','a0','b','Tm'};
    fprintf('\n%s:\n',name);
    for k = 1:numel(p)
        fprintf('   %-3s = %8.4f\n',nm{k},p(k));
    end
    fprintf('   free-run fit on estimation data: %s %%\n', ...
            num2str(m.Report.Fit.FitPercent(:).','%.1f  '));
end
