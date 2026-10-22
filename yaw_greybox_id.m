%% YAW_GREYBOX_ID  Grey-box system identification of a 3-DOF gimbal yaw axis
%  Estimation : standalone_chirp_yaw_log1.csv
%  Validation : standalone_chirp_yaw_log2.csv
%  Requires   : System Identification Toolbox + yawODE.m on the path
%
%  Log columns: time, tau_z_cmd, roll_rad, pitch_rad, yaw_rad,
%               yawspeed_rad_s, pwm1..pwm4
%
%  MODEL (see yawODE.m):
%      psi_ddot = -a1*psi_dot - a0*psi + b*u   [+ optional actuator lag]
%
%  READ THIS BEFORE TRUSTING THE RESULT
%  ------------------------------------
%  The yaw data is usable but clearly WEAKER than roll/pitch. Measured from
%  your own logs (peak-to-peak yaw vs. yaw-rate noise floor, per 5 s window):
%
%      0- 5 s : 19.6 deg   rate std 15.5 deg/s    <- good SNR
%      5-10 s :  5.6 deg   rate std  8.7 deg/s    <- marginal
%     10-15 s :  3.9 deg   rate std  6.0 deg/s    <- poor
%     15-30 s : 2.3-4.2 deg  rate std 4.1-5.0 deg/s  <- essentially noise
%
%  So only roughly the first 5-8 s of the sweep (about 0.3-0.6 Hz) carries
%  real information. Above ~1 Hz the yaw response is buried in gyro noise.
%  The script therefore offers TRIM_TIME to optionally fit on the informative
%  part only, and the fit statistics over the full record should be read with
%  that in mind. Report the identified model as valid over ~0.3-0.8 Hz, NOT
%  over the whole 0.3-2 Hz sweep.

clear; clc; close all;

%% ---------------- User settings ----------------
estFile   = 'standalone_chirp_yaw_log1.csv';
valFile   = 'standalone_chirp_yaw_log2.csv';

Ts        = 0.02;      % [s] uniform resampling period (logs are ~49.7 Hz)
USE_LAG   = false;     % Tm is NOT identifiable on this axis - leave false
USE_DELAY = true;      % true -> run delayest and put the delay in the model
DETREND_ORDER = 1;     % 1 = linear detrend. Yaw is gyro-integrated and DRIFTS
                       % (run1: -0.04 -> -2.15 deg, run2: 7.70 -> 8.53 deg),
                       % so order 1 is the sensible default here, unlike
                       % roll/pitch where order 0 sufficed.
TRIM_TIME = 0;         % 0 = use the whole record (as requested).
                       % Set e.g. 12 to fit only the first 12 s, where the
                       % signal-to-noise ratio is actually good. Compare both.
LPFILTER  = true;      % yaw rate is noticeably noisier than roll/pitch
fcut      = 5;         % [Hz] cut-off if LPFILTER = true

% Excitation band of the chirp: 0.3 ... 2 Hz, but see the SNR note above:
% the yaw model is only credible up to roughly 0.8 Hz.
fBand     = [0.3 0.8]; % [Hz] band over which to trust the model

%% ---------------- Load + resample onto a uniform grid ----------------
[uE,yE,ypE,tE] = loadLog(estFile,Ts);
[uV,yV,ypV,tV] = loadLog(valFile,Ts);

fprintf('Estimation run: mean yaw %.2f deg, last 5 s mean %.2f deg\n', ...
        rad2deg(mean(yE)), rad2deg(mean(yE(tE > tE(end)-5))));
fprintf('Validation run: mean yaw %.2f deg, last 5 s mean %.2f deg\n\n', ...
        rad2deg(mean(yV)), rad2deg(mean(yV(tV > tV(end)-5))));

if TRIM_TIME > 0
    kE = tE <= TRIM_TIME;  tE = tE(kE); uE = uE(kE); yE = yE(kE); ypE = ypE(kE);
    kV = tV <= TRIM_TIME;  tV = tV(kV); uV = uV(kV); yV = yV(kV); ypV = ypV(kV);
    fprintf('Using only the first %.0f s of each record.\n\n',TRIM_TIME);
end

if LPFILTER
    [bf,af] = butter(4,fcut/(0.5/Ts));
    yE  = filtfilt(bf,af,yE);   ypE = filtfilt(bf,af,ypE);
    yV  = filtfilt(bf,af,yV);   ypV = filtfilt(bf,af,ypV);
end

%% ---------------- iddata objects + detrending ----------------
zE_yaw  = iddata(yE ,uE,Ts,'OutputName','yaw','OutputUnit','rad', ...
                           'InputName','tau_z_cmd','InputUnit','-');
zE_rate = iddata(ypE,uE,Ts,'OutputName','yawrate','OutputUnit','rad/s', ...
                           'InputName','tau_z_cmd');
zE_both = iddata([yE ypE],uE,Ts,'OutputName',{'yaw','yawrate'});

zV_yaw  = iddata(yV ,uV,Ts,'OutputName','yaw','InputName','tau_z_cmd');
zV_rate = iddata(ypV,uV,Ts,'OutputName','yawrate','InputName','tau_z_cmd');
zV_both = iddata([yV ypV],uV,Ts,'OutputName',{'yaw','yawrate'});

zE_yaw  = detrend(zE_yaw ,DETREND_ORDER);  zV_yaw  = detrend(zV_yaw ,DETREND_ORDER);
zE_rate = detrend(zE_rate,DETREND_ORDER);  zV_rate = detrend(zV_rate,DETREND_ORDER);
zE_both = detrend(zE_both,DETREND_ORDER);  zV_both = detrend(zV_both,DETREND_ORDER);

%% ---------------- Input/output delay estimate ----------------
nk_yaw  = delayest(zE_yaw ,3,1,0,25);
nk_rate = delayest(zE_rate,3,1,0,25);
fprintf('delayest: yaw %d samples (%.0f ms), rate %d samples (%.0f ms)\n', ...
        nk_yaw,1e3*nk_yaw*Ts,nk_rate,1e3*nk_rate*Ts);
fprintf('(Usually resonance/noise phase rather than true dead time.)\n');

if USE_DELAY && ~USE_LAG
    delayYaw = nk_yaw*Ts;  delayRate = nk_rate*Ts;
else
    delayYaw = 0;          delayRate = 0;
end

%% ---------------- Grey-box model template ----------------
% Initial guesses from a quick simulation-error fit of this exact data set:
% a1 ~ 2.5, a0 ~ 0.9, b ~ 7.3. Note how much smaller b is than on roll (24.8)
% and pitch (43.4) - yaw torque comes from rotor drag, not thrust difference.
par0 = {'a1', 2.5 ; 'a0', 0.9 ; 'b', 7.3 ; 'Tm', 0.05};

parMin  = [0   0   -Inf 0.01];
parMax  = [Inf Inf  Inf 0.30];
parFree = [true true true USE_LAG];   % Tm frozen unless USE_LAG is on

opt = greyestOptions('InitialState','estimate', ...
                     'DisturbanceModel','none', ...
                     'SearchMethod','lsqnonlin', ...
                     'Display','on');

mkModel = @(outSel,delay) setDelay( ...
    setPars(idgrey(@yawODE,par0,'c',{outSel,USE_LAG},0, ...
                   'InputName','tau_z_cmd'), parMin,parMax,parFree), delay);

%% ---------------- Estimation ----------------
fprintf('\n--- yaw model (output = psi) ---\n');
mYaw  = greyest(zE_yaw, mkModel('yaw',delayYaw), opt);

fprintf('\n--- yaw-rate model (output = psi_dot) ---\n');
mRate = greyest(zE_rate, mkModel('rate',delayRate), opt);

fprintf('\n--- combined 2-output model ---\n');
optBoth = opt;
optBoth.OutputWeight = diag([1/var(zE_both.y(:,1)), 1/var(zE_both.y(:,2))]);
mBoth = greyest(zE_both, mkModel('both',0), optBoth);

showPars(mYaw ,'yaw model');
showPars(mRate,'rate model');
showPars(mBoth,'2-output model');

%% ---------------- Validation on the second run ----------------
figure('Name','Validation - yaw');      compare(zV_yaw , mYaw );  grid on
figure('Name','Validation - rate');     compare(zV_rate, mRate);  grid on
figure('Name','Validation - 2-output'); compare(zV_both, mBoth);  grid on

figure('Name','Residuals - yaw');       resid(zV_yaw , mYaw );
figure('Name','Residuals - rate');      resid(zV_rate, mRate);

[~,fitY] = compare(zV_yaw ,mYaw );   [~,fitR] = compare(zV_rate,mRate);
fprintf('\nValidation NRMSE fit: yaw %.1f %%, rate %.1f %%\n',fitY,fitR);

[~,fitYe] = compare(zE_yaw,mYaw);
fprintf('Estimation NRMSE fit (yaw): %.1f %%\n',fitYe);
fprintf('Expect the RATE fit to be poor here - the yaw gyro noise floor\n');
fprintf('(4-6 deg/s) is comparable to the response itself after ~10 s.\n');

%% ---------------- Frequency-domain cross-check ----------------
G = spa(zV_yaw);
figure('Name','Frequency response');
ssBoth = ss(mBoth);
bode(ss(mYaw),'b',ssBoth(1,1),'g--',G,'r'); grid on
legend('grey-box (yaw)','grey-box (2-output)','spa (data)');
xline(2*pi*fBand(1),'k:'); xline(2*pi*fBand(2),'k:');
title('Yaw: trust only between the dotted lines (~0.3-0.8 Hz)');

%% ---------------- Physical interpretation ----------------
p  = getpvec(mYaw);
a1 = p(1); a0 = p(2); b = p(3);
wn = sqrt(max(a0,0));
fprintf('\nTorsional natural frequency wn = %.2f rad/s (%.2f Hz)\n',wn,wn/2/pi);
if wn/2/pi < fBand(1)
    fprintf('  ** wn is BELOW the excited band - a0 is only weakly identified.\n');
    fprintf('     Report it as an order-of-magnitude estimate, not a firm value.\n');
end
fprintf('Damping ratio             zeta = %.3f\n',a1/(2*max(wn,eps)));
fprintf('Control effectiveness        b = %.2f rad/s^2 per unit command\n',b);
fprintf('  (compare: roll b ~ 24.8, pitch b ~ 43.4 - yaw authority is much lower,\n');
fprintf('   which is physically expected since yaw torque comes from rotor drag.)\n');

% -------- Recovering Iz (this data alone is NOT enough) --------
%   b = ku/Iz only. To split them:
%   (a) Measure the rotor drag torque constant kM on a stand (torque vs. rpm),
%       combine with the mixer scaling from tau_z_cmd to the four PWMs, then
%       Iz = ku/b.
%   (b) Get Iz from CAD or a bifilar-pendulum test, then ku = b*Iz.
%
%   Cross-check on a0: with the motors OFF, twist the rig in yaw and release.
%   If it slowly returns, that is the wiring harness acting as a torsion
%   spring; if it just stays where you put it, a0 is essentially zero and the
%   small fitted value is an artifact of the initial transient.
Iz = NaN;   % [kg*m^2] from route (b), if you have it
if ~isnan(Iz)
    fprintf('Implied ku = %.4f N*m per unit command\n', b*Iz);
end

%% ---------------- Model for the FREE-FLYING aircraft ----------------
% Unlike roll/pitch, there is no gravity restoring term about the vertical
% axis, so in free flight a0 -> 0 and your textbook model psi_ddot = tau_z/Iz
% is recovered directly. What carries over is b (and the damping a1, part of
% which is genuine rotor drag rather than a rig artifact).
s = tf('s');
if USE_LAG
    Tm = p(4);
    G_flight = b/(s^2) * 1/(Tm*s+1);   % psi(s)/u(s)
else
    G_flight = b/(s^2);                % psi(s)/u(s)
end
G_flight  %#ok<NOPTS>
G_flight_rate = s*G_flight;            % psi_dot(s)/u(s)

% A more honest free-flight yaw model keeps the aerodynamic drag damping,
% since rotor drag does act in flight (only the harness stiffness disappears):
G_flight_damped = b/(s*(s + a1));      %#ok<NASGU>

save('yaw_id_results.mat','mYaw','mRate','mBoth', ...
     'G_flight','G_flight_rate','G_flight_damped');

%% ================= local functions =================
function [u,y,yp,t] = loadLog(fname,Ts)
    T  = readtable(fname);
    tr = T.time - T.time(1);
    t  = (0:Ts:floor(tr(end)/Ts)*Ts).';
    u  = interp1(tr,T.tau_z_cmd     ,t,'linear');
    y  = interp1(tr,T.yaw_rad       ,t,'linear');
    yp = interp1(tr,T.yawspeed_rad_s,t,'linear');
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
