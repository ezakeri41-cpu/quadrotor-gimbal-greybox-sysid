%% VALIDATE_PITCH_MODEL_RUN3  Validate the identified mPitch model on a
%  THIRD, independent chirp run (never used in estimation or the first
%  validation).
%
%  Requires:
%    - pitch_id_results.mat  (produced by pitch_greybox_id.m, contains mPitch)
%    - standalone_chirp_pitch_log(0.1amp0.2nbt3).csv  in the current folder

clear; clc; close all;

%% ---------------- Settings ----------------
valFile = 'mimo_chirp_sysid_log1.csv';
Ts      = 0.02;   % must match the Ts used when mPitch was identified
DETREND_ORDER = 0; % must match the DETREND_ORDER used in pitch_greybox_id.m

%% ---------------- Load the identified model ----------------
load('pitch_id_results.mat','mPitch');   % also has mRate, mBoth, G_flight, G_flight_rate

%% ---------------- Load + resample the new run ----------------
[u3,y3,yp3,t3] = loadLog(valFile,Ts);

fprintf('Run 3: mean pitch %.2f deg, last 5 s mean %.2f deg\n', ...
        rad2deg(mean(y3)), rad2deg(mean(y3(t3 > t3(end)-5))));

z3_pitch = iddata(y3,u3,Ts,'OutputName','pitch','OutputUnit','rad', ...
                           'InputName','tau_y_cmd','InputUnit','-');
z3_pitch = detrend(z3_pitch,DETREND_ORDER);   % remove this run's own trim

%% ---------------- Validate ----------------
figure('Name','Validation on run 3 - pitch');
compare(z3_pitch, mPitch);
grid on
title('mPitch vs. standalone chirp run 3 (unseen data)');

figure('Name','Residuals on run 3 - pitch');
resid(z3_pitch, mPitch);

[~,fit3] = compare(z3_pitch, mPitch);
fprintf('mPitch validation NRMSE fit on run 3: %.1f %%\n', fit3);

%% ---------------- Local function ----------------
function [u,y,yp,t] = loadLog(fname,Ts)
    T  = readtable(fname);
    tr = T.time - T.time(1);
    t  = (0:Ts:floor(tr(end)/Ts)*Ts).';
    u  = interp1(tr,T.tau_y_cmd       ,t,'linear');
    y  = interp1(tr,T.pitch_rad       ,t,'linear');
    yp = interp1(tr,T.pitchspeed_rad_s,t,'linear');
end
