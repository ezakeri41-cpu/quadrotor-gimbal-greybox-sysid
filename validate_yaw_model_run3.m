%% VALIDATE_YAW_MODEL_RUN3  Validate the identified mYaw model on a THIRD,
%  independent chirp run (never used in estimation or the first validation).
%
%  Requires:
%    - yaw_id_results.mat  (produced by yaw_greybox_id.m, contains mYaw)
%    - standalone_chirp_yaw_log3.csv  in the current folder
%
%  NOTE: rename valFile below if your third yaw log uses a different name
%  (your roll/pitch files used "...nbt3..."; your yaw files so far are
%  numbered log1/log2, so log3 is assumed here).

clear; clc; close all;

%% ---------------- Settings ----------------
valFile = 'mimo_chirp_sysid_log1.csv';
Ts      = 0.02;      % must match the Ts used when mYaw was identified
DETREND_ORDER = 1;   % must match the DETREND_ORDER used in yaw_greybox_id.m
                      % (yaw drifts - linear detrend, unlike roll/pitch's 0)
TRIM_TIME = 0;        % must match the TRIM_TIME used in yaw_greybox_id.m
                      % (0 = whole record; set e.g. 12 if you trimmed there)
LPFILTER  = true;    % must match yaw_greybox_id.m
fcut      = 5;

%% ---------------- Load the identified model ----------------
load('yaw_id_results.mat','mYaw');   % also has mRate, mBoth, G_flight, ...

%% ---------------- Load + resample the new run ----------------
[u3,y3,yp3,t3] = loadLog(valFile,Ts);

fprintf('Run 3: mean yaw %.2f deg, last 5 s mean %.2f deg\n', ...
        rad2deg(mean(y3)), rad2deg(mean(y3(t3 > t3(end)-5))));

if TRIM_TIME > 0
    k3 = t3 <= TRIM_TIME;
    t3 = t3(k3); u3 = u3(k3); y3 = y3(k3); yp3 = yp3(k3);
    fprintf('Using only the first %.0f s of this record.\n',TRIM_TIME);
end

if LPFILTER
    [bf,af] = butter(4,fcut/(0.5/Ts));
    y3 = filtfilt(bf,af,y3);
    yp3 = filtfilt(bf,af,yp3); %#ok<NASGU> (kept for consistency; yaw uses angle output)
end

z3_yaw = iddata(y3,u3,Ts,'OutputName','yaw','OutputUnit','rad', ...
                         'InputName','tau_z_cmd','InputUnit','-');
z3_yaw = detrend(z3_yaw,DETREND_ORDER);   % remove this run's own drift

%% ---------------- Validate ----------------
figure('Name','Validation on run 3 - yaw');
compare(z3_yaw, mYaw);
grid on
title('mYaw vs. standalone chirp run 3 (unseen data)');

figure('Name','Residuals on run 3 - yaw');
resid(z3_yaw, mYaw);

[~,fit3] = compare(z3_yaw, mYaw);
fprintf('mYaw validation NRMSE fit on run 3: %.1f %%\n', fit3);
fprintf('Remember: only the first ~5-8 s of each yaw run (roughly 0.3-0.6 Hz)\n');
fprintf('has usable SNR - a modest fit number here is expected, not a red flag.\n');

%% ---------------- Local function ----------------
function [u,y,yp,t] = loadLog(fname,Ts)
    T  = readtable(fname);
    tr = T.time - T.time(1);
    t  = (0:Ts:floor(tr(end)/Ts)*Ts).';
    u  = interp1(tr,T.tau_z_cmd     ,t,'linear');
    y  = interp1(tr,T.yaw_rad       ,t,'linear');
    yp = interp1(tr,T.yawspeed_rad_s,t,'linear');
end
