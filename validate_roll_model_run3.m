%% VALIDATE_ROLL_MODEL_RUN3  Validate the identified mRoll model on a
%  THIRD, independent chirp run (never used in estimation or the first
%  validation).
%
%  Requires:
%    - roll_id_results.mat  (produced by roll_greybox_id.m, contains mRoll)
%    - standalone_chirp_roll_log(0.1amp0.2nbt3).csv  in the current folder

clear; clc; close all;

%% ---------------- Settings ----------------
valFile = 'mimo_chirp_sysid_log1.csv';
Ts      = 0.02;   % must match the Ts used when mRoll was identified

%% ---------------- Load the identified model ----------------
load('roll_id_results.mat','mRoll');   % also has mRate, mBoth, G_flight

%% ---------------- Load + resample the new run ----------------
[u3,y3,yp3,t3] = loadLog(valFile,Ts);

z3_roll = iddata(y3,u3,Ts,'OutputName','roll','OutputUnit','rad', ...
                          'InputName','tau_x_cmd','InputUnit','-');
z3_roll = detrend(z3_roll,0);   % remove this run's own mean bias

%% ---------------- Validate ----------------
figure('Name','Validation on run 3 - roll');
compare(z3_roll, mRoll);
grid on
title('mRoll vs. standalone chirp run 3 (unseen data)');

figure('Name','Residuals on run 3 - roll');
resid(z3_roll, mRoll);

[~,fit3] = compare(z3_roll, mRoll);
fprintf('mRoll validation NRMSE fit on run 3: %.1f %%\n', fit3);

%% ---------------- Local function ----------------
function [u,y,yp,t] = loadLog(fname,Ts)
    T  = readtable(fname);
    tr = T.time - T.time(1);
    t  = (0:Ts:floor(tr(end)/Ts)*Ts).';
    u  = interp1(tr,T.tau_x_cmd      ,t,'linear');
    y  = interp1(tr,T.roll_rad       ,t,'linear');
    yp = interp1(tr,T.rollspeed_rad_s,t,'linear');
end
