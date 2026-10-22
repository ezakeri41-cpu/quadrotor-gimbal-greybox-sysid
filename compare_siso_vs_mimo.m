%% COMPARE_SISO_VS_MIMO  Three SISO models (stacked) vs mBothMIMO_coupled on the MIMO chirp runs
%  SISO model  : roll, pitch and yaw mBoth models from the earlier SISO identification,
%                put side by side (tau_x->roll, tau_y->pitch, tau_z->yaw, no coupling,
%                with their own actuator lag / delay if they had one).
%  MIMO model  : mBothMIMO_coupled (full 3x3 B, no lag).
%  Validation  : the 6 SECOND-run files (not used to fit either model).
%  Pre-processing (same for both models): resample to 50 Hz, remove mean, zero-phase
%  band-pass BAND (stored with the MIMO model) on inputs and outputs.
%  Needs: mimo_full_models.mat, roll/pitch/yaw_id_results.mat, mimoFull_ODE.m,
%         rollODE.m, pitchODE.m, yawODE.m, the CSVs below.
clear; clc; close all;

valFiles = {'mimo_chirp_sysid_log2.csv'
            'mimo_chirp_sysid_log(1-1.5,0.3-0.8,1.8-2.3)1.csv'
            'mimo_chirp_sysid_log(0.3-0.8,1.8-2.3,1-1.5)1.csv'
            'mimo_chirp_sysid_log(1-1.5,1.8-2.3,0.3-0.8)1.csv'
            'mimo_chirp_sysid_log(1.8-2.3,1-1.5,0.3-0.8)1.csv'
            'mimo_chirp_sysid_log(1.8-2.3,0.3-0.8,1-1.5)1.csv'};
runTag   = {'x0.3-0.8 y1.0-1.5 z1.8-2.3'
            'x1.0-1.5 y0.3-0.8 z1.8-2.3'
            'x0.3-0.8 y1.8-2.3 z1.0-1.5'
            'x1.0-1.5 y1.8-2.3 z0.3-0.8'
            'x1.8-2.3 y1.0-1.5 z0.3-0.8'
            'x1.8-2.3 y0.3-0.8 z1.0-1.5'};

Ts = 0.02;
yN = {'roll','pitch','yaw','rollrate','pitchrate','yawrate'};
uN = {'tau_x','tau_y','tau_z'};

%% ---------------- models ----------------
R    = load('mimo_full_models.mat');  BAND = R.BAND;
mMIMO = R.mBothMIMO_coupled;

sR = load('roll_id_results.mat');  sP = load('pitch_id_results.mat');  sY = load('yaw_id_results.mat');
S  = append(ss(sR.mBoth), ss(sP.mBoth), ss(sY.mBoth));  % outputs: roll rollrate pitch pitchrate yaw yawrate
S  = S([1 3 5 2 4 6],:);                                  % -> roll pitch yaw rollrate pitchrate yawrate
S.InputName = uN;  S.OutputName = yN;
mSISO = idss(S);

pR = getpvec(sR.mBoth);  pP = getpvec(sP.mBoth);  pY = getpvec(sY.mBoth);  pM = getpvec(mMIMO);
fprintf('===== diagonal parameters: SISO vs MIMO coupled =====\n');
fprintf('%-6s %8s %8s %8s   %8s %8s %8s\n','axis','a1 SISO','a0 SISO','b SISO','a1 MIMO','a0 MIMO','b MIMO');
fprintf('%-6s %8.2f %8.2f %8.2f   %8.2f %8.2f %8.2f\n','roll' ,pR(1:3),pM(1),pM(2),pM(7));
fprintf('%-6s %8.2f %8.2f %8.2f   %8.2f %8.2f %8.2f\n','pitch',pP(1:3),pM(3),pM(4),pM(11));
fprintf('%-6s %8.2f %8.2f %8.2f   %8.2f %8.2f %8.2f\n','yaw'  ,pY(1:3),pM(5),pM(6),pM(15));

%% ---------------- validation ----------------
[zz,pp,kk] = butter(2,BAND/(0.5/Ts),'bandpass');  [sos,g] = zp2sos(zz,pp,kk);
optE = compareOptions('InitialCondition','e');
nR = numel(valFiles);  FS = zeros(6,nR);  FM = zeros(6,nR);
for i = 1:nR
    T  = readtable(valFiles{i},'VariableNamingRule','preserve');
    tr = T.time - T.time(1);
    t  = (0:Ts:floor(tr(end)/Ts)*Ts).';
    U  = [interp1(tr,T.tau_x_cmd,t), interp1(tr,T.tau_y_cmd,t), interp1(tr,T.tau_z_cmd,t)];
    Y  = [interp1(tr,T.roll_rad,t),  interp1(tr,T.pitch_rad,t), interp1(tr,T.yaw_rad,t), ...
          interp1(tr,T.rollspeed_rad_s,t), interp1(tr,T.pitchspeed_rad_s,t), interp1(tr,T.yawspeed_rad_s,t)];
    U  = filtfilt(sos,g,U - mean(U));
    Y  = filtfilt(sos,g,Y - mean(Y));
    z  = iddata(Y,U,Ts,'InputName',uN,'OutputName',yN);

    hS = compare(z,mSISO,optE);  hS = hS.OutputData;
    hM = compare(z,mMIMO,optE);  hM = hM.OutputData;
    FS(:,i) = 100*(1 - vecnorm(Y-hS)./vecnorm(Y-mean(Y))).';
    FM(:,i) = 100*(1 - vecnorm(Y-hM)./vecnorm(Y-mean(Y))).';

    figure('Name',sprintf('SISO vs MIMO, validation run %d: %s',i,runTag{i}));
    for j = 1:6
        subplot(6,1,j);
        plot(t,Y(:,j),'Color',[.6 .6 .6]); hold on;
        plot(t,hS(:,j),'b--'); plot(t,hM(:,j),'r','LineWidth',1);
        grid on; ylabel(yN{j});
        legend('measured',sprintf('SISO (%.1f%%)',FS(j,i)), ...
               sprintf('mBothMIMO\\_coupled (%.1f%%)',FM(j,i)),'Location','eastoutside');
    end
    xlabel('Time [s]');
    sgtitle(sprintf('Validation run %d (%s): SISO vs mBothMIMO\\_coupled',i,runTag{i}));
end

fprintf('\n===== Validation fit [%%], second-run files (SISO | MIMO coupled) =====\n');
fprintf('%-28s','run');  fprintf('%17s',yN{:});  fprintf('\n');
for i = 1:nR
    fprintf('%-28s',runTag{i});
    for j = 1:6, fprintf('%8.1f |%6.1f ',FS(j,i),FM(j,i)); end
    fprintf('\n');
end
fprintf('%-28s','MEAN');
for j = 1:6, fprintf('%8.1f |%6.1f ',mean(FS(j,:)),mean(FM(j,:))); end
fprintf('\n%-28s','MIMO - SISO');
for j = 1:6, fprintf('%17.1f',mean(FM(j,:))-mean(FS(j,:))); end
fprintf('\n\nOverall mean fit: SISO = %.1f %%, mBothMIMO_coupled = %.1f %%\n',mean(FS(:)),mean(FM(:)));

figure('Name','Mean validation fit: SISO vs MIMO');
bar(categorical(yN,yN),[mean(FS,2) mean(FM,2)]);
grid on; ylabel('Mean fit over 6 validation runs [%]');
legend('SISO models (stacked)','mBothMIMO\_coupled','Location','best');
title('SISO vs MIMO coupled on the second-run MIMO chirp files');
