%% MIMO_FULL_VALIDATE_ROLL  Validate mBothMIMO_coupled (and mBothMIMO) on roll-only runs
%  Data: tau_x only (tau_y = tau_z = 0). Logged outputs: roll, pitch, yaw, rollrate
%  (no pitchrate / yawrate in these files), so only these 4 outputs are scored.
%  Needs: mimo_full_models.mat (from mimo_full_greybox_id.m), mimoFull_ODE.m, the CSVs below.
%  Same pre-processing as the identification: resample to 50 Hz, remove mean,
%  zero-phase band-pass (BAND stored with the models) on inputs and outputs.
clear; clc; close all;

files = {'mimo_sine_validation_log3.csv'
         'mimo_sine_validation_log4.csv'
         'mimo_sine_validation_log5.csv'};
tags  = {'mimo_sine_validation_log3','mimo_sine_validation_log4','mimo_sine_validation_log5'};

R  = load('mimo_full_models.mat');
Ts = 0.02;  BAND = R.BAND;
yAll = {'roll','pitch','yaw','rollrate','pitchrate','yawrate'};
oIdx = 1:4;                          % outputs present in the files
yN   = yAll(oIdx);
uN   = {'tau_x','tau_y','tau_z'};

% 4-output versions of the identified models (same parameters)
mC = ss(R.mBothMIMO_coupled);  mC = mC(oIdx,:);
mU = ss(R.mBothMIMO);          mU = mU(oIdx,:);

[zz,pp,kk] = butter(2,BAND/(0.5/Ts),'bandpass');  [sos,g] = zp2sos(zz,pp,kk);
optE = compareOptions('InitialCondition','e');

nF = numel(files);  FC = zeros(4,nF);  FU = zeros(4,nF);
for i = 1:nF
    T  = readtable(files{i},'VariableNamingRule','preserve');
    tr = T.time - T.time(1);
    t  = (0:Ts:floor(tr(end)/Ts)*Ts).';
    U  = [interp1(tr,T.tau_x_cmd,t), zeros(numel(t),2)];        % tau_y = tau_z = 0
    Y  = [interp1(tr,T.roll_rad,t), interp1(tr,T.pitch_rad,t), ...
          interp1(tr,T.yaw_rad,t),  interp1(tr,T.rollspeed_rad_s,t)];
    U  = filtfilt(sos,g,U - mean(U));
    Y  = filtfilt(sos,g,Y - mean(Y));
    z  = iddata(Y,U,Ts,'InputName',uN,'OutputName',yN);

    hC = compare(z,mC,optE);  hC = hC.OutputData;
    hU = compare(z,mU,optE);  hU = hU.OutputData;
    FC(:,i) = 100*(1 - vecnorm(Y-hC)./vecnorm(Y-mean(Y))).';
    FU(:,i) = 100*(1 - vecnorm(Y-hU)./vecnorm(Y-mean(Y))).';

    figure('Name',['Roll validation: ' tags{i}]);
    for j = 1:4
        subplot(4,1,j);
        plot(t,Y(:,j),'Color',[.6 .6 .6]); hold on;
        plot(t,hC(:,j),'r','LineWidth',1); plot(t,hU(:,j),'b--');
        grid on; ylabel(yN{j});
        legend('measured',sprintf('mBothMIMO\\_coupled (%.1f%%)',FC(j,i)), ...
               sprintf('mBothMIMO (%.1f%%)',FU(j,i)),'Location','eastoutside');
    end
    xlabel('Time [s]');
    sgtitle(sprintf('%s - tau_y = tau_z = 0 (band-pass %.2f-%.2f Hz)',tags{i},BAND),'Interpreter','none');
end

fprintf('\n===== Fit [%%] on roll-only runs: mBothMIMO_coupled | mBothMIMO =====\n');
fprintf('%-30s','run'); fprintf('%16s',yN{:}); fprintf('\n');
for i = 1:nF
    fprintf('%-30s',tags{i});
    for j = 1:4, fprintf('%8.1f |%6.1f ',FC(j,i),FU(j,i)); end
    fprintf('\n');
end
fprintf('%-30s','MEAN');
for j = 1:4, fprintf('%8.1f |%6.1f ',mean(FC(j,:)),mean(FU(j,:))); end
fprintf('\n');
