%% MIMO_FULL_GREYBOX_ID  Full-coupling vs no-coupling MIMO grey-box on 6 rotated chirp runs
%  Estimation : the 6 FIRST-run files  (merged as one multi-experiment dataset)
%  Validation : the 6 SECOND-run files (same inputs, never used in estimation)
%  Both models use angles AND rates, no actuator lag.
%
%    mBothMIMO_coupled : full 3x3 B  (bxx bxy bxz / byx byy byz / bzx bzy bzz)
%    mBothMIMO         : diagonal B  (bxx byy bzz only)
%
%  Model (mimoFull_ODE.m), b<axis><input>:
%    phi_dd   = -a0r*phi   - a1r*p + bxx*u_x + bxy*u_y + bxz*u_z
%    theta_dd = -a0p*theta - a1p*q + byx*u_x + byy*u_y + byz*u_z
%    psi_dd   = -a0y*psi   - a1y*r + bzx*u_x + bzy*u_y + bzz*u_z
%
%  Pre-processing: resample to 50 Hz, zero-phase band-pass 0.25-2.4 Hz on ALL inputs and
%  outputs (removes yaw drift / offsets and noise outside the excited range; the same
%  filter on u and y keeps the input-output relation unchanged).
clear; clc; close all;

%% ---------------- settings ----------------
Ts   = 0.02;
BAND = [0.25 2.4];                  % [Hz] band-pass on all signals
SISO_DIR = pwd;                     % optional SISO results for start values

estFiles = {'mimo_chirp_sysid_log1.csv'
            'mimo_chirp_sysid_log(1-1.5,0.3-0.8,1.8-2.3).csv'
            'mimo_chirp_sysid_log(0.3-0.8,1.8-2.3,1-1.5).csv'
            'mimo_chirp_sysid_log(1-1.5,1.8-2.3,0.3-0.8).csv'
            'mimo_chirp_sysid_log(1.8-2.3,1-1.5,0.3-0.8).csv'
            'mimo_chirp_sysid_log(1.8-2.3,0.3-0.8,1-1.5).csv'};
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

yN = {'roll','pitch','yaw','rollrate','pitchrate','yawrate'};
uN = {'tau_x','tau_y','tau_z'};

%% ---------------- data ----------------
[zz,pp,kk] = butter(2,BAND/(0.5/Ts),'bandpass');  [sos,g] = zp2sos(zz,pp,kk);
nR = numel(estFiles);  zE = cell(1,nR);  zV = cell(1,nR);
for i = 1:nR
    zE{i} = loadRun(estFiles{i},Ts,sos,g,uN,yN);
    zV{i} = loadRun(valFiles{i},Ts,sos,g,uN,yN);
end
zEst = merge(zE{:});                              % 6-experiment estimation set
fprintf('Loaded %d estimation runs and %d validation runs (band-pass %.2f-%.2f Hz).\n',nR,nR,BAND);

%% ---------------- parameters & start values ----------------
names = {'a1r','a0r','a1p','a0p','a1y','a0y', ...
         'bxx','bxy','bxz','byx','byy','byz','bzx','bzy','bzz'};
pr = loadSiso(SISO_DIR,'roll_id_results.mat' ,{'mBoth','mRoll'} ,[2.6 18  16 ]);
pp_= loadSiso(SISO_DIR,'pitch_id_results.mat',{'mBoth','mPitch'},[2.2 30  35 ]);
py = loadSiso(SISO_DIR,'yaw_id_results.mat'  ,{'mBoth','mYaw'}  ,[2.5 0.9 7.3]);
p0 = [pr(1:2) pp_(1:2) py(1:2)  pr(3) 0 0  0 pp_(3) 0  0 0 py(3)];
iDiag = [7 11 15];  iOff = [8 9 10 12 13 14];

freeUnc = true(1,15);  freeUnc(iOff) = false;     % mBothMIMO
freeCpl = true(1,15);                             % mBothMIMO_coupled

opt = greyestOptions('InitialState','estimate','DisturbanceModel','none', ...
                     'SearchMethod','lsqnonlin','Display','off');
Yall = vertcat(zEst.OutputData{:});
opt.OutputWeight = diag(1./var(Yall));            % equal weight to every output

%% ---------------- estimation ----------------
fprintf('Estimating mBothMIMO (no coupling) ...\n');
mBothMIMO = greyest(zEst,buildModel(names,p0,freeUnc,uN,yN),opt);
fprintf('Estimating mBothMIMO_coupled (full B) ...\n');
p0c = getpvec(mBothMIMO).';                       % start from the uncoupled solution
mBothMIMO_coupled = greyest(zEst,buildModel(names,p0c,freeCpl,uN,yN),opt);
mBothMIMO.Name = 'mBothMIMO';  mBothMIMO_coupled.Name = 'mBothMIMO_coupled';

%% ---------------- parameters ----------------
[pU,sU] = getpvec(mBothMIMO);  [pC,sC] = getpvec(mBothMIMO_coupled);
fprintf('\n===== identified parameters =====\n');
fprintf('%-5s %22s   %22s %8s\n','','mBothMIMO','mBothMIMO_coupled','t (cpl)');
for i = 1:15
    if freeUnc(i), su = sprintf('%9.3f +/- %6.3f',pU(i),sU(i)); else, su = '0 (fixed)'; end
    fprintf('%-5s %22s   %9.3f +/- %6.3f %8.1f\n',names{i},su,pC(i),sC(i),pC(i)/sC(i));
end
Bc = reshape(pC(7:15),3,3).';                     % rows: roll pitch yaw, cols: tau_x tau_y tau_z
fprintf('\nCoupled B matrix  (rows roll/pitch/yaw, cols tau_x/tau_y/tau_z):\n');
disp(Bc);
fprintf('Cross terms relative to the diagonal term of the same axis [%%]:\n');
fprintf('%-6s %8s %8s %8s\n','','tau_x','tau_y','tau_z');
ax = {'roll','pitch','yaw'};
for r = 1:3, fprintf('%-6s %8.1f %8.1f %8.1f\n',ax{r},100*Bc(r,:)/Bc(r,r)); end
fprintf('(|t| < ~3 or a few %% of the diagonal: treat the cross term as not supported by the data)\n');

%% ---------------- validation on the second-run files ----------------
optE = compareOptions('InitialCondition','e');
FU = zeros(6,nR);  FC = zeros(6,nR);
for i = 1:nR
    z  = zV{i};  y = z.OutputData;
    hU = compare(z,mBothMIMO,optE);          hU = hU.OutputData;
    hC = compare(z,mBothMIMO_coupled,optE);  hC = hC.OutputData;
    FU(:,i) = nrmse(y,hU);  FC(:,i) = nrmse(y,hC);

    t = (0:size(y,1)-1)'*Ts;
    figure('Name',sprintf('Validation run %d: %s',i,runTag{i}));
    for j = 1:6
        subplot(6,1,j);
        plot(t,y(:,j),'Color',[.6 .6 .6]); hold on;
        plot(t,hU(:,j),'b'); plot(t,hC(:,j),'r');
        grid on; ylabel(yN{j});
        legend('measured', sprintf('mBothMIMO (%.1f%%)',FU(j,i)), ...
               sprintf('mBothMIMO\\_coupled (%.1f%%)',FC(j,i)),'Location','eastoutside');
    end
    xlabel('Time [s]');
    sgtitle(sprintf('Validation run %d (%s), band-pass %.2f-%.2f Hz',i,runTag{i},BAND));
end

fprintf('\n===== Validation fit [%%], second-run files (uncoupled | coupled) =====\n');
fprintf('%-28s','run');  fprintf('%17s',yN{:});  fprintf('\n');
for i = 1:nR
    fprintf('%-28s',runTag{i});
    for j = 1:6, fprintf('%8.1f |%6.1f ',FU(j,i),FC(j,i)); end
    fprintf('\n');
end
fprintf('%-28s','MEAN');
for j = 1:6, fprintf('%8.1f |%6.1f ',mean(FU(j,:)),mean(FC(j,:))); end
fprintf('\n%-28s','coupled - uncoupled');
for j = 1:6, fprintf('%17.1f',mean(FC(j,:))-mean(FU(j,:))); end
fprintf('\n\nOverall mean fit: mBothMIMO = %.1f %%, mBothMIMO_coupled = %.1f %%  (difference %+.1f)\n', ...
        mean(FU(:)),mean(FC(:)),mean(FC(:))-mean(FU(:)));

%% ---------------- summary figure ----------------
figure('Name','Mean validation fit per output');
bar(categorical(yN,yN),[mean(FU,2) mean(FC,2)]);
grid on; ylabel('Mean fit over 6 validation runs [%]');
legend('mBothMIMO (no coupling)','mBothMIMO\_coupled (full B)','Location','best');
title('Validation on second-run files');

save('mimo_full_models.mat','mBothMIMO','mBothMIMO_coupled','names','BAND','FU','FC','runTag');

%% ================= local functions =================
function z = loadRun(fname,Ts,sos,g,uN,yN)
    T  = readtable(fname,'VariableNamingRule','preserve');
    tr = T.time - T.time(1);
    t  = (0:Ts:floor(tr(end)/Ts)*Ts).';
    U  = [interp1(tr,T.tau_x_cmd,t), interp1(tr,T.tau_y_cmd,t), interp1(tr,T.tau_z_cmd,t)];
    Y  = [interp1(tr,T.roll_rad,t),  interp1(tr,T.pitch_rad,t), interp1(tr,T.yaw_rad,t), ...
          interp1(tr,T.rollspeed_rad_s,t), interp1(tr,T.pitchspeed_rad_s,t), interp1(tr,T.yawspeed_rad_s,t)];
    U  = filtfilt(sos,g,U - mean(U));
    Y  = filtfilt(sos,g,Y - mean(Y));
    z  = iddata(Y,U,Ts,'InputName',uN,'OutputName',yN,'ExperimentName',fname);
end

function m = buildModel(names,p0,free,uN,yN)
    par = [names(:) num2cell(p0(:))];
    m = idgrey(@mimoFull_ODE,par,'c',{},0,'InputName',uN,'OutputName',yN);
    for k = 1:numel(names)
        P = m.Structure.Parameters(k);
        P.Free = free(k);
        if k <= 6, P.Minimum = 0; end             % a0, a1 >= 0
        m.Structure.Parameters(k) = P;
    end
end

function p = loadSiso(dirName,file,vars,def)
    p = def;
    path = fullfile(dirName,file);
    if ~exist(path,'file'), path = which(file); end
    if isempty(path), fprintf('%s not found - default start values used.\n',file); return; end
    s = load(path);
    for k = 1:numel(vars)
        if isfield(s,vars{k})
            v = getpvec(s.(vars{k})).';  n = min(3,numel(v));  p(1:n) = v(1:n);
            fprintf('%-22s (%s): a1=%.3f a0=%.3f b=%.3f\n',file,vars{k},p);  return
        end
    end
end

function f = nrmse(y,s)
    f = 100*(1 - vecnorm(y-s)./vecnorm(y-mean(y))).';
end
