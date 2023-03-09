%% Manual run of REM-NREM model

cd 'L:\Lab_JamesR\lachlanW\REMswitch\REMswitch_2023\REM-NREM-extra_con_drive\'

%% Specify simulation settings

T = 3*24; % max time (h) to simulate
plotting = 1; % 1=yes for plotting figures
x0 = [0.6,-6.7,-5,-5,14.6,0.35,-0.86,-0.59,5]'; % initial conditions (Vm,Vv,Vr,Vn,Hv,n,x,xc,Hr)

%% Define defaut model parameter values
% Values from Skeldon et al. (2017)

sph = 3600; % seconds per hour (for unit conversions)

% sleep/wake switch parameters
p{1} = sph/10; % 1/tau_v=1/tau_m
p{2} = -1.8/sph;%-2.7/sph; % numv -1.8
p{3} = 1.3;%3.6;%1.3; % )  % orexin
p{4} = -2.1/sph;%-2.0/sph; % -3.50/sph; % nuvm -2.0
p{5} = 1.0; % nuvh
p{6} = -3.37; %-1.0; %-3.37; % nuvc
p{7} = -10.2; % D_0
p{8} = 1.0/45.0; %1.0/18.0; %1.0/45.0; % 1/chi_v
p{9} = 4.20/sph; %8.20/sph; %4.20/sph; % mu % mu_vm

% light processing parameters
p{10} = 0.16; % alpha0
p{11} = 9500; % I0
p{12} = 0.6; % p
p{13} = 0.4; %1.0; %0.4; % b
p{14} = 19.9; % G
p{15} = 60; % lambda
p{16} = 0.013; % beta

% circadian pacemaker parameters
tauc = 24.2; % Define intrinsic circadian period (h)
p{17} = pi/12; %1/kappa
p{18} = 0.23; %gamma
p{19} = 4.0/3.0; % coefficient for x^3 term
p{20} = (24.0/(0.99669*tauc))^2; % tau term
p{21} = 0.55; % k

% REM switch
p{22} = -1.0/sph; % nu_rn
p{23} = 2.0/sph; % nu_nv     2.0
p{24} = -1.0/sph; % nu_nr
p{25} = sph/10; % 1/tau_r=1/tau_n
p{26} = -2.0/sph; % nu_rm
p{27} = 0.5/sph; % nu_nm

% REM drive
p{28} = 1.0/sph; %1.10/sph; % mu % mu_rn
p{29} = 1.0/0.3; %1.0/0.5; % 1/chi_r6
p{30} = 1.0; % nu_rhr



% addition to sleep drive from REM 
p{31} = 0.00/sph; % mu % mu_vr

p{32} = 0.5/sph; % mu_rm 
p{33} = 0.0/sph; % mu_rv

% additional inputs into 
p{34} = -5; % DR_0
p{35} = 5; % DN_0


%% Solve the differential equations
[~,Y] = ode23s(@Phillips_REM_dblswitch_ecd,[0,2*7*24],x0',odeset,p); % This uses 2 weeks to remove transients (generally don't ever need longer than 4 weeks, unless conditions are really messed up)
Y=Y'; % Transpose of Y
disp('Finished transient')

% Remove transient and solve again
x0 = Y(:,end);
[t,Y] = ode23s(@Phillips_REM_dblswitch_ecd,[0,T],x0',odeset,p);
Y=Y'; % Transpose of Y
disp('Finished simulation')

t = t-24*floor(t(1)/24); % Rezero time

% Define a new time variable that goes a week longer for comparison of midsleeps with upcoming work times
tshift = [t(1):0.01:t(end)+168];

%% Define output variables

Vm = Y(1,:); % Voltage of MA population
Vv = Y(2,:); % Voltage of VLPO population
Vr = Y(3,:); % Voltage of REM-Brainstem(GABA) SLD/PC population
Vn = Y(4,:); % Voltage of NREM-Brainstem(GABA) vlPAG/LPT population
H = Y(5,:); % Homeostatic drive to VLPO
n = Y(6,:); % Photoreceptors
x = Y(7,:); % Pacemaker variable 1
xc = Y(8,:); % Pacemaker variable 2
Hr = Y(9,:); % REM Homeostatic drive 

C = 0.5*(1+0.80*xc-0.47*x); % Circadian drive
Dv = p{5}*H + p{6}*C + p{7}; % Overall drive to VLPO

%% states

state_n = repelem(1,length(Vm)); %m
state_n(Vv > Vm & Vr > Vn) = 2; %r
state_n(Vv > Vm & Vn >= Vr) = 3; %n
state_lab = repelem("W",length(Vm)); %m
state_lab(Vv > Vm & Vr > Vn) = "R"; %r
state_lab(Vv > Vm & Vn >= Vr) = "NR"; %n

primarystate_n = repelem(1,length(Vm)); %m
primarystate_n(Vv > Vm) = 2; %v
primarystate_lab = repelem("W",length(Vm)); %m
primarystate_lab(Vv > Vm) = "S"; %v

secondarystate_n = repelem(1,length(Vr)); %r
secondarystate_n(Vn >= Vr) = 2; %n
secondarystate_lab = repelem("R",length(Vr)); %r
secondarystate_lab(Vn >= Vr) = "NR"; %n

%% summarise sleep state cycle

sleepstate = Vm > Vv;
sleepon = t(diff(sleepstate)==-1);
sleepoff = t(diff(sleepstate)==1);

% Make sure sleep onsets and offsets are paired correctly from the end
% each pair should be onset and offset of a sleep bout
if sleepoff(1)<sleepon(1)
    sleepoff = sleepoff(2:end);
else
end

if sleepon(end)>sleepoff(end)
    sleepon = sleepon(1:end-1);
else
end
% this may need fixing in the many bout space

bout_lengths = sleepoff-sleepon;
avg_sleep_bout_lengths = mean(sleepoff-sleepon);

%% REM bouts

REM = Vr > Vn & Vv > Vm;

REMon = t(diff(REM)==1);
REMoff = t(diff(REM)==-1);

% Make sure sleep onsets and offsets are paired correctly from the end
% each pair should be onset and offset of a sleep bout
if REMoff(1)<REMon(1)
    REMoff = REMoff(2:end);
else
end

if REMon(end)>REMoff(end)
    REMon = REMon(1:end-1);
else
end
% this may need fixing in the many bout space
REM_bout_lengths = REMoff-REMon;
avg_REM_bout_lengths = mean(REMoff-REMon);

% rough estimate at cycle length - will need to be updated
cyclelengthsREM = diff(REMon);
cyclelengthsREM = cyclelengthsREM(cyclelengthsREM < avg_sleep_bout_lengths);
mincyclelengthsREM = min(cyclelengthsREM);
maxcyclelengthsREM = max(cyclelengthsREM);

%% NREM bouts

NREM = Vn > Vr & Vv > Vm;

NREMon = t(diff(NREM)==1);
NREMoff = t(diff(NREM)==-1);

% Make sure sleep onsets and offsets are paired correctly from the end
% each pair should be onset and offset of a sleep bout
if NREMoff(1)<NREMon(1)
    NREMoff = NREMoff(2:end);
else
end

if NREMon(end)>NREMoff(end)
    NREMon = NREMon(1:end-1);
else
end
% this will probably need fixing in the many bout space
NREM_bout_lengths = NREMoff-NREMon;
avg_NREM_bout_lengths = mean(NREMoff-NREMon);

cyclelengthsNREM = diff(NREMon);
cyclelengthsNREM = cyclelengthsNREM(cyclelengthsNREM < avg_sleep_bout_lengths);
mincyclelengthsNREM = min(cyclelengthsNREM);
maxcyclelengthsNREM = max(cyclelengthsNREM);

%% Proportion of each state

% using interpolation
interp_times = 0:0.005:(max(t)-0.005);
state_n_interp = round(interp1(t, double(state_n),interp_times));

prop_w = sum(state_n_interp == 1)/length(state_n_interp);
prop_r = sum(state_n_interp == 2)/length(state_n_interp);
prop_n = sum(state_n_interp == 3)/length(state_n_interp);

% amount of sleep that is REM or NREM using times

wb_NREM = sum(NREMoff(NREMoff >= sleepon(1) & NREMoff <= sleepoff(end)) - NREMon(NREMon >= sleepon(1) & NREMon <= sleepoff(end)));
wb_REM  = sum( REMoff(REMoff  >= sleepon(1) & REMoff  <= sleepoff(end)) -  REMon(REMon  >= sleepon(1) & REMon  <= sleepoff(end)));
wb_sleep = sum(sleepoff - sleepon);

prop_REM = wb_REM/wb_sleep;
prop_NREM = wb_NREM/wb_sleep;

%% length of first state 

firststate = state_n(find(ismember(t', sleepon))+1);
firststatedesc = zeros(length(firststate),2);
for slpon = 1:length(sleepon)
    minREM = min(REMoff(REMoff > sleepon(slpon)) - sleepon(slpon));
    minNREM = min(REMon(REMon > sleepon(slpon)) - sleepon(slpon));
    if minREM < minNREM && firststate(slpon) == 2
        firststatedesc(slpon,1) = 2;
        firststatedesc(slpon,2) = minREM;
    elseif minNREM < minREM && firststate(slpon) == 3
        firststatedesc(slpon,1) = 3;
        firststatedesc(slpon,2) = minNREM;
    else
        firststatedesc(slpon,1) = 99;
        firststatedesc(slpon,2) = 999;
    end 
end


%% quick plots

figure;
subplot(4,1,1)
plot(t,Y(1,:),"r",t,Y(2,:),"b");xlabel("Time");ylabel("Resting Potential");
title("Sleep-Wake Switch");xticks(0:6:72);legend("wake","sleep")
subplot(4,1,2);
plot(t,Y(3,:),"y",t,Y(4,:),"g");xlabel("Time");ylabel("Resting Potential");
title("REM-NREM Switch");xticks(0:6:72);legend("REM","NREM")
subplot(4,1,3);
plot(t,state_n);set(gca,'YDir','reverse');yticks([1,2,3]);yticklabels({'Wake', 'REM','NREM'});
xlabel("Time");ylabel("Sleep State");title("Hypnogram (Sleep State over Time)");xticks(0:6:72)
subplot(4,1,4);
%plot(t,Y(5,:),'Color',"#7E2F8E",t,Y(9,:),"m");xlabel("Time");ylabel("Homeostatic Drives");
p = plot(t,Y(5,:),t,Y(9,:)/10,"m");xlabel("Time");ylabel("Homeostatic Drives");
p(1).Color = "#7E2F8E";
title("REM-NREM Switch");xticks(0:6:72);legend("H[v] (Sleep)","H[r]/10 (REM)")
            