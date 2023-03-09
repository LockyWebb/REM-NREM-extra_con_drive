function dy = Phillips_REM_dblswitch_ecd(t,y,p)
% Define differential equations for combined REM-NREM sleep-circadian model.
% Equations originally from Skeldon et al. (2017)
% updated to include REM-NREM switching

% Variables: y(1) = Vm, y(2) = Vv, y(3) = Vr, y(4) = Vn, y(5) = Hv, 
% y(6) = n, y(7) = x, y(8) = xc, y(9) = Hr

% define circadian input to VLPO
%C = 0.5*(1+0.80*y(6)-0.47*y(5));
C = 0.5*(1+0.80*y(8)-0.47*y(7)); % possible update for REM
% Define total drive to VLPO (combination of circadian and homeostatic)
%Dv = p{5}*y(3) + p{6}*C + p{7};

I = light_func(t).*(y(1)>y(2)); % define light level (filtering by sleep/wake state)

alpha = p{10}*(I/p{11}).^p{12}; % define alpha (retinal activation rate)
B = (1-p{13}*y(7)).*(1-p{13}*y(8))*p{14}*alpha*(1-y(6)); % define B (retinal drive to pacemaker)

Qmax = 100*3600;

% Define the differential equations:
% some coefficnets in p are set to 0 to remove terms
dy = [p{1}*(-y(1) + (p{2}*sigmoid(y(2)) + p{3}  ));%Vm
    % 1/tau * (-Vm + nu_mv*Q_v + O)
    
    p{1}*(-y(2) + p{4}*sigmoid(y(1)) + p{5}*y(5) + p{6}*C + p{7});%Vv
    %1/tau * (-Vv + nu_vm*Q_m + D[: nu_vh*H + nu_vc*C+ D_0])
    
    p{25}*(-y(3) + p{22}*sigmoid(y(4)) + p{26}*sigmoid(y(1)) + p{30}*y(9) + p{34}); % Vr
    %1/tau * (-Vr + nu_rn*Q_n + nu_rm*Q_m + nu_rhr*Hr + DR_0)
    
    p{25}*(-y(4) + p{23}*sigmoid(y(2)) + p{24}*sigmoid(y(3)) + p{27}*sigmoid(y(1)) + p{35}); % Vn
    %1/tau * (-Vn + nu_nv*Q_v + nu_nr*Q_r + nu_nm*Q_m + DN_0)
    
    p{8}*(-y(5) + p{9}*sigmoid(y(1)) + p{31}*sigmoid(y(3)));%Hv
    %1/chi_v(-H_v + mu_vm*Q_m + mu_vr*Q_r)    
    
    %
    p{15}*(alpha*(1-y(6))-p{16}*y(6));%n
    %lambda*(alpha*(1-n)-beta*n)
    
    p{17}*(p{18}*(y(7)-p{19}*y(7).^3)-y(8)*(p{20}+p{21}*B));%x
    %(1/kappa)*(gamma*(x-4/3*x^3)-y*(tau+k*B))
    
    p{17}*(y(7)+B);%xc or y
    %(1/kappa)*(x+B)
    
    p{29}*(-y(9) + p{28}*sigmoid(y(4)) + p{32}*sigmoid(y(1)) + p{33}*sigmoid(y(2)))];%Hr
    %1/chi_r(-H_r + mu_rn*Q_n + mu_rm*Q_m + mu_rv*Q_v)    
end