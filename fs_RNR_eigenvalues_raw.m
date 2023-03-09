function ev = fs_RNR_eigenvalues_raw(fsp,Vrp,Vnp)

% sph = 3600; % seconds per hour (for unit conversions)
% fsp{1} = sph/10; % 1/tau_r=1/tau_n
% fsp{2} = -2.1/sph; % nu_rn
% fsp{3} = -1.8/sph; % nu_nr 
% fsp{4} = 16; % A
% fsp{5} = 0.90/sph; % B
% fsp{6} = -0.4/sph; % nu_rHr    
% fsp{7} = -0.5/sph; % Hr     

pm = fsp{1}*sqrt(fsp{2}*fsp{3}*dSdV(Vrp)*dSdV(Vnp));

if pm > 0
    ev = [-fsp{1}+pm -fsp{1}-pm];
elseif pm == 0
    ev = [-fsp{1} NA];
elseif pm < 0
    ev = [999 999];
end 


% ev = roots([1 lin con]);


end