function dS = dSdV(V)


Qmax = 100*3600;
theta = 10;
sigma = 3;

dS = (Qmax*exp((theta-V)/sigma))/(sigma *(exp((theta-V)/sigma)+1)^2);


end