function Q = sigmoid(V)

% Sigmoidal firing rate function, output in units of per hour

Qmax = 100*3600;
theta = 10;
sigma = 3;

Q = Qmax./(1+exp((theta-V)/sigma));

end