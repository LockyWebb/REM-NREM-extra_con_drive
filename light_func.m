function L = light_func(t)
%% Light input as a function of time (hours) in units of lux
% From Andrew Phillips, Equations from Skeldon et al. (2017)

dshift = 0; % Amount by which to shift the solar curve (hours)
evelight = 0;%40; % Set level (lux) of evening light (home lighting during wake, outside solar day) 



dawn=(8+dshift)*3600;     %Clock time when light is approx (day+evening)/2
dusk=(17+dshift)*3600;    %Clock time when light is approx (day+evening)/2
daylight=500; % Max daylight level (lux)
steepness=1/6000; % Steepness of solar curve
a=evelight;
b=(daylight-evelight)/2;
% Combine two sigmoid curves to approximate the solar curve
lightf=@(a,b,dawn,dusk,time)a+b*(tanh(steepness*(time*3600-dawn))-tanh(steepness*(time*3600-dusk)));
L = lightf(a,b,dawn,dusk,mod(t,24));

end