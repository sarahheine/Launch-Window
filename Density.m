function[Ro] = Density (h, Fx, Fb, d, t, Lat, Kp)

% Thermospheric Density (TD88) model from 
% L. Sehnal, Bull. Astron. Inst. Czechosl. 39 (1988), 120-127

% Micro-X Launch Window Calculations
% Function Density
% Version 1.0
    
% The formulas for the g values were changed such that they give
% periodicity for annual and diurnal cycles (it wasn't clear in the paper 
% what the argument of the sin function is, and I chose the most logical
% way, i.e. added the /365*2*pi and /24*2*pi terms).

% The values obtained with this function agree very well with the graphs in
% the above paper.

% Only valid between 150 and 750 km!!
% The exosphere > 750 km is not taken into account.


% These are the input parameters
% h ... Altitude in km
% Fx ... Solar radio flux at 10.7 cm for previous day (in 1E-22 Wm^-2Hz^-1) (60 to 220)
% Fb ... Mean solar radio flux at 10.7 cm for 81 days (in 1E-22 Wm^-2Hz^-1) (60 to 220)
% d ... Day in the year (0 to 365)
% t ... Local solar time in hours (0 to 24)
% Lat ... Latitude (-pi to +pi)
% Kp ... Geomagnetic index before date (0 to 10)


% Model parameters

H = 40;        % Scale height 40 km
hbar = 120;    % Reference altitude in km

a = [0.007 0.2875 0.04762 0.0471 7 7 0.3333 15];

p = [0 0 263 -263 -29.41 8.0913 10.0813];

Kn0 = [-2.1041E-12 -7.8628E-13 1.4228E-14 1.9794E-13 2.2388E-13 -2.1152E-12 -1.7174E-13];
Kn1 = [2.5051E-9 1.1295E-10 5.2364E-11 -7.6532E-11 -1.1571E-10 9.0702E-10 7.6048E-11];
Kn2 = [-1.5026E-10 -6.8497E-11 1.3776E-11 1.5305E-11 1.6575E-11 -2.2008E-10 -1.8880E-11];
Kn3 = [7.9498E-11 6.3063E-11 -3.2733E-12 -7.6109E-12 -9.8514E-12 9.7056E-11 8.3009E-12];


% Calculate f-values and k-value
fm = (Fb - 60) / 160;
fx = 1 + a(1) * (Fx - Fb);
f0 = a(2) + fm;
k0 = 1 + a(3)*(Kp-3);


% Calculate g-values
g = [1 (fm)+a(4) sin((d-p(3))/365*2*pi)*sin(Lat) (a(5)*fm+1)*sin((d-p(4))/365*2*pi) (a(6)*fm+1)*sin((2*(d-p(5)))/365*2*pi) (a(7)*fm+1)*sin((t-p(6))/24*2*pi)*cos(Lat) (a(8)*fm+1)*sin((2*(t-p(7)))/24*2*pi)*(cos(Lat)^2)];


% Calculate h-values
hn = Kn0 + Kn1*(exp((hbar - h)/H)) + Kn2*(exp((hbar - h)/(2*H))) + Kn3*(exp((hbar - h)/(3*H)));


% Calculate Density
Ro = fx * f0 * k0 * (g(1)*hn(1) + g(2)*hn(2) + g(3)*hn(3) + g(4)*hn(4) + g(5)*hn(5) + g(6)*hn(6) + g(7)*hn(7));          % Unit is kg/m3

