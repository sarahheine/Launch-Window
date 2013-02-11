function[Fx, Fb] = SolarFlux(Year, Month, Day)

% 10.7 cm solar radio flux taken from http://solarscience.msfc.nasa.gov/predict.shtml
% This data is valid from 1995 to 2015.

% Micro-X Launch Window Calculations
% Function Solar Flux
% Version 1.0
    
% The data is given in variable FluxData as a function of time ( Variable
% MonthData). Monthdata is the number of months which have passed since
% January 1997.

MonthData = [1 : 1 : 216];
FluxData = [93.5 93.5 93.6 94 94.6 95.5 96.6 98 99.7 101.6 103.9 106.3 109.1 112 115.2 118.5 122 125.7 129.4 133.3 137.2 141.2 145.1 149.1 153 156.9 160.8 164.6 168.2 171.8 175.2 178.6 181.7 184.8 187.6 190.3 192.9 195.2 197.4 199.4 201.2 202.8 204.2 205.4 206.5 207.3 208 208.4 208.7 208.8 208.8 208.5 208.1 207.6 206.8 206 204.9 203.8 202.5 201.1 199.5 197.9 196.1 194.2 192.3 190.3 188.2 186 183.7 181.4 179.1 176.7 174.3 171.9 169.4 166.9 164.4 162 159.5 157 154.6 152.1 149.7 147.4 145 142.7 140.5 138.3 136.1 134 131.9 129.9 128 126.1 124.3 122.5 120.8 119.1 117.5 116 114.6 113.2 111.8 110.5 109.3 108.2 107 106 105 104 103.1 102.3 101.5 100.8 100 99.4 98.8 98.2 97.6 97.1 96.6 96.2 95.8 95.4 95 94.7 94.4 94.1 93.8 93.6 93.3 93.1 93.1 93.2 93.6 94.4 95.5 97.2 99.3 101.9 105 108.6 112.6 117.1 121.9 127.1 132.6 138.2 144.1 150.1 156.1 162.2 168.2 174.2 180.1 185.8 191.3 196.7 201.8 206.7 211.2 215.5 219.5 223.2 226.5 229.5 232.1 234.4 236.4 238 239.2 240.1 240.7 240.9 240.9 240.5 239.8 238.8 237.5 236 234.2 232.1 229.8 227.4 224.7 221.8 218.8 215.7 212.4 208.9 205.4 201.8 198.2 194.5 190.7 186.9 183.2 179.4 175.6 171.9 168.2 164.5 160.9 157.4 153.9 150.6 147.3 144.1 141 138.1 135.2 132.4 129.8 127.2 124.8 122.5];

ActualDay = 12*(Year - 1997) + Month + Day/30.0;         % in months
AverageDay = [ActualDay : -1/30 : ActualDay - (81/30)];  % in months

Fx = interp1(MonthData,FluxData,ActualDay);              % Interpolate to get value for specific day

AverageFlux = interp1(MonthData,FluxData,AverageDay);    % Interpolate and calculate 81 day average
Fb = mean(AverageFlux);

