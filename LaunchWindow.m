% Micro-X Launch Window Calculations
% Version 1.0

% Enter the start day and time, and the time for which the calculation
% shall be performed here. UT is Universal Time. Alt is the rocket's 
% observation altitude in km.

% This is the start time and data

Day = 1;
Month = 1;
Year = 2013;
UT = 0.0;

% Calculation time and time step specified here:

CalcTime = 24*365;     % in hours
TimeStep = 1./60;       % in hours

% Minimum observation altitude:

Alt = 100;%160.0;            % in kilometers


% Launch Criteria

%Puppis-A
MinDistMoonPup = 30;         % in degrees
MaxAltSun = -18;             % in degrees  -12 < theta < -18 = astronomical twilight
MinAltPup = 15;               % in degrees
MaxAbsorption = 100;          % in percent

% It is impossible the predict the geomagnetic index Kp that far in
% advance. Kp ranges from 0 to 9, with 9 being a very strong magnetic
% storm. Actual Kp values are published on http://www.swpc.noaa.gov/rt_plots/kp_3d.html
% Kp = 3 seems to be a pretty reasonable assumption.

Kp = 3;   


% Location (latitude, longitude) of White Sands Missile Range (WSMR) 
% from the NASA Sounding Rocket Handbook

Lat = 32.5;                 % in degrees
Long = -106.5;              % in degrees

Lat = Lat/180*pi;           % in rad
Long = Long/180*pi;         % in rad


% Right ascension (RA) and declination (Dec) of Puppis-A (from NED)

% Variable names here are a poor choice - should be switched

% RAPup = 125.7;            % for PuppisA, in degrees (08h 22min 48s)
% DecPup = -42.76;          % for PuppisA, in degrees (-42° 45' 36")

 RAPup = 350.8583;            % for CasA, in degrees (23h 23min 26s)
 DecPup = 58.8;               % for CasA, in degrees (58° 48')

% RAPup = 187.7058;            % for M87, in degrees (12h 30min 49.4s)
% DecPup = 12.3911;            % for M87, in degrees (12° 23' 28")

RAPup = RAPup/180*pi;       % in rad
DecPup = DecPup/180*pi;     % in rad


% Here, the Julian day number is calculated from the date given above:

if (Month < 3)
    Month = Month + 12;
    Year = Year - 1;
end

A = floor(Year / 100);
B = floor(A / 4);
C = 2 - A + B;
E = floor(365.25 * (Year + 4716));
F = floor(30.6001 * (Month + 1));
DayNumber = C + Day + E + F - 1524.5;

% Add Time of Day
DayNumber = DayNumber + UT/24.0;

% Write header into File

save_filename = 'casA_launchwindow2.txt';

dlmwrite(save_filename,'Launch Criteria', 'delimiter', '')
dlmwrite(save_filename,' ','-append')

dlmwrite(save_filename,'Minimum Distance Between Moon and Pup-A (deg)','-append', 'delimiter', '')
dlmwrite(save_filename,MinDistMoonPup,'-append', 'delimiter', '')
dlmwrite(save_filename,' ','-append')

dlmwrite(save_filename,'Maximum Altitude of the Sun (deg)','-append', 'delimiter', '')
dlmwrite(save_filename,MaxAltSun,'-append', 'delimiter', '')
dlmwrite(save_filename,' ','-append')

dlmwrite(save_filename,'Minimum Altitude of Puppis A (deg)','-append', 'delimiter', '')
dlmwrite(save_filename,MinAltPup,'-append', 'delimiter', '')
dlmwrite(save_filename,' ','-append')

dlmwrite(save_filename,'Maximum Allowable Absorption (percent)','-append', 'delimiter', '')
dlmwrite(save_filename,MaxAbsorption,'-append', 'delimiter', '')
dlmwrite(save_filename,' ','-append')

dlmwrite(save_filename,'JDN;ActualDay;ActualMonth;ActualYear;Hour;Min;Sec;AltSun;AzSun;DistSunPup;AltPup;AzPup;AltMoon;AzMoon;MoonPhase;DeltaV;DistMoonPup;Absorption;LaunchWindow;', '-append', 'delimiter', '')

% ***************************** MAIN LOOP ******´************************

ActualDayNumber = 0;
Window = 0;
OldWindow = 0;
x = 10;

for ActualDayNumber = DayNumber : TimeStep/24.0 : DayNumber + CalcTime/24.0
 
    % Calculate ActualDay, ActualMonth, ActualYear from ActualDayNumber
    
      Z = floor(ActualDayNumber + 0.5);
      W = floor((Z - 1867216.25)/36524.25);
      X = floor(W/4);
      A = Z+1+W-X;
      B = A+1524;
      C = floor((B-122.1)/365.25);
      D = floor(365.25 * C);
      E = floor((B-D)/30.6001);
      F = floor(30.6001 * E);
      ActualDay = B-D-F;    
    
      ActualMonth = E-1;
      if (ActualMonth > 12)
         ActualMonth = E - 13;
      end
      
      ActualYear = C - 4716;
      if (ActualMonth < 3)
         ActualYear = C - 4715;
      end
      
      DayoftheYear = floor(ActualDayNumber) - (367*Year - floor(7*(Year+floor((10)/12))/4) + floor(275/9) + 1 - 730530);
      Time = (((ActualDayNumber) - floor(ActualDayNumber))-0.5)*24;

      if (Time < 0)
         Time = Time + 24;
      end
      
    % Results are calculated here
    % 2451545.5 is subtracted from JDN to convert into J2000.
    
    [AltSun, AzSun, RASun, DecSun] = SunPosition(ActualDayNumber-2451545.5, Lat, Long);                                            % Calculate the Sun's position in the sky at WSMR
    %[AltMoon, AzMoon, RAMoon, DecMoon, DistanceMoonkm] = MoonPosition(ActualDayNumber-2451545.5, Lat, Long, Alt);                                    % Calculate the Moon's position in the sky at WSMR
    %[AltMoon, AzMoon] = MoonPosition(ActualDayNumber-2451545.5, Lat, Long, Alt);                                    % Calculate the Moon's position in the sky at WSMR
    [RAMoon, DecMoon, AzMoon, AltMoon, DistanceMoonkm] = LunarAzEl(ActualDayNumber, Lat*180/pi, Long*180/pi, Alt);                                   % Calculate the Moon's position in the sky at WSMR
    RAMoon = RAMoon * pi / 180.;
    DecMoon = DecMoon * pi / 180.;
    AltMoon = AltMoon * pi / 180.;
    AzMoon = AzMoon * pi / 180.;
    
    
    [AltPup, AzPup] = FixPosition(ActualDayNumber-2451545.5, Lat, Long, RAPup, DecPup);                             % Calculate the Puppis-A's position in the sky at WSMR
    DistMoonPup = acos(sin(AltMoon)*sin(AltPup)+cos(AltMoon)*cos(AltPup)*cos(AzPup-AzMoon));                        % Calculate the angular distance between Pup-A and the Moon
    DistSunPup = acos(sin(AltSun)*sin(AltPup)+cos(AltSun)*cos(AltPup)*cos(AzPup-AzSun));                        % Calculate the angular distance between Pup-A and the Moon
    MoonPhase = mphase(ActualDayNumber, RASun, RAMoon, DecSun, DecMoon, DistanceMoonkm);
    [DeltaV,MoonElongation,ObjMoonDist,MoonPhase2] = moon_sky_brightness(ActualDayNumber,[RAPup,DecPup],[Long,Lat,Alt*1000.]);

    AzSun = AzSun/pi*180;
    AltSun = AltSun/pi*180;

    AzMoon = AzMoon/pi*180;
    AltMoon = AltMoon/pi*180;

    AzPup = AzPup/pi*180;
    AltPup = AltPup/pi*180;

    DistMoonPup = DistMoonPup/pi*180;
    DistSunPup = DistSunPup/pi*180;


    % Use function SolarFLux to predict the 10.7 cm solar radio flux as a
    % function of time. Fx is the solar radio flux at 10.7 cm for the previous day (in
    % 1E-22 Wm^-2Hz^-1) (should be between 60 and 220) and Fb is the mean solar
    % radio flux at 10.7 cm for 81 days (in 1E-22 Wm^-2Hz^-1) (should also be between 60 and 220)
  
    [Fx, Fb] = SolarFlux (ActualYear, ActualMonth, ActualDay);


    % Calculate Column Density (Density integrated along line of sight through
    % the atmsophere)

    % Integrate (sum and divide) the column density for a line pointing towards the zenith
    
    StepSize = 1;
    N = 0;
    for i = Alt : StepSize: 750
        D = Density(i, Fx, Fb, DayoftheYear, (Time+24*Long/360), Lat, Kp)/1000;     % Original unit is kg/m^3, converted into g/cm^3
        N = N + D * StepSize;
    end

    % Calculate Form Factor. This form factor takes into account that the
    % atmosphere is approximately spherical, and that Micro-X does not look
    % straight to zenith, but to a point which is at some angle not equal 90° 
    % over the horizon.
    
    Alpha = (6371 + Alt)/(6371 + 750)*sin(AltPup + pi/2);
    Beta = pi - Alpha - (AltPup + pi/2);
    L = (Alt + 6371)/sin(Alpha)*sin(Beta);
    FormFactor = L/(750-Alt);

    N = N * FormFactor;


    % Calculate Absorption 

    Sigma = 3310;                              % Mass attenuation coefficient of air and 1 keV photons in cm2/g - this is a maximum (from High Energy Astrophysics Handbook 13-2)
    Absorption = (1 - exp(-Sigma*N)) * 1E8;    % Absorption in 10^6 percent
  
    
    % This block compares actual values with the allowable thresholds. If
    % all parameters are within these thresholds, the variable Window is
    % set to 1 (which means that the window ist open), otherwise Window is
    % set to zero (i.e. closed).
    
    if (DistMoonPup > MinDistMoonPup) && (AltSun < MaxAltSun) && (AltPup > MinAltPup) %%&& (Absorption < MaxAbsorption*1E6)
        Window = 1;
    else
        Window = 0;
        x = 10;
    end

    
    % This block converts the variable Time into a more convenient format of
    % hours, minutes and seconds

    Hour = floor(Time);
    Min = floor((Time-Hour)*60);
    Sec = floor((Time-Hour-Min/60)*3600);
    
    
    % This writes the data into the file, if the Launch Window is open
    % (i.e. if Window == 1) and only writes every tenth value (i.e. x has to
    % be larger than 9 at the same time).
    
    % if ((Window == 1) && (x > 9)) || (Window ~= OldWindow)
    %    Data = [ActualDayNumber ActualDay ActualMonth ActualYear Hour Min Sec AltSun AzSun AltPup AzPup AltMoon AzMoon DistMoonPup Absorption Window];
    %    dlmwrite('launchwindow.txt',Data, '-append', 'delimiter', ';', 'precision', 12)
    %    x = 0;
    % end

    % This writes the data into the file, if the Launch Window status
    % changes (Window ~= OldWindow)
    
    %if (Window ~= OldWindow)
    if Window==1
       Data = [ActualDayNumber ActualDay ActualMonth ActualYear Hour Min Sec AltSun AzSun DistSunPup AltPup AzPup AltMoon AzMoon MoonPhase DeltaV DistMoonPup Absorption Window];
       dlmwrite(save_filename,Data, '-append', 'delimiter', ';', 'precision', 12)
       x = 0;
    end
    
    x = x + 1;
    OldWindow = Window;
    
    clc
    disp(['Date (DD-MM-YYYY): ', num2str(ActualDay), '-', num2str(ActualMonth), '-', num2str(ActualYear), '            Time (HH:MM:SS): ', num2str(Hour),':', num2str(Min),':', num2str(Sec)])
    disp(['JDN;ActualDay;ActualMonth;ActualYear;Hour;Min;Sec;AltSun;AzSun;DistSunPup;AltPup;AzPup;AltMoon;AzMoon;MoonPhase;DeltaV;DistMoonPup;Absorption;LaunchWindow;', '-append', 'delimiter', '']);
    disp([DistMoonPup,AltSun,Absorption,AltPup]);
end
