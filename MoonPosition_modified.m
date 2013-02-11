function[AltMoon, AzMoon, RAMoon, DecMoon, MoonDistancekm] = MoonPosition(Time, Lat, Long, Alt)

    % Function Moon Position
    % Version 1.0
    
    % The orbital elements of the Sun and tilt of the ecliptic plane
    % are calculated in the following section:

    wSun = (282.9404 + 0.0000470935 * Time)/180*pi;             % Argument of Perihelion in rad
    MSun = Reduce((356.0470 + 0.9856002585 * Time)/180*pi);     % Mean Anomaly in rad
    ecl = (23.4393 - 0.0000003563 * Time)/360*2*pi;             % Tilt of eccliptic plane in rad


    % Calculation of the Greenwich Mean Sideral Time (GMST) and the Local Sideral Time
    % (LST) at WSMR:

    GMST = Reduce(wSun + MSun + pi);                                % Calculated in rad, reduced to below 2pi
    LST = GMST + (Time-floor(Time))*24*(15.0/360*2*pi) + Long;      % This is the sideral time at WSMR in rad


    % These the orbital elements of the Moon:

    NMoon = (125.1228 - 0.0529538083 * Time)/360*2*pi;                          % Longitude of ascending node in rad
    iMoon = 5.1454/360*2*pi;                                                    % Orbit inclination in rad
    wMoon = (318.0634 + 0.1643573223 * Time)/360*2*pi;                          % Argument of perihelion in rad
    aMoon = 60.2666;                                                            % Semi major axis in earth radii
    eMoon = 0.054900;                                                           % Eccentricity
    MMoon = Reduce((115.3654 + 13.0649929509 * Time)/360*2*pi);                 % Mean Anomaly in rad
    EMoon = MMoon + eMoon .* sin(MMoon) .* (1.0 + eMoon .* cos(MMoon));
    
    
    % Moons distance r (in Earth radii) and true anomaly v are calculated
    % from the Moon's position in its orbit:

    xvMoon = aMoon * (cos(EMoon) - eMoon);                                      
    yvMoon = aMoon * (sqrt(1.0 - eMoon^2) .* sin(EMoon));
    vMoon = atan2(yvMoon,xvMoon);                                      % True anomaly in rad
    rMoon = sqrt(xvMoon.*xvMoon + yvMoon.*yvMoon);                     % Distance Moon - Earth
    
    % This is now converted into ecliptic, rectangular, geocentric coordinates:
    
    xgMoon = rMoon .* (cos(NMoon) .* cos(vMoon+wMoon) - sin(NMoon) .* sin(vMoon+wMoon) .* cos(iMoon));
    ygMoon = rMoon .* (sin(NMoon) .* cos(vMoon+wMoon) + cos(NMoon) .* sin(vMoon+wMoon) .* cos(iMoon));
    zgMoon = rMoon .* (sin(vMoon+wMoon) .* sin(iMoon));
    
    
    % Lonecl and latecl are the ecliptic longitude and latitude. They need
    % to be calculated here since the perturbation is given as a change in
    % lonec and lated
     
    lonecl = atan2(ygMoon, xgMoon);                                 % Ecliptic longitude in rad
    latecl = atan2(zgMoon, sqrt(xgMoon*xgMoon+ygMoon*ygMoon));      % Ecliptic latitudue in rad
    
    
    % Perturbation needs to be taken into account for an accuracy better
    % then 2 degrees. With the terms below, 2 arcmins should be feasible:

    LSun = MSun + wSun;                     % Mean latitude for the Sun in rad
    LMoon = MMoon + wMoon + NMoon;          % Mean latitude for the Moon in rad
    D = LMoon - LSun;                       % Mean elongation of the Moon in rad
    F = LMoon - NMoon;                      % Argument of latitude for the Moon in rad

    % Changes to ecliptic longitude:
    lonecl = lonecl - (1.274 * sin(MMoon - 2*D))/180*pi;             
    lonecl = lonecl + (0.658 * sin(2*D))/180*pi;             
    lonecl = lonecl - (0.186 * sin(MSun))/180*pi; 
    lonecl = lonecl - (0.059 * sin(2*MMoon - 2*D))/180*pi;
    lonecl = lonecl - (0.057 * sin(MMoon - 2*D + MSun))/180*pi;
    lonecl = lonecl + (0.053 * sin(MMoon + 2*D))/180*pi;
    lonecl = lonecl + (0.046 * sin(2*D - MSun))/180*pi;
    lonecl = lonecl + (0.041 * sin(MMoon - MSun))/180*pi;
    lonecl = lonecl - (0.035 * sin(D))/180*pi;              
    lonecl = lonecl - (0.031 * sin(MMoon + MSun))/180*pi;
    lonecl = lonecl - (0.015 * sin(2*F - 2*D))/180*pi;
    lonecl = lonecl + (0.011 * sin(MMoon - 4*D))/180*pi;

    % Changes to ecliptic latitude:
    latecl = latecl - (0.173 * sin(F - 2*D))/180*pi;
    latecl = latecl - (0.055 * sin(MMoon  - F - 2*D))/180*pi;
    latecl = latecl - (0.046 * sin(MMoon  + F - 2*D))/180*pi;
    latecl = latecl + (0.033 * sin(F + 2*D))/180*pi;
    latecl = latecl + (0.017 * sin(2*MMoon + F))/180*pi;
 
    % Changes to distance from earth:
    rMoon = rMoon - 0.58 * cos(MMoon - 2*D) - 0.46 * cos(2*D);

    
    % And now, the position is again calculated in ecliptic, rectangular, geocentric
    % coordinates:
    
    xgMoon = rMoon * cos(lonecl) * cos(latecl);
    ygMoon = rMoon * sin(lonecl) * cos(latecl);
    zgMoon = rMoon * sin(latecl);

    
    % This is converted into equatorial, rectangular, geocentric coordinates:

    xeMoon = xgMoon;
    yeMoon = ygMoon .* cos(ecl) - zgMoon .* sin(ecl);
    zeMoon = ygMoon .* sin(ecl) + zgMoon .* cos(ecl);

    
    % Compute the Moon's Right Ascension (RAMoon) and Declination (DecMoon):

    RAMoon  = atan2(yeMoon, xeMoon);
    DecMoon = atan2(zeMoon, sqrt(xeMoon.*xeMoon+yeMoon.*yeMoon));
    EarthRadius = 6378.1; %approximate value in km
    MoonDistancekm = rMoon * 6378.1; %Added by Phil 1/10/2012

    
    % Since the moon is so close to earth, the observation altitude has to be taken into account.
    % RA and Dec have to be converted to their values as seen from rocket altitude:

    mpar = asin(1/rMoon);               % Moons parallax in rad
    gclat = Lat;                        % This is the geocentric latitude. We assume that the rocket launches straigt up.
    rho = 1.0+(Alt/6371);               % The oberservation altitude is calculated (distance from Earth's center in earth radii)
    HA = LST - RAMoon;                  % The moons hour angle in rad
    g = atan(tan(gclat)/cos(HA));       % Auxiliary angel in rad
    
    topRAMoon   = RAMoon - mpar * rho * cos(gclat) * sin(HA) / cos(DecMoon);            % RA of the Moon in rad as seen from rocket altitude abouve WSMR
    topDecMoon = DecMoon - mpar * rho * sin(gclat) * sin(g - DecMoon) / sin(g);         % Dec of the Moon in rad as seen from rocket altitude abouve WSMR
    LHAMoon = LST - topRAMoon;                                                          % LHA of the Moon in rad as seen from rocket altitude abouve WSMR

   
    % Convert these results into Altitude and Azimuth as seen at rocket altitude above WSMR.
    % These are the two values returned by this function.
    
    AltMoon = asin(sin(Lat)*sin(topDecMoon)+cos(Lat)*cos(topDecMoon)*cos(LHAMoon));             % Altitude at WSMR and rocket altitude in rad
    AzMoon = Reduce(pi+atan2(sin(LHAMoon),(sin(Lat)*cos(LHAMoon)-cos(Lat)*tan(topDecMoon))));   % Azimuth at WSMR
    

