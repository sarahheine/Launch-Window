function [AltSun, AzSun, RASun, DecSun] = SunPosition(Time, Lat, Long)

    % Micro-X Launch Window Calculations
    % Function Sun Position
    % Version 1.0

    % The orbital elements of the Sun and the tilt of the eccliptic plane
    % are calculated in the following section:

    w = (282.9404 + 0.0000470935 * Time)/360*2*pi;           % Argument of Perihelion in rad
    M = Reduce((356.0470 + 0.9856002585 * Time)/360*2*pi);   % Mean Anomaly in rad
    e = 0.016709 - 0.000000001151 * Time;                    % Eccentricity
    E = M + e*sin(M).*(1+e*cos(M));                          % Eccentric Anomaly
    ecl = (23.4393 - 0.0000003563 * Time)/360*2*pi;          % Tilt of eccliptic plane in rad


    % Calculation of the Greenwich Mean Sideral Time (GMST) and the Local Sideral Time
    % (LST) at WSMR:

    GMST = Reduce(w + M + pi);                                    % Calculated in rad, reduced to below 2pi
    LST = GMST + (Time-floor(Time))*24*(15.0/360*2*pi) + Long;    % This is the sideral time at WSMR in rad


    % The Sun's distance r and true anomaly v are calcuated from its position in its orbit:

    xv = cos(E) - e;
    yv = sqrt(1.0 - e^2) * sin(E);
    v = atan2(yv, xv);                     % This is the sun's true anomaly
    r = sqrt(xv^2 + yv^2);                 % The suns distance to the Earth (in AU)

    
    % Calculate the sun's true longitude and its position in the ecliptic rectangular
    % geocentric coordinate system. zs is zero, as the sun is in the ecliptic.
    
    lonsun = v + w;                           % This is the sun's true longitude
    xs = r * cos(lonsun);                     % x coordinate in ecliptic rectangular geocentric system
    ys = r * sin(lonsun);                     % y coordinate in ecliptic rectangular geocentric system

    
    % Convert that into equatorial, rectangular, geocentric coordinates:

    xe = xs;
    ye = ys * cos(ecl);
    ze = ys * sin(ecl);

    
    % Now, the Sun's Right Ascension, Declination and Local Hour Angle

    RASun  = atan2(ye,xe);
    DecSun = atan2(ze,sqrt(xe^2+ye^2));
    LHASun = LST - RASun;  
    
    
    % Convert these results into Altitude and Azimuth as seen from WSRM.
    % These are the two values returned by this function.
    
    AltSun = asin(sin(Lat)*sin(DecSun)+cos(Lat)*cos(DecSun)*cos(LHASun));                % Altitude at WSMR 
    AzSun = Reduce(pi+atan2(sin(LHASun),(sin(Lat)*cos(LHASun)-cos(Lat)*tan(DecSun))));   % Azimuth at WSMR

