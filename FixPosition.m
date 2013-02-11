function[Alt, Az] = FixPosition(Time, Lat, Long, RA, Dec)

    % The orbital elements of the Sun and the tilt of the eccliptic plane
    % are calculated in the following section:

    % Micro-X Launch Window Calculations
    % Function Fix Position
    % Version 1.0

    w = (282.9404 + 0.0000470935 * Time)/360*2*pi;           % Argument of Perihelion in rad
    M = Reduce((356.0470 + 0.9856002585 * Time)/360*2*pi);   % Mean Anomaly  in rad
    
    
    % Calculation of the Greenwich Mean Sideral Time (GMST), the Local Sideral Time
    % (LST), and the objects LHA at WSMR:

    GMST = Reduce(w + M + pi);                                      % Greenwich Mean Sideral Time in rad
    LST = GMST + (Time-floor(Time))*24*(15.0/360*2*pi) + Long;      % Local Sideral Time in rad
    LHA = LST - RA;                                                 % Local Hour Angle in rad
    
    
    % Convert these results into Altitude and Azimuth as seen from WSRM.
    % These are the two values returned by this function.
    
    Alt = asin(sin(Lat)*sin(Dec)+cos(Lat)*cos(Dec)*cos(LHA));                % Altitude at WSMR 
    Az = Reduce(pi+atan2(sin(LHA),(sin(Lat)*cos(LHA)-cos(Lat)*tan(Dec))));   % Azimuth at WSMR


