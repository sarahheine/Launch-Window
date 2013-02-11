function [DeltaV,D,ObjMoonDist,K]=moon_sky_brightness(Date,ObjCoo,GeodPos,C_Ext,Vsky);
%-----------------------------------------------------------------------
% moon_sky_brightness function     Calculate sky brightness due to the
%                                moon for a given date and sky position,
%                                by taking into acount the distance from
%                                the moon and zenith distance (V mag).
% Input  : - Date [day, month, year, frac_day] or JD.
%            if one element is given than assumed to be JD.
%          - Object apparent equatorial
%            coordinates [RA, Dec] in radians.
%          - Observer geodetic position [East_Long, Lat, Height],
%            radians and meters above ref ellips.
%            default is wise observatory position.
%          - Extinction coef. in V. (default is 0.3mag/airmass).
%          - Sky brightness in V. (default is 21.7 mag/sq. arcsec.).
% Output : - The change in the V-band sky brightness caused by moonlight.
%          - Moon elongation, radians.
%          - Object-Moon distance, radians.
%          - Moon illuminated fraction.
% Reference : Krisciunas, K. and Schaefer, B. 1991 PASP 103, 1033.
% Tested : Matlab 5.3
%     By : Eran O. Ofek            August 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------
RAD = 57.29577951308232;


if (nargin==2),
   GeodPos = [34.7630./RAD, 30.5960./RAD, 0];
   C_Ext      = 0.3;
   Vsky       = 21.7;
elseif (nargin==3),
   C_Ext      = 0.3;
   Vsky       = 21.7;
elseif (nargin==4),
   Vsky       = 21.7;
elseif (nargin==5),
   % do nothing
else
   error('Illigal number of input arguments');
end


% Julian Day
if (length(Date(1,:))>1),
   JD = julday([Date(1),Date(2),Date(3),Date(4)]);
else
   JD = Date;
end

% Geodetic to Geocentric
[GeocPos]=geod2geoc(GeodPos(2),GeodPos(3),'WGS84');

% Moon & Sun Position
[MoonRA,MoonDec,MoonHP] = mooncool(JD,[GeodPos(1),GeocPos]);
[SunRA,SunDec,SunR]     = suncoo(JD,'a');
MoonR = asin(MoonHP).*6378.137./149597870.0; %AU

% Moon Elongation
D = distsp(MoonDec,SunDec,MoonRA,SunRA);
Alpha = RAD.*(pi - D);


% Selenographic elongation of the Earth from the Sun
I = atan2((SunR.*sin(D)),(MoonR - SunR.*cos(D)));

% Illuminated fraction
K = 0.5.*(1 + cos(I));




% convert obj coo. to horiz coo.
ObjHor = horiz_coo(ObjCoo,JD,GeodPos,'h');
%ObjHor = ObjCoo;
Z = pi./2 - ObjHor(2);

% Moon horizontal coo.
MoonHor = horiz_coo([MoonRA,MoonDec],JD,GeodPos,'h');
Z_Moon = pi./2 - MoonHor(:,2);


% Object-Moon distance
%ObjMoonDist = distsp(ObjCoo(2),MoonDec,ObjCoo(1),MoonRA);
ObjMoonDist = distsp(ObjHor(2),MoonHor(2),ObjHor(1),MoonHor(1));



% moon ilumination
I_star = 10.^(-0.4.*(3.84 + 0.026.*abs(Alpha) + 4e-9.*Alpha.^4));
F_Rho  = (10.^5.36).*(1.06 + cos(ObjMoonDist).^2) + 10.^(6.15 - RAD.*ObjMoonDist./40);

Xz     = (1 - 0.96.*sin(Z).^2).^(-0.5);
XzMoon = (1 - 0.96.*sin(Z).^2).^(-0.5);

% moon sky brigtness in nanoLamberts
%Bmoon  = F_Rho.*I_star.*(1 - 10.^(-0.4.*C_Ext.*Xz)).*10.^(-0.4.*C_Ext.*Xz);
Bmoon  = F_Rho.*I_star.*(1 - 10.^(-0.4.*C_Ext.*Xz)).*10.^(-0.4.*C_Ext.*XzMoon);
IbelowHor = find(MoonHor(:,2)<0);
Bmoon(IbelowHor) = 0;

% convert nanLamberts to mag/sq. arcsec.
%Bsky   = -(log(Bmoon./34.08) - 20.7233)./0.92104


% convert sky brightness to nanoLamberts
VskyNL = 34.08.*exp(20.7233 - 0.92104.*Vsky);

DeltaV = -2.5.*log10((Bmoon + VskyNL)./VskyNL);









% distance function
function Dist=distsp(D1,D2,R1,R2)
Dist = acos(sin(D1).*sin(D2) + cos(D1).*cos(D2).*cos(R1-R2));

function [RA,Dec,HP]=mooncool(Date,EarthPos,Algo);
%------------------------------------------------------------------------
% mooncool function  Calculate low-accuracy Moon Topocentric Equatorial
%                  coordinates (Equinox of date).
% input  : - matrix od dates, [D M Y frac_day] per line,
%            or JD per line. In TT time scale.
%          - [East_Long, North_Lat] of observer in radians.
%            If NaN then calculate geocentric position.
%          - Algorithm:
%            'l' : very low accuracy (default).
%                  0.3 deg in pos. (apparent coordinates).
%                  0.003 deg. in horizontal parallax.
%            'b' : low accuracy ~1' in position.
% output : - vector of RA, in radians.
%          - vector of Dec. in radians.
%          - Vector of horizontal parallax.
%            r = 1/sin(HP)  SD = 0.2725.*HP
%    By  Eran O. Ofek           August 1999
%------------------------------------------------------------------------
RAD = 180./pi;

%FunTPI = inline('(X./(2.*pi) - floor(X./(2.*pi))).*2.*pi','X');


if (nargin==2),
   Algo = 'l';
elseif (nargin==3),
   % do nothing
else
   error('Illigal number of input arguments');
end


SizeDate = size(Date);
N        = SizeDate(1);
ColN     = SizeDate(2);

if (ColN==4),
   JD = julday(Date).';
elseif (ColN==1),
   JD = Date;
else
   error('Illigal number of columns in date matrix');
end

T   = (JD - 2451545.0)./36525.0;

switch Algo
 case 'l'

    n  = JD - 2451545.0;
    if (max(n)>50.*365 | min(n)<-50.*365),
       error('This formulae give good results only between 1950-2050');
    end

    SA1 = sin((134.9 + 477198.85.*T)./RAD);
    SA2 = sin((259.2 - 413335.38.*T)./RAD);
    SA3 = sin((235.7 + 890534.23.*T)./RAD);
    SA4 = sin((269.9 + 954397.70.*T)./RAD);
    SA5 = sin((357.5 +  35999.05.*T)./RAD);
    SA6 = sin((186.6 + 966404.05.*T)./RAD);
    Lam = (218.32 + 481267.883.*T + 6.29.*SA1 - 1.27.*SA2 + 0.66.*SA3 + 0.21.*SA4 - 0.19.*SA5 - 0.11.*SA6)./RAD;
    
    BA1 = sin((93.3  + 483202.03.*T)./RAD);
    BA2 = sin((228.2 + 960400.87.*T)./RAD);
    BA3 = sin((318.3 +   6003.18.*T)./RAD);
    BA4 = sin((217.6 - 407332.20.*T)./RAD);
    Bet = (5.13.*BA1 + 0.28.*BA2 - 0.28.*BA3 - 0.17.*BA4)./RAD;
    
    CA1 = cos((134.9 + 477198.85.*T)./RAD);
    CA2 = cos((259.2 - 413335.38.*T)./RAD);
    CA3 = cos((235.7 + 890534.23.*T)./RAD);
    CA4 = cos((269.9 + 954397.70.*T)./RAD);
    HP  = (0.9508 + 0.0518.*CA1 + 0.0095.*CA2 + 0.0078.*CA3 + 0.0028.*CA4)./RAD;    
    r = 1./sin(HP);
    
    l = cos(Bet).*cos(Lam);
    m = 0.9175.*cos(Bet).*sin(Lam) - 0.3978.*sin(Bet);
    n = 0.3978.*cos(Bet).*sin(Lam) + 0.9175.*sin(Bet);
    
    x = r.*l;
    y = r.*m;
    z = r.*n;

 case 'b'

    [Lon,Lat,Rad,HP] = moonecool(JD);

    R = 1./sin(HP);
    
    Obl = obliquity(JD);
    Nt = length(JD);
    x  = zeros(Nt,1);
    y  = zeros(Nt,1);
    z  = zeros(Nt,1);
    for I=1:1:Nt,

       L = cos(Lat(I)).*cos(Lon(I));
       M = cos(Obl(I)).*cos(Lat(I)).*sin(Lon(I)) - sin(Obl(I)).*sin(Lat(I));
       N = sin(Obl(I)).*cos(Lat(I)).*sin(Lon(I)) + cos(Obl(I)).*sin(Lat(I));
    
       x(I) = R(I).*L;
       y(I) = R(I).*M;
       z(I) = R(I).*N;   
    end

 otherwise
    error('Unknown algorithm');
end


if (isnan(EarthPos)==1),
   % geocentric coordinates
   xt = x;
   yt = y;
   zt = z;
else
   LST = lst(JD,EarthPos(1),'m').*2.*pi;
   Lat = EarthPos(2);

   xt = x - cos(Lat).*cos(LST);
   yt = y - cos(Lat).*sin(LST);
   zt = z - sin(Lat);
end

RA  = atan2(yt,xt);
Dec = asin(zt./sqrt(xt.^2 + yt.^2 + zt.^2));

function [RA,Dec,R,SL,EquationTime]=suncoo(JD,EquinoxType);
%--------------------------------------------------------------------
% suncoo function                                              ephem
% Description: Calculate the Sun equatorial coordinates using low
%              accuracy formale. Accuracy : 0.01 deg. in long.
% Input  : - Vector of JDs.
%          - Equinox for output coordinates:
%            'g' - True Geometric referred to the mean equinox of date.
%            'a' - Apparent, referred to the true equinox of date. (default).
%            'j' - J2000, referred to the J2000.0 equinox.
% Output : - vector of RA, in radians.
%          - vector of Dec. in radians.
%          - Vector of radius vectors, in AU.
%          - Solar longitude in the same ref. frame as RA/Dec. (radians).
%          - Equation of Time [Minuts of time]
% See Also : suncoo1; mooncoo
% Tested : matlab 5.3
%     By : Eran O. Ofek           September 1999
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [RA,Dec,R,SL,EquationTime]=suncoo(2451545+[0:1:10]','j');
%--------------------------------------------------------------------
if (nargin==1),
   EquinoxType = 'a';
elseif (nargin==2),
   % no default
else
   error('Illigal number of input arguments');
end
RAD = 180./pi;

T   = (JD - 2451545.0)./36525;
L0  = (280.46645 + 36000.76983.*T + 0.0003032.*T.*T)./RAD;
M   = (357.52910 + 35999.05030.*T - 0.0001559.*T.*T - 0.00000048.*T.*T.*T)./RAD;
e   = 0.016708617 - 0.000042037.*T - 0.0000001236.*T.*T;

C   = (1.914600 - 0.004817.*T - 0.000014.*T.*T).*sin(M) + (0.019993 - 0.000101.*T).*sin(2.*M) + 0.000290.*sin(3.*M);
C   = C./RAD;
% Sun longitude
SL  = L0 + C;

% the sun true Anomaly:
Ni  = M + C;

% the sun radius vector
R   = 1.000001018.*(1-e.^2)./(1+e.*cos(Ni));

if (EquinoxType=='a'),
   Om = (125.04 - 1934.136.*T)./RAD;
   SL = SL - (0.00569 - 0.00478.*sin(Om))./RAD;
elseif (EquinoxType=='j'),
   SL = SL - (0.01397.*T.*100)./RAD;
elseif (EquinoxType=='g'),
   % Allready geometric longitude
else
   error('Illigal equinox type');
end


Obl = obliquity(JD);

SL     = (SL./(2.*pi) - floor(SL./(2.*pi))).*2.*pi;
RA     = atan2(cos(Obl).*sin(SL),cos(SL));
Dec    = asin(sin(Obl).*sin(SL));


EquationTime = L0 - 0.0057183./RAD - RA;
EquationTime = 1440.*(EquationTime./(2.*pi) - floor(EquationTime./(2.*pi)));
if (EquationTime>720),
   EquationTime = EquationTime - 1440;
end

function LST=lst(JD,EastLong,STType);
%--------------------------------------------------------------------
% lst function         Local Sidereal Time, (mean or apparent),
%                    for vector of JD's and a given East Longitude.
% input  : - Vector of JD, in UT1 time scale.
%          - East Longitude in radians.
%          - Sidereal Time Type,
%            'm' - Mean (default).
%            'a' - apparent.
% output : - vector of LST in fraction of day.
%    By  Eran O. Ofek           August 1999
%--------------------------------------------------------------------
RAD = 57.29577951308232;

if (nargin==2),
   STType = 'm';
elseif (nargin==3),
   % do nothing
else
   error('Illigal number of input arguments');
end


% convert JD to integer day + fraction of day
TJD = floor(JD - 0.5) + 0.5;
DayFrac = JD - TJD;

T = (TJD - 2451545.0)./36525.0;

GMST0UT = 24110.54841 + 8640184.812866.*T + 0.093104.*T.*T - 6.2e-6.*T.*T.*T;

% convert to fraction of day in range [0 1)
GMST0UT = GMST0UT./86400.0;

GMST0UT = GMST0UT - floor(GMST0UT);
LST = GMST0UT + 1.0027379093.*DayFrac + EastLong./(2.*pi);
LST = LST - floor(LST);


switch STType
 case {'m'}
    % do nothing
 case {'a'}
    % calculate nutation
    NutMat = nutation(JD);
    Obl    = obliquity(JD);
    EquationOfEquinox = (RAD.*3600).*NutMat(:,1).*cos(Obl)./15;
    LST = LST + EquationOfEquinox./86400;    
 otherwise
    error('Unknown sidereal time type');
end

function Obl=obliquity(JulianDay,Type);
%---------------------------------------------------------
% obliquity function        calculating obliquity of
%                         ecliptic (with respect to the
%                         mean equator od date)
%                         for a given julian day.
% Input  : - vector of Julian Days.
%          - Caqlculation type:
%            'L' - IAU 1976, good from 1000-3000 AD,
%                  default.
%            'H' - Laskar expression, more accurate.
% Output : - obliquity of ecliptic of date in radians.
% Tested : Matlab 5.3
%     By : Eran O. Ofek         August 1999
%                  Last Update: October 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------
RADIAN = 180./pi;

if (nargin==1),
   Type = 'L';
elseif (nargin==2),
   % do nothing
else
   error('Illigal number of input arguments');
end

switch Type
 case 'L'
    T   = (JulianDay - 2451545.0)./36525.0;
    Obl = 23.439291 - 0.0130042.*T - 0.00000016.*T.*T + 0.000000504.*T.*T.*T;
    Obl = Obl./RADIAN;
 case 'H'
    T   = (JulianDay - 2451545.0)./36525.0;
    U   = T./100;
    Obl = 23.44484666666667 ...
       -   (4680.93.*U ...
          - 1.55.*U.^2 ...
          + 1999.25.*U.^3 ...
            - 51.38.*U.^4 ...
           - 249.67.*U.^5 ...
            - 39.05.*U.^6 ...
             + 7.12.*U.^7 ...
            + 27.87.*U.^8 ...
             + 5.79.*U.^9 ...
             + 2.45.*U.^10)./3600;
    Obl = Obl./RADIAN;
 otherwise
    error('Unknown calculation type in obliquity.m');  
end
  

function OutCoo=horiz_coo(InCoo,JD,TopoPos,Direction);
%--------------------------------------------------------------------
% horiz_coo function        Horizontal coordinates conversion
%                       Converting from equatorial coordinates
%                       to Horizontal coordinates and visa
%                       versa.
% input  : - Two columns matrix of coordinates,
%            (Long & Lat) | (Az & Alt) in radians
%          - vector of JDs + UT fraction of day,
%            if scalar value is given, then it duplicated
%            for all coordinates.
%          - Geodetic Coordinates, east long & north lat in radians
%            if scalar value is given then it is duplicate for
%            all coordinates.
%          - Direction,
%            'h' - from equatorial to horizontal (default).
%            'e' - from horizontal to equatorial.
% output : - two column matrix of output coordinates.
% Tested : matlab 5.3
%     By : Eran O. Ofek           August 1999
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------
if (nargin==3),
   Direction = 'h';
elseif (nargin==4),
   % no default
else
   error('Illigal number of input arguments');
end

N = length(InCoo(:,1));
if (length(JD)==1),
   JD = JD.*ones(N,1);
end

if (length(TopoPos(:,1))==1),
   TopoPos = ones(N,1)*TopoPos;
end


% Don't convert Geodetic latitude to Geocentric.
%GeodLat = TopoPos(:,2);
%GeocLatTemp = geod2geoc([TopoPos(:,1),GeodLat,zeros(N,1)],'WGS84');
%GeocLat     = GeocLatTemp(:,2);
%TopoPos(:,2) = GeocLat;


% calculating Local Mean Sidereal Time
LST=lst(JD,TopoPos(:,1),'m');


if (Direction=='h'),
   % convert equatorial to horizontal

   % calculate the Hour Angle
   HA  = LST.*2.*pi - InCoo(:,1);
   Dec = InCoo(:,2);
   Lat = TopoPos(:,2);

   SinAlt = sin(Dec).*sin(Lat) + cos(Dec).*cos(HA).*cos(Lat);
   CosAlt = sqrt(1-SinAlt.*SinAlt);

   SinAz  = (-cos(Dec).*sin(HA))./CosAlt;
   CosAz  = (sin(Dec).*cos(Lat) - cos(Dec).*cos(HA).*sin(Lat))./CosAlt;

   Az     = atan2(SinAz, CosAz);
   Alt    = asin(SinAlt);

   I      = find(Az<0);
   Az(I)  = 2.*pi+Az(I);

   OutCoo = [Az, Alt];
elseif (Direction=='e'),
   Az     = InCoo(:,1);
   Alt    = InCoo(:,2);
   Lat = TopoPos(:,2);

   SinDec = sin(Alt).*sin(Lat) + cos(Alt).*cos(Az).*cos(Lat);
   CosDec = sqrt(1 - SinDec.*SinDec);

   SinHA  = (-cos(Alt).*sin(Az))./CosDec;
   CosHA  = (sin(Alt).*cos(Lat) - cos(Alt).*cos(Az).*sin(Lat))./CosDec;
   HA     = atan2(SinHA, CosHA);

   RA     = LST.*2.*pi - HA;
   Dec    = asin(SinDec);

   % converting to range [0,1)
   RA     = 2.*pi.*(RA./(2.*pi) - floor(RA./(2.*pi)));

   OutCoo = [RA, Dec];
else
   error('Illigal type of conversion');
end









