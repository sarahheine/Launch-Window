function[nut_long,nut_obliq] = nutate(jd)
% ;+
% ; NAME:
% ;       NUTATE
% ; PURPOSE:
% ;       Return the nutation in longitude and obliquity for a given Julian date
% ;
% ; CALLING SEQUENCE:
% ;       NUTATE, jd, Nut_long, Nut_obliq
% ;
% ; INPUT:
% ;       jd - Julian ephemeris date, scalar or vector, double precision  
% ; OUTPUT:
% ;       Nut_long - the nutation in longitude, same # of elements as jd
% ;       Nut_obliq - nutation in latitude, same # of elements as jd
% ;
% ; EXAMPLE:
% ;       (1) Find the nutation in longitude and obliquity 1987 on Apr 10 at Oh.
% ;              This is example 22.a from Meeus
% ;        IDL> jdcnv,1987,4,10,0,jul
% ;        IDL> nutate, jul, nut_long, nut_obliq
% ;             ==> nut_long = -3.788    nut_obliq = 9.443
% ;            
% ;       (2) Plot the large-scale variation of the nutation in longitude 
% ;               during the 20th century
% ;
% ;       IDL> yr = 1900 + indgen(100)     ;Compute once a year        
% ;       IDL> jdcnv,yr,1,1,0,jul          ;Find Julian date of first day of year
% ;       IDL> nutate,jul, nut_long        ;Nutation in longitude
% ;       IDL> plot, yr, nut_long
% ;
% ;       This plot will reveal the dominant (18.6 year) period, but a finer
% ;       grid is needed to display the shorter periods in the nutation.
% ; METHOD:
% ;       Uses the formula in Chapter 22 of ``Astronomical Algorithms'' by Jean 
% ;       Meeus (1998, 2nd ed.) which is based on the 1980 IAU Theory of Nutation
% ;       and includes all terms larger than 0.0003".
% ;
% ; PROCEDURES CALLED:
% ;       POLY()                       (from IDL User's Library)
% ;       CIRRANGE, ISARRAY()          (from IDL Astronomy Library)
% ;
% ; REVISION HISTORY:
% ;       Written, W.Landsman (Goddard/HSTX)      June 1996       
% ;       Converted to IDL V5.0   W. Landsman   September 1997
% ;       Corrected minor typos in values of d_lng W. Landsman  December 2000
% ;       Updated typo in cdelt term              December 2000
% ;       Avoid overflow for more than 32767 input dates W. Landsman January 2005
% ;-
%  compile_opt idl2
%  On_error,2
%  
%  if N_params() LT 2 then begin
%         print,'Syntax - NUTATE, jd, nut_long, nut_obliq'
%         return
%  endif

 dtor = pi/180.0;
 %  form time in Julian centuries from 1900.0

 T = (jd - 2451545.0)/36525.0;


% Mean elongation of the Moon

   coeff1 = [297.85036,  445267.111480, -0.0019142, 1./189474d0 ];
  d = polyval(fliplr(coeff1),T)*dtor;
  d = cirrange(d,1);

% Sun's mean anomaly

   coeff2 = [357.52772, 35999.050340, -0.0001603, -1./3d5 ];
   M = polyval(fliplr(coeff2),T)*dtor;
   M = cirrange(M,1);

% Moon's mean anomaly

   coeff3 = [134.96298, 477198.867398, 0.0086972, 1.0/5.625d4 ];
   Mprime = polyval(fliplr(coeff3),T)*dtor;
   Mprime = cirrange(Mprime,1);

% Moon's argument of latitude

    coeff4 = [93.27191, 483202.017538, -0.0036825, -1.0/3.27270d5 ];
    F = polyval(fliplr(coeff4),T )*dtor ;
    F=cirrange(F,1);

% Longitude of the ascending node of the Moon's mean orbit on the ecliptic,
%  measured from the mean equinox of the date

  coeff5 = [125.04452, -1934.136261, 0.0020708, 1./4.5d5];
  omega = polyval(fliplr(coeff5),T)*dtor;
  omega=cirrange(omega,1);

 d_lng = [0,-2,0,0,0,0,-2,0,0,-2,-2,-2,0,2,0,2,0,0,-2,0,2,0,0,-2,0,-2,0,0,2,...
   -2,0,-2,0,0,2,2,0,-2,0,2,2,-2,-2,2,2,0,-2,-2,0,-2,-2,0,-1,-2,1,0,0,-1,0,0, ...
     2,0,2];

 m_lng = [0,0,0,0,1,0,1,0,0,-1,zeros(1,17),2,0,2,1,0,-1,0,0,0,1,1,-1,0, ...
  0,0,0,0,0,-1,-1,0,0,0,1,0,0,1,0,0,0,-1,1,-1,-1,0,-1];

 mp_lng = [0,0,0,0,0,1,0,0,1,0,1,0,-1,0,1,-1,-1,1,2,-2,0,2,2,1,0,0,-1,0,-1, ...
   0,0,1,0,2,-1,1,0,1,0,0,1,2,1,-2,0,1,0,0,2,2,0,1,1,0,0,1,-2,1,1,1,-1,3,0];

 f_lng = [0,2,2,0,0,0,2,2,2,2,0,2,2,0,0,2,0,2,0,2,2,2,0,2,2,2,2,0,0,2,0,0, ...
   0,-2,2,2,2,0,2,2,0,2,2,0,0,0,2,0,2,0,2,-2,0,0,0,2,2,0,0,2,2,2,2];

 om_lng = [1,2,2,2,0,0,2,1,2,2,0,1,2,0,1,2,1,1,0,1,2,2,0,2,0,0,1,0,1,2,1, ...
   1,1,0,1,2,2,0,2,1,0,2,1,1,1,0,1,1,1,1,1,0,0,0,0,0,2,0,0,2,2,2,2];

 sin_lng = [-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217, ...
    -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29, 29, 26, -22, ...
     21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7, -7, -7, ...
     6,6,6,-6,-6,5,-5,-5,-5,4,4,4,-4,-4,-4,3,-3,-3,-3,-3,-3,-3,-3];
 
 sdelt = [-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0, -0.5, 0, 0.1, ...
     0,0,0.1, 0,-0.1,zeros(1,10), -0.1, 0, 0.1, zeros(1,33)] ;


 cos_lng = [ 92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95,0,-70,-53,0, ...
    -33, 26, 32, 27, 0, -24, 16,13,0,-12,0,0,-10,0,-8,7,9,7,6,0,5,3,-3,0,3,3,...
     0,-3,-3,3,3,0,3,3,3, zeros(1,14)];

 cdelt = [8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3, zeros(1,53)];


% Sum the periodic terms 

 n = 1;%N_elements(jd);
 nut_long = 0;
 nut_obliq = 0;
 arg = d_lng*d + m_lng*M +mp_lng*Mprime + f_lng*F +om_lng*omega;
 sarg = sin(arg);
 carg = cos(arg);
 for i=1:n
        nut_long(i) =  0.0001*sum( (sdelt*T(i) + sin_lng)*sarg(1:end,i) );
        nut_obliq(i) = 0.0001*sum( (cdelt*T(i) + cos_lng)*carg(1:end,i) );
 end
 %if ~isarray(jd) then begin
        nut_long = nut_long(1);
        nut_obliq = nut_obliq(1);
 %endif

 end