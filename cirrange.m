function [ output_angle ] = cirrange( input_angle,radians )
% ;+
% ; NAME:
% ;       CIRRANGE
% ; PURPOSE:
% ;       To force an angle into the range 0 <= ang < 360.
% ; CALLING SEQUENCE:
% ;       CIRRANGE, ang, [/RADIANS]
% ;
% ; INPUTS/OUTPUT:
% ;       ang     - The angle to modify, in degrees.  This parameter is
% ;                 changed by this procedure.  Can be a scalar or vector.
% ;                 The type of ANG is always converted to double precision
% ;                 on output.
% ;
% ; OPTIONAL INPUT KEYWORDS:
% ;       /RADIANS - If present and non-zero, the angle is specified in
% ;                 radians rather than degrees.  It is forced into the range
% ;                 0 <= ang < 2 PI.
% ; PROCEDURE:
% ;       The angle is transformed between -360 and 360 using the MOD operator.   
% ;       Negative values (if any) are then transformed between 0 and 360
% ; MODIFICATION HISTORY:
% ;       Written by Michael R. Greason, Hughes STX, 10 February 1994.
% ;       Get rid of WHILE loop, W. Landsman, Hughes STX, May 1996
% ;       Converted to IDL V5.0   W. Landsman   September 1997
% ;-



%  Determine the additive constant.

 %if keyword_set(RAD) then 
 %OUTPUT TO RADIANS
 if radians == 1
     cnst = pi * 2.0;
 end
 if radians == 0
    cnst = 360.;
 end
 
% Deal with the lower limit.

 %ang = input_angle mod cnst
 ang = rem(input_angle,cnst);
 
% Deal with negative values, if any
 
 neg = find(ang < 0.);
 ang(neg) = ang(neg) + cnst;

 output_angle=ang;
 
 end
