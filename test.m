jd = 2.456210603472222e+06;
altmoon=zeros(3650,1);
azmoon=altmoon;


for k=1:3650
    [altmoon(k),azmoon(k)] = MoonPosition(jd+k/10,Lat,Long,0);
end

