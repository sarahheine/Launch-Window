file = '/Users/Oakley/Documents/Work/MATLAB/launch_window/casA_launchwindow_2015.txt';
window = dlmread(file,';',19,0);

julian      = window(1:end,1);
day         = window(1:end,2);
month       = window(1:end,3);
year        = window(1:end,4);
hour        = window(1:end,5);
minute      = window(1:end,6);
second      = window(1:end,7);
altsun      = window(1:end,8);
azsun       = window(1:end,9);
distsunpup  = window(1:end,10);
altpup      = window(1:end,11);
azpup       = window(1:end,12);
altmoon     = window(1:end,13);
azmoon      = window(1:end,14);
moonphase   = window(1:end,15);
deltav      = window(1:end,16);
distmoonpup = window(1:end,17);
absorption  = window(1:end,18);
launch      = window(1:end,19);

firstof2015 = 2457023.500000;
firstof2013 = 2456293.500000;

%firstof2015 = firstof2013;

[year month day hour min sec] = jd2date(julian);

testdate = [num2str(year) num2str(month) num2str(year) num2str(year) num2str(year)]; 

brightness_ratio = 10.^(deltav/-2.5);

plot(julian-firstof2015,altpup,'.');
xlim([0,365]);

hold on

%[AX,H1,H2] = plotyy(julian-firstof2012,altpup,julian-firstof2012,brightness_ratio);
%set(AX(2),'Ylim',[-20,20]); 
%set(H1,'Marker','.');
%set(H2,'LineStyle','.');
%set(AX(1),'Ylim',[-100,150]) ;


plot(julian-firstof2015,altsun,'r.');
plot(julian-firstof2015,distmoonpup,'m.');
plot(julian-firstof2015,distsunpup,'g.');
plot(julian-firstof2015,brightness_ratio,'c.');

scatter(julian - firstof2015,altmoon,30,[1-moonphase 1-moonphase 1-moonphase]);


%month_numbers = [306, 336, 1+366, 32+366, 60+366, 91+366, 121+366];
month_numbers = [1,32,60,91,121,152,182,213,244,274,305,335];
%month_names = {'November','December','January','February','March','April','May'};
month_names = {'January','February','March','April','May','June','July','August','September','October','November','December'};
y_boundaries = ylim;
month_text_y = month_numbers * 0. - 90;
month_text_x = month_numbers + 1;
for month = 1:size(month_numbers,2)
    plot([month_numbers(month),month_numbers(month)],[y_boundaries(1),y_boundaries(2)],'k--');
end

text(month_text_x,month_text_y,month_names);




xlabel('Days since January 1, 2013');
ylabel('Altitude, Angular Separation or Brightness Ratio [Degrees]');
legend('Altitude of Cas-A','Altitude of the Sun','Angular Separation of the Moon/Cas-A','Angular Separation of the Sun/Cas-A','(Moon + Sky Flux) / Dark Sky Flux','Altitude of the Moon');

ah = gca;
tah = addtxaxis(ah,'x-2457023.5',month_numbers+firstof2015,'Julian Date [January 1, 2013 - December 31, 2013]');

