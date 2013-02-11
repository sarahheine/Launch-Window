file = '/Users/Oakley/Documents/Work/MATLAB/launch_window/launchwindow.txt';
[jdn, day, month, year, hours, minute, second,altsun,azsun, distsunpup,altpup,azpup,altmoon,azmoon,moonphase,deltav,distmoonpup,absorption,window] = textread(file,'%f %u %u %u %u %u %u %f %f %f %f %f %f %f %f %f %f %f %u','delimiter',';','headerlines',15);


%mountain_hours = floor((jdn-floor(jdn)-.5)*24)-7;   %convert from UTC to Mountain Time


mountain_hours = hours - 7;
mountain_day = day;
mountain_day(find(mountain_hours < 0)) = mountain_day(find(mountain_hours < 0)) - 1;

mountain_hours(find(mountain_hours < 0)) = mountain_hours(find(mountain_hours < 0)) + 24;

DST = find((month == 11 & mountain_day < 4) | (month == 10) | (month == 3 & day > 10 | month == 4));
mountain_hours(DST) = mountain_hours(DST) + 1;



%convert_to_military = find(hours > 6);  %All times will be <6AM or >6PM
%hours(convert_to_military) = hours(convert_to_military) + 12

jdn_floor = floor(jdn);

optimal_alt = 14;
%acceptable_alt = 10; %Not needed as every value is >10 already
optimal_deltav = -1;
%acceptable_deltav = -4; %Not needed as every value is >-4 already

window_jdn_floored = floor(jdn);
window_day = zeros(max(window_jdn_floored)-min(window_jdn_floored)+1);

dlmwrite('launch_table.txt','Month,Day,Year,Acceptable Start Hour,Acceptable Start Minute,Optimal Start Hour,Optimal Start Minute,Optimal Stop Hour,Optimal Stop Minute,Acceptable Stop Hour,Acceptable Stop Minute','delimiter','');

for alldays = min(window_jdn_floored):max(window_jdn_floored)
    optimal_location = find(jdn > (alldays+4./24) & jdn < (alldays + 1 + 4./24) & altpup > optimal_alt & deltav > optimal_deltav);
    acceptable_location = find(jdn > (alldays+4./24) & jdn < (alldays + 1+4./24));
    
    %optimal_location = find(window_jdn_floored == alldays & altpup > optimal_alt & deltav == 0);
    %acceptable_location = find(window_jdn_floored == alldays);

    
    %     if (hours(max(acceptable_location)) == 12 & minute(max(acceptable_location))) == 0
%         optimal_location = find(jdn > alldays & jdn > alldays+.5 & altpup > optimal_alt & deltav == 0);
%         acceptable_location = find(jdn > alldays & jdn < alldays + .5);
%     end
       
    if size(acceptable_location,1) > 1
        start_acceptable = min(acceptable_location);
        start_optimal = min(optimal_location);
        stop_optimal = max(optimal_location);
        stop_acceptable = max(acceptable_location);
    
        if day(start_acceptable) == 3
            disp('stop')
        end
        
        if size(optimal_location,1) > 1  
%            data = [month(start_acceptable),day(start_acceptable),year(start_acceptable),hours(start_acceptable),minute(start_acceptable),hours(start_optimal),minute(start_optimal),hours(stop_optimal),minute(stop_optimal),hours(stop_acceptable),minute(stop_acceptable)];
            data = [month(start_acceptable),mountain_day(start_acceptable),year(start_acceptable),mountain_hours(start_acceptable),minute(start_acceptable),mountain_hours(start_optimal),minute(start_optimal),mountain_hours(stop_optimal),minute(stop_optimal),mountain_hours(stop_acceptable),minute(stop_acceptable)];
        else
%            data = [month(start_acceptable),day(start_acceptable),year(start_acceptable),hours(start_acceptable),minute(start_acceptable),0,0,0,0,hours(stop_acceptable),minute(stop_acceptable)];            
            data = [month(start_acceptable),mountain_day(start_acceptable),year(start_acceptable),mountain_hours(start_acceptable),minute(start_acceptable),0,0,0,0,mountain_hours(stop_acceptable),minute(stop_acceptable)];            
        end
        dlmwrite('launch_table.txt',data, '-append', 'delimiter', ',', 'precision', 12);
    end
    
end




    
    