% JB2008 drag model requires as an input the day of the year
% Not the day of month
% This is an util function that gives day of year

function dayofyr = dayofyear(year, month, day)
    % detect leap year
    daysuptomnth = 30*(month-1)+ceil((month-1)/2);
    if month > 2
        if mod(year,4)==0
            daysuptomnth = daysuptomnth-1;
        else
            daysuptomnth = daysuptomnth-2;
        end
    end
    dayofyr = daysuptomnth + day;
end