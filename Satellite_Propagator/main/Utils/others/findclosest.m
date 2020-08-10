function closest_index = findclosest(Mjd_UTC,vector_mjdutc)
% Input: the target mjd time and a vector of mjd time from CPF
% Output: the indices of the vector_mjdutc which is just less than Mjd_UTC  
m  = length(vector_mjdutc);
error = Mjd_UTC*ones(m,1) - vector_mjdutc;
for i = 1:m
    if error(i)<0
        closest_index = i-1;
        break
    elseif i == m
        error('the first argument is larger than the elements in the second argument')
    end
end
end