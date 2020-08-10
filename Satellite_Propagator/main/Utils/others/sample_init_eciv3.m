function pos_vel = sample_init_eciv3(full_state,samples, radius_limit,station_eci)
% We assume main error along track in this work
% Cross track error is small
% First sample a torus around the base orbit for position
% values

    % Calculate the vector connecting station with the satellite 
    vec_c = full_state(1:3,1) - station_eci;
    unitvec_c = vec_c/norm(vec_c);
    % Uniformly sample from 0 to radius_limit
    sample_delvec = radius_limit*rand(1, samples);
    position_vec = repmat(full_state(1:3,1),1,samples)+repmat(unitvec_c,1,samples).*repmat(sample_delvec,3,1);
    % no error in velocity
    pos_vel = [position_vec;repmat(full_state(4:6,1),1,samples)]; 
end