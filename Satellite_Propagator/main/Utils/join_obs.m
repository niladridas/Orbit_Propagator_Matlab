function observation = join_obs(pass_numbers,passes)
    i = pass_numbers(1); j = pass_numbers(2);
    observation.mjdutc = [passes(i).mjdutc;passes(j).mjdutc];
    observation.range = [passes(i).range;passes(j).range];
    observation.P = [passes(i).P;passes(j).P];
    observation.T = [passes(i).T;passes(j).T];
    observation.Rh = [passes(i).Rh;passes(j).Rh];
    