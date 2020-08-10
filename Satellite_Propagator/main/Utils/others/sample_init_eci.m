function pos_vel = sample_init_eci(full_state,samples, radius_limit,h_maglimit, h_anglimit)
angles = 2*pi*rand(2,samples)-pi*ones(2,samples);
    radius = 2*radius_limit*rand(1,samples)-radius_limit*ones(1,samples);
    x = radius.*cos(angles(2,:)).*cos(angles(1,:));
    y = radius.*sin(angles(2,:)).*sin(angles(1,:));
    z = radius.*sin(angles(2,:));
    xyz_prior = [x;y;z];
    xyz_samples = full_state(1:3,1).*ones(3,samples)+xyz_prior;
    pos_vel = zeros(6,samples);
    h_0 = cross(full_state(1:3,1),full_state(4:6,1));
    h_0_norm = h_0/norm(h_0);
    i = 1;
    count = 0;  
    while 1 
      del_v = 1*rand(3,1)-0.5*ones(3,1); 
      vel = full_state(4:6,1)+del_v;
      if i > samples
          break
      end
      new_h = cross(xyz_samples(:,i),vel);
      new_h_norm = new_h/norm(new_h);
      if abs(acos(dot(h_0_norm,new_h_norm)))<= h_anglimit && abs(norm(new_h)-norm(h_0))<= h_maglimit
          pos_vel(:,i) = [xyz_samples(:,i);vel];
          i = i+1;
      else
          count = count+1;
          if count >= 300
              error('too many iterations');
          end
          continue
      end
      count = 0;
    end


%     pos_vel = [xyz_samples;repmat(full_state(4:6,1),1,samples)];
end