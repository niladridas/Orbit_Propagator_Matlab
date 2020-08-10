% Converting ECI elements to COE
function coe_el = eci2coe(eci_el,mu)
    % Format of ECI elememts
    r_vec = eci_el(1,1:3)';
    v_vec = eci_el(1,4:6)';
    % Algorithm 9 from book of Vallado
    % Step 1:
    h_vec = cross(r_vec,v_vec);
    h = norm(h_vec);
    % Step 2:
    n_vec = cross([0;0;1],h_vec);
    e_vec = ((dot(v_vec,v_vec)-(mu/norm(r_vec)))*r_vec - dot(r_vec,v_vec)*v_vec)/mu;
    e = norm(e_vec);
    % Step 3:
    zeta = (dot(v_vec,v_vec)/2)-(mu/norm(r_vec));
    if e ~= 1
        a = -mu/(2*zeta);
        p = a*(1-e^2);
    else
        disp('THE ORBIT IS CIRCULAR')
        p = (h^2)/mu;
        a = INF;    
    end
    % Step 4:
    i = acos(dot(h_vec,[0;0;1])/norm(h_vec));
    Ohm = acos(dot(n_vec,[1;0;0])/norm(n_vec));
    if dot(n_vec,[0;1;0])<0
        Ohm = 2*pi-Ohm;
    end
    % Step 5:
    w = acos(dot(n_vec,e_vec)/(norm(n_vec)*norm(e_vec)));
    if dot(e_vec,[0;0;1])<0
        w = 2*pi - w;
    end
    % Step 6:
    nu = acos(dot(e_vec,r_vec)/(norm(e_vec)*norm(r_vec)));
    if dot(r_vec,v_vec) <0
        nu = 2*pi - nu;
    end
%     disp('THE COE elements are in radian where applicable');
    coe_el = [a,e,i,w,Ohm,nu];
end