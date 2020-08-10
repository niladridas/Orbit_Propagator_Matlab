% Classical Orbital Elements to Modified Equanoctial Elements
function mee_el = coe2mee(coe_el)
    % coe_el =  [a,e,i,w,Ohm,nu];
    a = coe_el(1,1);
    e = coe_el(1,2);
    i = coe_el(1,3);
    w = coe_el(1,4);
    Ohm = coe_el(1,5);
    nu = coe_el(1,6);
    
    if e == 1
        error('The orbit is circular, the code will not work')
    end
    p = a*(1-e^2);
    f = e*cos(w+Ohm);
    g = e*sin(w+Ohm);
    h = tan(i/2)*cos(Ohm);
    k = tan(i/2)*sin(Ohm);
    L = Ohm + w + nu;
    
    mee_el = [p,f,g,h,k,L];
end