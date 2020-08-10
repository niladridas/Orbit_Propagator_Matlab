%--------------------------------------------------------------------------
% Version  Date         Author      Mod History
%--------------------------------------------------------------------------
% 1        25 March 2018  Vedang D.   Inital version
%--------------------------------------------------------------------------

% non linear equality constraint function - ellipse equation with given 
% center, semi major and semi minor axes

function [c, ceq] = ellipseEqn(closePt,center,a,b)
    c = [];
    ceq = (closePt(1)-center(1))^2/a^2 + (closePt(2)-center(2))^2/b^2 - 1;
end