% Note: This code is taken from MathWorks File Exchange.
% Author: Meysam Mahooti
%--------------------------------------------------------------------------
%
%  Inputs:
%    vec1         vector 1
%    vec2         vector 2
%
%  Output:
%    theta        angle between the two vectors [-pi,pi]
%
%--------------------------------------------------------------------------
function theta = angl(vec1, vec2)

small     = 0.00000001;
undefined = 999999.1;

magv1 = norm(vec1);
magv2 = norm(vec2);

if (magv1*magv2 > small^2)
    temp= dot(vec1,vec2)/(magv1*magv2);
    if (abs(temp) > 1)
        temp= sign(temp);
    end
    theta= acos(temp);
else
    theta= undefined;
end