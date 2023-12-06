function [flag, u, v, t] = fastRayTriangleIntersection(o, d, p0, p1, p2) 
% Ray/triangle intersection using the algorithm proposed by Möller and Trumbore (1997). 
% 
% Input: 
% IMPORTANT: ALL INPUTS NEED TO BE 3 ELEMENT ROW VECTORS 
% o : origin. 
% d : direction. 
% p0, p1, p2: vertices of the triangle. 
% Output: 
% flag: (0) Reject, (1) Intersect. 
% u,v: barycentric coordinates. 
% t: distance from the ray origin. 
% Author: 
% Originally written by Jesus Mena, edited by David Berman (dberm22@gmail.com)
% Edited by Aman Aberra and Ruochen Wang, only flags intersection if occurs within line
% segment
epsilon = 1e-5;
e1 = p1-p0;
e2 = p2-p0;
q = [d(2)*e2(3)-d(3)*e2(2), d(3)*e2(1)-d(1)*e2(3), d(1)*e2(2)-d(2)*e2(1)]; %cross product
a = e1*q'; % determinant of the matrix M
if (a>-epsilon && a<epsilon)
    % the vector is parallel to the plane (the intersection is at infinity)
    flag=0;
    u=0;
    v=0;
    t=0;
    return;
end

f = 1/a;
s = o-p0;
u = f*s*q';

if (u<0.0)
    % the intersection is outside of the triangle
    flag=0;
    u=0;
    v=0;
    t=0;
    return;
end

r = [s(2)*e1(3)-s(3)*e1(2), s(3)*e1(1)-s(1)*e1(3), s(1)*e1(2)-s(2)*e1(1)];
v = f*d*r';

if (v<0.0 || u+v>1.0)
    % the intersection is outside of the triangle
    flag=0;
    u=0;
    v=0;
    t=0;
    return;
end

t = f*e2*r'; % verified!
% flag = 1; 
% return;
if t > 0 && t < 1 % intersection occurs within line segment (ray)
   flag = 1; 
   return;
else % intersection occurs beyond end of line segment (ray)
    flag = 0;
    return; 
end
% if use_vec_mag % check that intersection occurs within line segment
%     if (t > 0 && t < 1)
%         flag = 1;
%         return; 
%     else
%         flag = 0;
%         return; 
%     end
% else
%     flag = 1;
%     return;
% end
end