function p3 = closestPointOnLine(p1,p2,q)
%CLOSESTPOINTONLINE For query point q, finds closest point on line defined
%by two points p1 and p2
%  
%   Inputs 
%   ------ 
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 

u = p2 - p1; 
pq = q - p1; 
w2 = pq - dot(pq,u)*u/norm(u)^2; 
p3 = q - w2; 