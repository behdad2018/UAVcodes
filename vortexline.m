function [vin,vn]=vortexline(n,xa,ya,za,x1,y1,z1,x2,y2,z2,g)

% N : number of pannels in each half wing
% n : nomral vectors for the pannels
% xa,ya,za  location of the collocartion point
% x1,y1,z1 start of a vortex line
% x2,y2,z2 end of a vortex line
% once g=1, it find the coefficient that mutiplies to Gamma

    r12=[x2-x1 y2-y1 z2-z1];   % vortex vector
    r1a=[xa-x1 ya-y1 za-z1];   % distance vector from the begining of the vortex to the collocation point
    r2a=[xa-x2 ya-y2 za-z2];   % distance vector from the end of the vortex to the collocation point

% cosine of the angles
b1=dot(r12,r1a)/(norm(r12)*norm(r1a));
b2=dot(r12,r2a)/(norm(r12)*norm(r2a));

% normal distance of a line vortex to the collocation point
h=abs(norm(r1a)*sin(pi-acos(b1)));

dir1=cross(r12,r1a);
dir=dir1/norm(dir1);         % direction of the induced velocity

vin=g/(4*pi*h)*abs((b1-b2))*dir;  % induced velocity vector
%vin=g/(4*pi*h)*(b1-b2)*dir;      % induced velocity vector

% normal velocity to the pannel, note n(:,1) is the pannel normal of -y
% half wing and n(:,N+1) is the pannel normal of +y half wing

vn=dot(vin,n');


