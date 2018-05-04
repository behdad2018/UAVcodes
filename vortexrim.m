function [Vin,Vn]= vortexrim(n,xcol1,ycol1,zcol1,x,y,z,r0,vryaw,vrpitch,g)
% written by B. Davoudi, University of Michigan
% this code will find the induced velocity of a ring of the panel 
% the inputs are locations of the control points: xcol2, ycol1 and zcol1
% the locations of the center of the ring are x,y and z
% r0 is the radiaus of the ring
% g is the magnitude of the vorticicty which is constant along the vortex
% ring the vortex ring moves in x -- ******** in mannuaver, we need the
% angle at which ring is released and can rotate the induced velocity with
% that angle to result in the "real" coordinate
% vryaw is the vortex ring yaw at the time of relaease
% vrpitch is the vortex ring pitch at the time of relaease
% the vector from the center of the ring to the coll point

cen=[xcol1-x,ycol1-y,zcol1-z];

% introduced for integration
th=0:0.1:2*pi;

% the vector from the element of the vortex ring to the coll point
r=repmat(cen,size(th,2),1)+r0*[zeros(size(th,2),1),-cos(th'),-sin(th')];

% the direction of the vortex element chosen on the ring
dl=[zeros(size(th,2),1),-sin(th'),cos(th')];

% absolute value of the r
absr=sqrt(sum(abs(r).^2,2));

% induced velocity vector on the panel
Vv1=g*r0/(4*pi)*trapz(th',bsxfun(@rdivide,cross(dl,r),absr.^3))

Ry=[cos(vrpitch),0,sin(vrpitch);0,1,0;-sin(vryaw),0,cos(vryaw)];

Rz=[cos(vryaw),-sin(vryaw),0;sin(vryaw),cos(vryaw),0;0,0,1];

vin=Ry*Ry*Vv1'

% the normal induced velocity
vn=dot(vin,n)

end