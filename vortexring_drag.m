
function [vin,vn]=vortexring_drag(n,dlx,dly,aoa,Lam,dih,xa,ya,za,x,y,z,g,theta)

% this code only considerss the trailing vortices that are used in the
% drag calculation

% n : nomral vectors for the pannels
% xa ya za -- location of the collocation point
% x y z -- locartion of quarter of a panel or l=the lading line of a vortex
% ring

% right vortex
y1=y-0.5*dly*cos(Lam)*cos(dih);

% left vortex
y2=y+0.5*dly*cos(Lam)*cos(dih);
% 1 and 2 are the left and right neighbors (j-direction) and 3 and 4 are the ones
% downstream (i-direction)

if y<0
    x1=x+0.5*dly*sin(Lam)*cos(aoa);
    z1=z+0.5*dly*(sin(dih)-sin(aoa)*sin(Lam));
    x2=x-0.5*dly*sin(Lam)*cos(aoa);
    z2=z-0.5*dly*(sin(dih)-sin(aoa)*sin(Lam));
    x3=x1+dlx*cos(aoa);
    y3=y1;
    x4=x2+dlx*cos(aoa);
    z3=z1-dlx*sin(aoa)*cos(dih);
    z4=z2-dlx*sin(aoa)*cos(dih);
    y4=y2;
    
elseif y>=0
    x1=x-0.5*dly*sin(Lam)*cos(aoa);
    z1=z-0.5*dly*(sin(dih)-sin(aoa)*sin(Lam));
    x2=x+0.5*dly*sin(Lam)*cos(aoa);
    z2=z+0.5*dly*(sin(dih)-sin(aoa)*sin(Lam));
    x3=x1+dlx*cos(aoa);
    y3=y1;
    x4=x2+dlx*cos(aoa);
    z3=z1-dlx*sin(aoa)*cos(dih);
    z4=z2-dlx*sin(aoa)*cos(dih);
    y4=y2;
end

%% if the pannel makes an angle with the axis of airplane this angle
% is defined with positive y -- which is the left wing looking at
% the UAV from the front
% also note that collocation points (xa,ya,za) have already been updated
% outside this function
if theta~=0
    
    A=[1,0,0;0,cosd(theta),-sind(theta);0,sind(theta),cosd(theta)]*[x1;y1;z1];
    x1=A(1);y1=A(2);z1=A(3);
    
    B=[1,0,0;0,cosd(theta),-sind(theta);0,sind(theta),cosd(theta)]*[x2;y2;z2];
    x2=B(1);y2=B(2);z2=B(3);
    
    C=[1,0,0;0,cosd(theta),-sind(theta);0,sind(theta),cosd(theta)]*[x3;y3;z3];
    x3=C(1);y3=C(2);z3=C(3);
    
    D=[1,0,0;0,cosd(theta),-sind(theta);0,sind(theta),cosd(theta)]*[x4;y4;z4];
    x4=D(1);y4=D(2);z4=D(3);
  
end
%% for drag: only streamwise contributions are important, so 1-2 ans 3-4 must be excluded.

%[vin1,vn1]=vortexline(n,xa,ya,za,x1,y1,z1,x2,y2,z2,g);
[vin2,vn2]=vortexline(n,xa,ya,za,x3,y3,z3,x1,y1,z1,g);
[vin3,vn3]=vortexline(n,xa,ya,za,x2,y2,z2,x4,y4,z4,g);
%[vin4,vn4]=vortexline(n,xa,ya,za,x4,y4,z4,x3,y3,z3,g);

%%
vn1=0;vn4=0;vin1=0;vin4=0;

vn=vn1+vn2+vn3+vn4;
vin=vin1+vin2+vin3+vin4;

end