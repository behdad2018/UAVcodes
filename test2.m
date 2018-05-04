%%

% vortexrim([0,0,1],0,0,0,1,0,0,1,0,0,1)
%%
tic
clc; clear all;
%% wing config
Nxw=8;Nxt=2;               % chordwise panel per hald wing
Nyw=10;Nyt=Nyw/2;               % spanwise panel per hald wing
ARw=176.5^2/7150;             % Aspect ratio
ARt=3.5;
bw=1;bt=bw*64.4/176.5;         % wing span
%bt=1;ARt=ARw;
trw=1.8;trt=2;              % taper ratio
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoa=degtorad(7);             % angle of attack
uinf=1;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector
aoaf_Lw=degtorad(0);aoaf_Lt=degtorad(0);
aoaf_Rw=degtorad(0);aoaf_Rt=degtorad(0);

%% propeller config

nb=2;                   % number of blade
rpm=5000;               % rpm
R=0.1;                  % propeller radius
c=0.025;                % propeller chord
ro=1;                   % density
uinf=1;                 % free stream velocity

A=pi*R^2;               % disk area
sig=nb*c/(pi*R);        % solidity
om=rpm/60;            % rev per second
gam=ro*2*pi*c*R^4/1;    % lock number
Vt=om*R;                % tip velocity
alpha_p=degtorad(90);   % for a propeller it is close to 90 deg
num1=6;                 % this is the number of sections in each circle
num2=4;                 % number of rings (sqrt must be an integer)
rnum=10;                % number of cylinders

% location of the propeller
xp=-0.15;yp=0;zp=0;

[xs,ys,zs,r,gx,gt,gd,vinf]=propeller(bw,xp,yp,zp,num1,num2,rnum,nb,rpm,R,c,uinf,alpha_p,uinf);

% lifting line location
xll=xp-3*c/4;
yll(:)=ys(2,:,end);
zll(:)=zs(2,:,end);

% control point locations for the propeller
xcolp=xll+c/2;
ycolp=yll+c/2;
zcolp=zll+c/2;

%%
for i=1:num1-1     % number of points in a cricle (diff azimuth angle) , i.e. number of collacation points

    xcol1=xcolp;    ycol1=ycolp(i);    zcol1=zcolp(i);
    
    for j=1:num1-1 
       
         [~,Awing(i,j)]=vortexring
        
    end
    
    for k=1:num2
        
        for j=1:rnum
            
            [Vin,Vn]=vortexrim(n,xcol1,ycol1,zcol1,x,y,z,r0,vryaw,vrpitch,om)
             
        end
    end   
end


