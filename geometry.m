%% simplified unsteady 3D wing lifting line method by B. Davoudi
% Aerospace Engineering Department, University of Michigan 7/27/2016
function [x,y,z,xcol,ycol,zcol,n,dl_x,dly,S,alpha,cr]=geometry(AR,b,tr,Nx,Ny,Lam,dih,aoa,aoaf_L,aoaf_R)

%Nx=4;                        % chordwise panel per hald wing
%Ny=6;                        % spanwise panel per hald wing
%AR=8;                        % Aspect ratio
%b=1;                         % wing span
%tr=2;                        % taper ratio
%Lam                        % sweep angle, backward swept, positive
%dih                        % dihedral angle defined at the c/4
%aoa                        % angle of attack of the wing
%aoaf_L                     % angle of attack of the left flap
%aoa_R                      % angle of attack of the left flap

S=b^2/AR;                    % wing surface
wl=b/cos(Lam)/cos(dih);      % span length normal to the unswept wing(length of two wings) %%%%%%%%%%%%%%%
cr=2*S/wl*(1/(1+1/tr));      % chord length at the root
uinf=1;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector
hh=1/(tr-1)*wl*0.5;          % hh+b = triangle height formed by a half wing extented from the tip
dly=0.5*wl/Ny;               % pannel length span wise


% location of nodes at the quarter chord -y to y
% note leng starts from root to tip!!
if rem(2*Ny,2)>0
    leng=linspace(dly,0.5*wl-0.5*dly,Ny-0.5);
else
    leng=linspace(0.5*dly,0.5*wl-0.5*dly,Ny);
end

% chord lengths arranged from -y to y
if tr~=1
    %  c=[fliplr(cr*((hh+b/2)-linspace(0.5*dly,0.5*b-0.5*dly,Ny))/(hh+b/2)*cos(aoa)) cr*((hh+b/2)-linspace(0.5*dly,0.5*b-0.5*dly,Ny))/(hh+b/2)*cos(aoa)];
    if rem(2*Ny,2)>0
        c=[fliplr(cr*((hh+b/2)-leng)/(hh+b/2)*cos(aoa)) cr*cos(aoa) cr*((hh+b/2)-leng)/(hh+b/2)*cos(aoa)];
    else
        c=[fliplr(cr*((hh+b/2)-leng)/(hh+b/2)*cos(aoa)) cr*((hh+b/2)-leng)/(hh+b/2)*cos(aoa)];
    end
else
    % airfoil at the root
    c=repmat(cr,1,2*Ny);
end


if rem(2*Ny,2)>0
    leng2=[fliplr(-leng) 0 leng];
    leng3=[fliplr(leng) 0 leng];
else
    leng2=[fliplr(-leng) leng];
    leng3=[fliplr(leng) leng];
end

dl_x=c/Nx;                                        % pannel length stream wise
dlx=min(dl_x);                                    % used for delt t
alpha=repmat(aoa,Nx,2*Ny);

% locations of each pannel quarter chord on the wing
for i=1:Nx
    y(i,:)=leng2*cos(Lam)*cos(dih);
    % z(i,:)=leng3*sin(dih)*cos(aoa)-dl_x*(i-round(Nx/4+1))*sin(aoa)*cos(dih)
    % x(i,:)=leng3*sin(Lam) + dl_x*(i-round(Nx/4+1))*cos(aoa);
    x(i,:)=leng3*sin(Lam)-(c-dl_x)/2 *cos(aoa) + dl_x*(i-1)*cos(aoa);
    z(i,:)=leng3*sin(dih)*cos(aoa) + (c-dl_x)/2 * sin(aoa)*cos(dih) - dl_x*(i-1)*sin(aoa)*cos(dih);
end

%% collocation points only differs in x and z
xcol=x+0.5*repmat(dl_x,Nx,1).*cos(alpha);zcol=z-0.5*repmat(dl_x,Nx,1).*sin(alpha)*cos(dih);
ycol=y;

%% Fontana
% for Fontana aircraft in each wing, the flap is extended 80% started from tip.

Nf=round(0.85*Ny);
Nf=Ny;

% Nov 2017
% note that we put the bound vortex on the last wing panel, th eextension
% would be out of the wing if there is a flap, so we can modify the last
% wing panel (breal it at %75 which is hard to implement )or start the flap 
% panel at the end of the extended wing panel. this mean nothing to do for
% the flap paanel.. note that the deflection willbe taken care of in the
% fast steady code!!

for i=1:Nf
    
    
    %% flap angles of attack
    alpha(Nx,2*Ny-i+1)=aoaf_R+aoa;
    alpha(Nx,i)=aoaf_L+aoa;
     
    %% flap position

%  new way -- only collocation points will be affected

     % left
     xcol(Nx,i)=x(Nx,i)+0.5*dl_x(i)*cos(aoaf_L);
     zcol(Nx,i)=z(Nx,i)-0.5*dl_x(i)*sin(aoaf_L)*cos(dih);
     % right
     xcol(Nx,2*Ny-i+1)=x(Nx,2*Ny-i+1)+0.5*dl_x(2*Ny-i+1)*cos(aoaf_R);
     zcol(Nx,2*Ny-i+1)=z(Nx,2*Ny-i+1)-0.5*dl_x(2*Ny-i+1)*sin(aoaf_R)*cos(dih);
 

%  old way
    %right ... not neccssary
%     x(Nx,2*Ny-i+1)=x(Nx,2*Ny-i+1)+0.25*dl_x(2*Ny-i+1)*(cos(aoaf_R+aoa)-1)*cos(aoa);
%     y(Nx,2*Ny-i+1)=y(Nx,2*Ny-i+1);
%     z(Nx,2*Ny-i+1)=z(Nx,2*Ny-i+1)-0.25*dl_x(2*Ny-i+1)*(sin(aoaf_R+aoa)-sin(aoa))*cos(dih);
    
%     xcol(Nx,2*Ny-i+1)=xcol(Nx,2*Ny-i+1)+0.25*dl_x(2*Ny-i+1)*(cos(aoaf_R+aoa)-1)*cos(aoa);
%     zcol(Nx,2*Ny-i+1)=zcol(Nx,2*Ny-i+1)-0.25*dl_x(2*Ny-i+1)*(sin(aoaf_R+aoa)-sin(aoa))*cos(dih);
     
    %left ... not neccssary    
%     x(Nx,i)=x(Nx,i)+0.25*dl_x(i)*(cos(aoaf_L+aoa)-1)*cos(aoa);
%     y(Nx,i)=y(Nx,i);
%     z(Nx,i)=z(Nx,i)-0.25*dl_x(i)*(sin(aoaf_L+aoa)-sin(aoa))*cos(dih);
%     
%     xcol(Nx,i)=xcol(Nx,i)+0.25*dl_x(i)*(cos(aoaf_L+aoa)-1)*cos(aoa);
%     zcol(Nx,i)=zcol(Nx,i)-0.25*dl_x(i)*(sin(aoaf_L+aoa)-sin(aoa))*cos(dih);
    
    
end


%% used only for verification against 2D

% alpha(Nx,:)=aoaf_L+aoa;
%% normal vecotors


edger=[x(1,round(Ny+2))-x(1,round(Ny+1)) y(1,round(Ny+2))-y(1,round(Ny+1)) z(1,round(Ny+2))-z(1,round(Ny+1))];    % the c/4 vector line for the left half (+y)
edgel=[x(1,1)-x(1,2) y(1,1)-y(1,2) z(1,1)-z(1,2)];                      % the c/4 vector line for the right half (-y)
edgeb=[cr/Nx*cos(aoa) 0 -cr/Nx*sin(aoa)];                                 % root base section

% right and left normal vectors
nl=cross(edgel,edgeb)/norm(cross(edgel,edgeb));
nr=cross(edgeb,edger)/norm(cross(edgeb,edger));

%  normal to the pannels normal from -y to +y

if rem(2*Ny,2)>0
    n=repmat([repmat(nl',1,Ny-0.5) [sin(aoa); 0 ;cos(aoa)] repmat(nr',1,Ny-0.5)],1,Nx);
else
    n=repmat([repmat(nl',1,Ny) repmat(nr',1,Ny)],1,Nx);
end

% don't need n for visualization
if rem(2*Ny,2)<=0
    
    % effect of flap, note that from the TE row, any two points have the same
    % edge as the flaps. only the base edge changes....
    
    % Left
    
 edgebf=[dl_x(2*Ny-1)*cos(alpha(Nx,2*Ny-1)) 0 -dl_x(2*Ny-1)*sin(alpha(Nx,2*Ny-1))];

    
    edgelf=[x(Nx,Ny+3)-x(Nx,Ny+2) ...
        y(Nx,Ny+3)-y(Nx,Ny+2) ...
        z(Nx,Ny+3)-z(Nx,Ny+2)];
    %    n(:,2*Nx*Ny-1)=cross(edgebf,edgelf)/norm(cross(edgebf,edgelf));
    
    % Right

  edgebfr=[dl_x(2)*cos(alpha(Nx,2)) 0 -dl_x(2)*sin(alpha(Nx,2))];
      
    
    edgerf=[x(Nx,3)-x(Nx,4) ...
        y(Nx,3)-y(Nx,4) ...
        z(Nx,3)-z(Nx,4)];
    %    n(:,2*(Nx-1)*Ny+2)=cross(edgerf,edgebfr)/norm(cross(edgebfr,edgerf));
    
    %% Fontana
    
    for i=1:Nf
        n(:,2*Nx*Ny-i+1)=cross(edgebf,edgelf)/norm(cross(edgebf,edgelf));
        n(:,2*(Nx-1)*Ny+i)=cross(edgerf,edgebfr)/norm(cross(edgebfr,edgerf));
    end
    
    %% used only for verification against 2D
    %     for i=2*(Nx-1)*Ny+1:2*Nx*Ny
    %      n(:,i)=cross(edgerf,edgebfr)/norm(cross(edgebfr,edgerf));
    %     end
end
end
