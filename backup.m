% 3D wing lifting line method
clc;clear all;close all;

N=20;                        % number of pannels in each half wing
Lam=degtorad(0);           % sweep angle, backward swept, positive
dih=degtorad(0);           % dihedral angle defined at the c/4
aoa=degtorad(5)/cos(Lam);   % angle of attack 
tr=4;                       % taper ratio
cr=1;                       % chord length at the root
b=10/cos(Lam);                       % span length
AR=b^2/(2*(cr+cr/tr)/2*b/2);% Aspect ratio
uinf=1*[cos(aoa) 0 sin(aoa)];               % incidence velocity
hh=1/(tr-1)*b*0.5;

dl=0.5*b/N;                 % pannel length

% chord lengths arranged from -y to y

c=[fliplr(cr*((hh+b/2)-linspace(0.5*dl,0.5*b-0.5*dl,N))/(hh+b/2)*cos(aoa)) cr*((hh+b/2)-linspace(0.5*dl,0.5*b-0.5*dl,N))/(hh+b/2)*cos(aoa)];

% location of nodes at the quarter chord -y to y

y=[fliplr(-linspace(0.5*dl,0.5*b-0.5*dl,N)*cos(Lam)*cos(dih)) linspace(0.5*dl,0.5*b-0.5*dl,N)*cos(Lam)*cos(dih)];
z=[fliplr(linspace(0.5*dl,0.5*b-0.5*dl,N)*sin(dih)) linspace(0.5*dl,0.5*b-0.5*dl,N)*sin(dih)]-c/4*sin(aoa);
x=[fliplr(0.25*cr+linspace(0.5*dl,0.5*b-0.5*dl,N)*sin(Lam)*cos(aoa)) 0.25*cr+linspace(0.5*dl,0.5*b-0.5*dl,N)*sin(Lam)*cos(aoa)];

% collocation point only differs in x and z
xcol=x+0.5*c*cos(aoa);zcol=z-0.5*c*sin(aoa);

edger=[x(N+2)-x(N+1) y(N+2)-y(N+1) z(N+2)-z(N+1)];    % the c/4 vector line for the right half (+y)
edgel=[x(1)-x(2) y(1)-y(2) z(1)-z(2)];                % the c/4 vector line for the left half (-y)
edgeb=[cr/4*cos(aoa) 0 -cr/4*sin(aoa)];               % root base section

nr=cross(edgeb,edger)/norm(cross(edgeb,edger));
nl=cross(edgel,edgeb)/norm(cross(edgeb,edger));

% pannels normal from -y to +y
n=[repmat(nl',1,N) repmat(nr',1,N)];

% RHS
RHS=dot(-repmat(uinf',1,2*N),n);

for i=1:2*N             % collocation index
    for j=1:2*N         % pannel index
        j
        % left trailing vortex (shorter)
        x1=x(j)+0.5*dl*sin(Lam)*cos(aoa);
        z1=z(j)+0.5*dl*sin(dih)-c(i)/4
        0.5*dl*cos(Lam)*sin(aoa);
        % right trailing vortex (longer)
        x2=x(j)-0.5*dl*sin(Lam)*cos(aoa);
        z2=z(j)-0.5*dl*sin(dih)+0.5*dl*cos(Lam)*sin(aoa);
        if y(j)<0
            y1=y(j)-0.5*dl*cos(Lam)*cos(dih);
            y2=y(j)+0.5*dl*cos(Lam)*cos(dih);
        elseif y(j)>0
            y1=y(j)+0.5*dl*cos(Lam)*cos(dih);
            y2=y(j)-0.5*dl*cos(Lam)*cos(dih);
%             temp1=x2;
%             x2=x1;
%             x1=temp1;
%             temp2=z2;
%             z2=z1;
%             z1=temp2;
        end


%         if y(j)<0
%             y1=y(j)-0.5*dl*cos(Lam)*cos(dih);
%             y2=y(j)+0.5*dl*cos(Lam)*cos(dih);
%                     x1=x(j)+0.5*dl*sin(Lam)*cos(aoa);
%         z1=z(j)+0.5*dl*sin(dih)-0.5*dl*cos(Lam)*sin(aoa);
%         % right trailing vortex (longer)
%         x2=x(j)-0.5*dl*sin(Lam)*cos(aoa);
%         z2=z(j)-0.5*dl*sin(dih)+0.5*dl*cos(Lam)*sin(aoa);
%         elseif y(j)>0
%             y1=y(j)+0.5*dl*cos(Lam)*cos(dih);
%             y2=y(j)-0.5*dl*cos(Lam)*cos(dih);
%         x1=x(j)-0.5*dl*sin(Lam)*cos(aoa);
%         z1=z(j)-0.5*dl*sin(dih)+0.5*dl*cos(Lam)*sin(aoa);
%         % right trailing vortex (longer)
%         x2=x(j)+0.5*dl*sin(Lam)*cos(aoa);
%         z2=z(j)+0.5*dl*sin(dih)-0.5*dl*cos(Lam)*sin(aoa);
% 
%         end
        A_b(i,j)=vortexline(N,n,xcol(i),y(i),zcol(i),x1,y1,z1,x2,y2,z2);
        A_tv1(i,j)=vortexline(N,n,xcol(i),y(i),zcol(i),20*b,y1,z1,x1,y1,z1);
        A_tv2(i,j)=vortexline(N,n,xcol(i),y(i),zcol(i),x2,y2,z2,20*b,y2,z2);
        
    end
end

A= A_tv1+ A_tv2+ A_b;
G=A^-1*RHS';

plot(y,G)
   
%         xtest2(nnn)=x2;
%         ztest2(nnn)=z2;
%         ztest1(nnn)=z1;
%         xtest1(nnn)=x1;
%         ytest1(nnn)=y1;
%         ytest2(nnn)=y2;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5/23 before wake velocity

% steady 3D wing lifting line method by B. Davoudi 5/12/2016

clc; clear all;
Nx=4;                        % chordwise panel per hald wing
Ny=7;                        % spanwise panel per hald wing
AR=4;
b=10;
tr=1;
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoa=degtorad(5);             % angle of attack
S=b^2/AR;                    % wing surface
wl=b/cos(Lam);               % span length normal to the unswept wing(length of two wings) %%%%%%%%%%%%%%%
cr=2*S/wl*(1/(1+1/tr));      % chord length at the root
uinf=10;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector
hh=1/(tr-1)*wl*0.5;          % hh+b = triangle height formed by a half wing extented from the tip

dly=0.5*wl/Ny;               % pannel length span wise

% chord lengths arranged from -y to y
if tr~=1
    c=[fliplr(cr*((hh+b/2)-linspace(0.5*dly,0.5*b-0.5*dly,Ny))/(hh+b/2)*cos(aoa)) cr*((hh+b/2)-linspace(0.5*dly,0.5*b-0.5*dly,Ny))/(hh+b/2)*cos(aoa)];
else
    c=repmat(cr,1,2*Ny);
end

% location of nodes at the quarter chord -y to y
leng=linspace(0.5*dly,0.5*wl-0.5*dly,Ny);
dl_x=c/Nx;               % pannel length stream wise
dlx=min(dl_x);           % this should be modified for tr~=1

% locations of each pannel quarter chord on the wing
for i=1:Nx
    y(i,:)=[fliplr(-leng) leng]*cos(Lam)*cos(dih);
    z(i,:)=[fliplr(leng) leng]*(sin(dih)-sin(aoa)*sin(Lam)) - dlx*(i-1)*sin(aoa);
    x(i,:)=[fliplr(leng) leng]*sin(Lam)*cos(aoa) + dlx*(i-1)*cos(aoa);
end

% collocation point only differs in x and z
xcol=x+0.5*dlx*cos(aoa);zcol=z-0.5*dlx*sin(aoa);

edger=[x(1,Ny+2)-x(1,Ny+1) y(1,Ny+2)-y(1,Ny+1) z(1,Ny+2)-z(1,Ny+1)];    % the c/4 vector line for the right half (+y)
edgel=[x(1,1)-x(1,2) y(1,1)-y(1,2) z(1,1)-z(1,2)];                % the c/4 vector line for the left half (-y)
edgeb=[cr/4*cos(aoa) 0 -cr/4*sin(aoa)];               % root base section

nr=cross(edgeb,edger)/norm(cross(edgeb,edger));
nl=cross(edgel,edgeb)/norm(cross(edgeb,edger));

% pannels normal from -y to +y
n=repmat([repmat(nl',1,Ny) repmat(nr',1,Ny)],1,Nx);

% RHS
tmax=20;
% firt wake row locations

dt=10*dlx/uinf;
% xw=zeros(20,2*Ny);
% xw(1,:)=x(Nx,:)+dlx*cos(aoa);
% yw(1,:)=y(Nx,:);
% zw(1,:)=z(Nx,:)-dlx*cos(aoa);
xw=repmat(x(Nx,:)+dlx*cos(aoa),tmax-1,1);
yw=repmat(y(Nx,:),20,1);
zw=repmat(z(Nx,:)-dlx*sin(aoa),tmax-1,1);

wake=0; %u=[0 0 0];

for t=1:tmax
    t
    % RHS=dot(-repmat(u',1,2*Nx*Ny),n)'-wake';
    u=uinf*[1 0 0];
    %     if t==10
    %         u=uinf*[0 0 0];
    %     end
    %  time=(t-1)*dt
    if t>2
        xw(1:t-2,:)=xw(1:t-2,:)+uinf*dt;
    end
    for i=1:2*Nx*Ny
        
        % panels interations on the wing
        for j=1:2*Nx*Ny
            
            index1col=ceil(i/(2*Ny));
            index2col=i-(2*Ny)*(index1col-1);
            index1=ceil(j/(2*Ny));
            index2=j-(2*Ny)*(index1-1);
            x1=x(index1,index2);
            y1=y(index1,index2);
            z1=z(index1,index2);
            xcol1=xcol(index1col,index2col);
            ycol1=y(index1col,index2col);
            zcol1=zcol(index1col,index2col);
            
            A(i,j)=vortexring(n(:,i),dlx,dly,aoa,Lam,dih,xcol1,ycol1,zcol1,x1,y1,z1,1);       % Effect of ring vortices on the wing
            
        end
        
        % panels on the wake
        if t>1
            Nw=t-1;
            wake2=0;
            for j=1:2*Nw*Ny
                index1w=ceil(j/(2*Ny));
                index2w=j-(2*Ny)*(index1w-1);
                xw1=xw(index1w,index2w);      yw1=yw(1,index2w);      zw1=zw(1,index2w);
                
                wake1=vortexring(n(:,i),dt*uinf,dly,0,0,0,xcol1,ycol1,zcol1,xw1,yw1,zw1,Gw(j));
                wake2=wake2+wake1;              
            end
            wake(i)=wake2;
        end
    end
    for k=1:2*Nw*Ny
        index1w=ceil(k/(2*Ny));
        index2w=k-(2*Ny)*(index1w-1);
        for j=1:2*Nw*Ny
            index1w=ceil(j/(2*Ny));
            index2w=j-(2*Ny)*(index1w-1);
            xw1=xw(index1w,index2w);      yw1=yw(1,index2w);      zw1=zw(1,index2w);
            
            index1w=ceil(j/(2*Ny));
            index2w=j-(2*Ny)*(index1w-1);
            xw1=xw(index1w,index2w);      yw1=yw(1,index2w);      zw1=zw(1,index2w);
            
            vel_wake(k)=vortexring(n(:,i),dt*uinf,dly,0,0,0,xcol1,ycol1,zcol1,xw1,yw1,zw1,Gw(j))
            
        end
    end
    
    % RHS=RHS-wake';
    RHS=dot(-repmat(u',1,2*Nx*Ny),n)'-wake';
    G=linsolve(A,RHS);
    %    Gw(1+2*(t-1)*Ny:2*t*Ny)=flipud(G(1+2*(Nx-1)*Ny:2*Nx*Ny));
    Gw(1+2*(t-1)*Ny:2*t*Ny)=G(1+2*(Nx-1)*Ny:2*Nx*Ny);
    
    %cL(t)=sum(2*G(1:2*Ny)-G(1+2*(Nx-1)*Ny:2*Nx*Ny));
    % cL(t)=2*sum(G)*dly/(uinf*S);
    %cL(t)=2*trapz(y,G)/(uinf*S);
    
    cL(t)=2*sum(G)*dlx*dly*wl/S    /(uinf*S);
    
    sum(G);
    %   cl_cL=2*G./(uinf*c')./cL(t);
    
end
plot(G);hold on
%plot(cL);xlabel('time'); ylabel('C_L');grid minor;hold on;
% plot(y(1,:),G(1:2*Ny)')

%%%%%%%%%%%  Force calculation
%  dp= dot([u_ind + repmat(uinf,Nx,2*Ny) ; v_ind; w_ind], [repmat(cos(aoa),Nx,Ny) repmat(cos(aoa),Nx,Ny);...
%          repmat(0,Nx,Ny) repmat(0,Nx,Ny);...
%          repmat(0,Nx,Ny) repmat(0,Nx,Ny)]) * dGdx + ...
%        dot([u_ind + repmat(uinf,Nx,2*Ny) ; v_ind; w_ind], [repmat(-sin(aoa),Nx,Ny) repmat(-sin(aoa),Nx,Ny);...
%          repmat(0,Nx,Ny) repmat(0,Nx,Ny);...
%          repmat(0,Nx,Ny) repmat(0,Nx,Ny)]) *  dGdy +...
%         diff(G(:,t))/dt;
