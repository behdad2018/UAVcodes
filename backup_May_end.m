% steady 3D wing lifting line method by B. Davoudi 5/12/2016

clc; clear all;
Nx=4;                        % chordwise panel per hald wing
Ny=6;                        % spanwise panel per hald wing
AR=4;
b=10;
tr=1;
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoa=degtorad(5);             % angle of attack
S=b^2/AR;                    % wing surface
wl=b/cos(Lam);               % span length normal to the unswept wing(length of two wings) %%%%%%%%%%%%%%%
cr=2*S/wl*(1/(1+1/tr));      % chord length at the root
uinf=1;                      % incidence velocity
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
dl_x=c/Nx;                                        % pannel length stream wise
dlx=min(dl_x);                                    % this should be modified for tr~=1

% locations of each pannel quarter chord on the wing
for i=1:Nx
    y(i,:)=[fliplr(-leng) leng]*cos(Lam)*cos(dih);
    z(i,:)=[fliplr(leng) leng]*(sin(dih)-sin(aoa)*sin(Lam))-dlx*(i-1)*sin(aoa);
    x(i,:)=[fliplr(leng) leng]*sin(Lam)*cos(aoa) + dlx*(i-1)*cos(aoa);
end

% collocation point only differs in x and z
xcol=x+0.5*dlx*cos(aoa);zcol=z-0.5*dlx*sin(aoa);

edger=[x(1,Ny+2)-x(1,Ny+1) y(1,Ny+2)-y(1,Ny+1) z(1,Ny+2)-z(1,Ny+1)];    % the c/4 vector line for the right half (+y)
edgel=[x(1,1)-x(1,2) y(1,1)-y(1,2) z(1,1)-z(1,2)];                      % the c/4 vector line for the left half (-y)
edgeb=[cr/4*cos(aoa) 0 -cr/4*sin(aoa)];                                 % root base section

% roght and left normal vectors
nr=cross(edgeb,edger)/norm(cross(edgeb,edger));
nl=cross(edgel,edgeb)/norm(cross(edgeb,edger));

% pannels normal from -y to +y
n=repmat([repmat(nl',1,Ny) repmat(nr',1,Ny)],1,Nx);

% RHS
% maximum simulation time
tmax=25;

% time step
dt=1/16*dlx*Nx/uinf;
% dlw=dt*uinf;

% firt wake row locations -- the leading element of the wake vortex ring measure from TE pannel
xw=repmat(x(Nx,:)+dlx*cos(aoa),tmax-1,1);
yw=repmat(y(Nx,:),tmax-1,1);
zw=repmat(z(Nx,:)-dlx*sin(aoa),tmax-1,1);
wake=0; % u=[0 0 0];

for t=1:tmax
    t
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
            [a9,A(i,j)]=vortexring(n(:,i),dlx,dly,aoa,Lam,dih,xcol1,ycol1,zcol1,x1,y1,z1,1);       % Effect of ring vortices on the wing
        end
        % panels on the wake on the pannel on the win
        
        if t>1
            Nw=t-1;
            wake2=0;a2=0;
            for j=1:2*Nw*Ny
                index1w=ceil(j/(2*Ny));
                index2w=j-(2*Ny)*(index1w-1);
                xw1=xw(index1w,index2w);      yw1=yw(1,index2w);      zw1=zw(1,index2w);
                [a1,wake1]=vortexring(n(:,i),dlw1(index1w,index2w),dly,0,0,0,xcol1,ycol1,zcol1,xw1,yw1,zw1,Gw(j));
                wake2=wake2+wake1;
                a2=a2+a1;
            end
            wake(i)=wake2;
            a(:,i)=a2;
        end
    end
    %% RHS = RHS-wake'
   
    RHS=dot(-repmat(u',1,2*Nx*Ny),n)'-wake';
    G(:,t)=linsolve(A,RHS);                                    % vorticity on the wing
    Gw(1+2*(t-1)*Ny:2*t*Ny)=G(1+2*(Nx-1)*Ny:2*Nx*Ny,t);        % vorticity on the wake i.e., its strength the same as TE pannel
    
    if t>1
        uw=wakevelocity(Nx,Ny,Nw,dlw1,dlx,dly,aoa,Lam,dih,x,y,z,xw,yw,zw,G(:,t),Gw);
        uw1=reshape(uw(1,:,:),Nw,2*Ny);
        vw1=reshape(uw(2,:,:),Nw,2*Ny);
        ww1=reshape(uw(3,:,:),Nw,2*Ny);
        dlw1(1:t-1,:)=(uw1+uinf)*dt;
       xw(1:t-1,:)=xw(1:t-1,:)+(uw1+uinf)*dt;
    %  xw(1:t-1,:)=bsxfun(@minus,xw(1:t-1,:)+(uw1+uinf)*dt,0.7*uinf*dt);
        yw(1:t-1,:)=yw(1:t-1,:)+vw1*dt;
        zw(1:t-1,:)=zw(1:t-1,:)+ww1*dt;
      
    end
    
    if t==1
        a=zeros(3,2*Nx*Ny);
        uw11(1,:)=zeros(1,2*Ny);
        dlw1(t,:)=(uw11+uinf)*dt;
    else
      dlw1(t,:)=(uw11+uinf)*dt;
    end
    %% Force calculations
    
    for j=1:2*Nx*Ny
        
        index1=ceil(j/(2*Ny));
        index2=j-(2*Ny)*(index1-1);
        
        G1(index1,index2)=G(j,t)';
        
        if t>1
            G2(index1,index2)=G(j,t-1)';
        else
            G2(index1,index2)=0;
        end
        u_ind(index1,index2)=a(1,j);
        v_ind(index1,index2)=a(2,j);
        w_ind(index1,index2)=a(3,j);
    end
    
    dGdx=[G1(1,:);diff(G1,1,1)]/dlx;
    dGdy=[zeros(Nx,1),diff(G1(:,Ny+1:2*Ny),1,2)/dly];
    dgdt=(G1-G2)/dt;
    
    dp= (u_ind +repmat(uinf,Nx,2*Ny)).*[repmat(cos(aoa),Nx,Ny) repmat(cos(aoa),Nx,Ny)] .* dGdx + ...
        repmat(v_ind(:,Ny+1:2*Ny).*repmat(cos(Lam)*cos(dih),Nx,Ny).*  dGdy,1,2) +...
        dgdt;
    
    F=sum(sum((dp*dlx*dly)))
    cLnew(t)=F*cos(aoa)/(0.5*S*uinf^2);
    
    cL(t)=2*sum(G(:,t))*dlx*dly*wl/S    /(uinf*S);
    
    u=uinf*[1 0 0];
end

figure(1)
% plot(G); hold on;
% figure
plot(linspace(0,tmax*dt*uinf/cr,tmax),cLnew);xlabel('time'); ylabel('C_L');grid minor;hold on;axis([0 9 0 0.55]);
hold on;
% plot(y(1,:),G(1:2*Ny)')
%legend('AR = 4','AR = 8','AR = 12','AR = 20','AR = \inf');
figure(3)
surf(x,y,G2)