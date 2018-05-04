% steady 3D wing lifting line method by B. Davoudi 5/12/2016
% This code uses Ligting line theory and compute vorticity
% distribution for a thin trapezoidal wing. Number of pannels
% Taper ratio, total wing length, incidence velocity magnitude and
% chord length at the root are the inputs as well as angle of
% attack, dihedral and sweep angles. The vortexline function
% calculates the induced velcoity normal to each collocation
% point assuming the Gamma is one

function [y,cL,cl_cL,G]=L_method(N,tr,b,AR,aoa,Lam,dih)

% N                            % number of pannels in each half wing
% tr                           % taper ratio
Lam=degtorad(Lam);             % sweep angle, backward swept, positive
dih=degtorad(dih);             % dihedral angle defined at the c/4
aoa=degtorad(aoa);
%/cos(Lam);  % angle of attack
S=b^2/AR;                       % wing surface
wl=b/cos(Lam);                 % span length normal to the unswept wing(length of two wings) %%%%%%%%%%%%%%%
cr=2*S/wl*(1/(1+1/tr));        % chord length at the root
% S=2*(cr+cr/tr)/2*wl/2;       % wing surface
% cr=S/b                      % wing surface, not correct ... Bertin's book used for verification

uinf=1;                       % incidence velocity
u=uinf*[1 0 0];               % incidence velocity vector
hh=1/(tr-1)*wl*0.5;           % hh+b = triangle height formed by a half wing extented from the tip

dl=0.5*wl/N;                  % pannel length

% chord lengths arranged from -y to y
if tr~=1
    c=[fliplr(cr*((hh+b/2)-linspace(0.5*dl,0.5*b-0.5*dl,N))/(hh+b/2)*cos(aoa)) cr*((hh+b/2)-linspace(0.5*dl,0.5*b-0.5*dl,N))/(hh+b/2)*cos(aoa)];
else
    c=repmat(cr,1,2*N);
end

% location of nodes at the quarter chord -y to y
leng=linspace(0.5*dl,0.5*wl-0.5*dl,N);
y=[fliplr(-leng) leng]*cos(Lam)*cos(dih);
z=[fliplr(leng) leng]*(sin(dih)-sin(aoa)*sin(Lam));
x=[fliplr(leng) leng]*sin(Lam)*cos(aoa);

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
RHS=dot(-repmat(u',1,2*N),n);
nnn=0;

for i=1:2*N             % collocation index
    
    for j=1:2*N         % pannel index
        nnn=nnn+1;
        % left trailing vortex
        y1=y(j)-0.5*dl*cos(Lam)*cos(dih);
        
        % right trailing vortex
        y2=y(j)+0.5*dl*cos(Lam)*cos(dih);
        
        if y(j)<0
            x1=x(j)+0.5*dl*sin(Lam)*cos(aoa);
            z1=z(j)+0.5*dl*(sin(dih)-sin(aoa)*sin(Lam));
            x2=x(j)-0.5*dl*sin(Lam)*cos(aoa);
            z2=z(j)-0.5*dl*(sin(dih)-sin(aoa)*sin(Lam));
%             x3=x1+0.75*c(j)*cos(aoa);
%             x4=x2+0.75*c(j)*cos(aoa);
%             z3=z1-0.75*c(j)*sin(aoa);
%             z4=z2-0.75*c(j)*sin(aoa);
             x3=x1+c(j)*cos(aoa);
             x4=x2+c(j)*cos(aoa);
             z3=z1-c(j)*sin(aoa);
             z4=z2-c(j)*sin(aoa);
        elseif y(j)>0
            x1=x(j)-0.5*dl*sin(Lam)*cos(aoa);
            z1=z(j)-0.5*dl*(sin(dih)-sin(aoa)*sin(Lam));
            x2=x(j)+0.5*dl*sin(Lam)*cos(aoa);
            z2=z(j)+0.5*dl*(sin(dih)-sin(aoa)*sin(Lam));
%             x3=x1+0.75*c(j)*cos(aoa);
%             x4=x2+0.75*c(j)*cos(aoa);
%             z3=z1-0.75*c(j)*sin(aoa);
%             z4=z2-0.75*c(j)*sin(aoa);
             x3=x1+c(j)*cos(aoa);
             x4=x2+c(j)*cos(aoa);
             z3=z1-c(j)*sin(aoa);
             z4=z2-c(j)*sin(aoa);
        end
        % trailing wakes lie on the wing and then become parallel to x after leaving the TE
        
        [a, A_b(i,j)]=vortexline(n(:,i),xcol(i),y(i),zcol(i),x1,y1,z1,x2,y2,z2,1);       % Effect of bound vortex
        [a, A_tv1w(i,j)]=vortexline(n(:,i),xcol(i),y(i),zcol(i),x3,y1,z3,x1,y1,z1,1);    % Effect of left trailing vortex on the wing
        [a, A_tv2w(i,j)]=vortexline(n(:,i),xcol(i),y(i),zcol(i),x2,y2,z2,x4,y2,z4,1);    % Effect of right trailing vortex on the wing
        [a, A_tv1(i,j)]=vortexline(n(:,i),xcol(i),y(i),zcol(i),20*b,y1,z3,x3,y1,z3,1);   % Effect of left trailing vortex
        [a, A_tv2(i,j)]=vortexline(n(:,i),xcol(i),y(i),zcol(i),x4,y2,z4,20*b,y2,z4,1);   % Effect of right trailing vortex
        
  %          A_test(i,j)=vortexline(n(:,i),xcol(i),y(i),zcol(i),x4,y2,z4,x3,y1,z3,1);    % Effect of right trailing vortex on the wing
%         Atest1(i,j)=vortexring(n(:,i),cr,dl,aoa,Lam,dih,xcol(i),y(i),zcol(i),x(j),y(j),z(j),1);        
%         Atest2(i,j)=vortexring(n(:,i),20*b,dl,0,0,0,xcol(i),y(i),zcol(i),0.5*(x3+x4),0.5*(y1+y2),0.5*(z3+z4),1);
        % trailing wake is always parallel to the x
        
%         A_b(i,j)=vortexline(N,n,xcol(i),y(i),zcol(i),x1,y1,z1,x2,y2,z2);
%         A_tv1(i,j)=vortexline(N,n,xcol(i),y(i),zcol(i),20*b,y1,z1,x1,y1,z1);
%         A_tv2(i,j)=vortexline(N,n,xcol(i),y(i),zcol(i),x2,y2,z2,20*b,y2,z2);     
    end
end

A   = A_tv1 + A_tv2 + A_b + A_tv1w + A_tv2w ;
% Ain = A_tv1 + A_tv2 + A_b + A_tv1w + A_tv2w ;
% A= A_tv1+ A_tv2+ A_b; 
% A= Atest1+ Atest2;
% A   = A_tv1w + A_tv2w + A_b + A_test ;
% vorticity distrubution on the wing..note that the trailing vortices and
% bound vortex of a horseshow vortex element have the same magnitude\


G=linsolve(A,RHS');

cL=2*trapz(y,G)/(uinf*S);
% cL=2*sum(G)*dl*cos(Lam)/(uinf*S);  Bertin's book used for verification
cl_cL=2*G./(uinf*c')./cL;

end
