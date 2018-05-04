% propeller code
% using lifting line to find vorticity distruntion

function [x,y,z,r,gx,gt,gd,cl_asitav]=propeller(b,xp,yp,zp,num1,num2,rnum,nb,rpm,R,cm,uinf,alpha,V,F)

% b=10;                   % win span
% nb=2;                   % number of blade
% rpm=5000;               % rpm
% R=0.2;                  % propeller radius
% c=0.05;                 % propeller chord
% V=10;                   % the velocity should be given
% uinf=1;                 % free stream velocity
% alpha=degtorad(90);     % for a propeller it is close to 90 deg
% num1=5;                 % this is the number of sections in each circle
% num2=4;                 % number of rings (swrt must be an integer)
% rnum                    % number of cylinders
% xp,yp,zp                % location of the propeller

A=pi*R^2;               % disk area
sig=nb*cm/(pi*R);        % solidity
om1=rpm/60;             % rev per second
om=rpm*2*pi/60;         % rad per second
gam=R*2*pi*cm*R^4/1;     % lock number
Vt=om*R;                % tip velocity

lamc=V/Vt;              % climb inflow ratio
mu=V*cos(alpha)/Vt;     % advance ratio

%c=[linspace(cm/2,cm,floor(rnum/2)) linspace(cm,cm/2,rnum-floor(rnum/2))];

r=linspace(1/rnum,1,rnum);        % propeller radial
%cla=repmat(1.7*pi,1,rnum);        % lift curve slope
%cla=linspace(0.8*pi,1.5*pi,rnum);        % lift curve slope
% th=repmat(degtorad(30),1,rnum);   % the pitch angle
% th=[repmat(degtorad(60),1,floor(rnum/8)) linspace(degtorad(50),degtorad(30),rnum-floor(rnum/8))];
% th=linspace(degtorad(70),degtorad(50),rnum);
% from the thesis info

% th0=[10,5];
% r0=[.65,1];the0=[repmat(10,1,12) interp1(r0,th0,r(13:end))];
% cl0=[1.7/degtorad(20) 1.2/degtorad(12)];
% cla=[repmat(1.7/degtorad(20),1,12) interp1(r0,cl0,r(13:end))];

% used in conference 2017 aviation 
% th0=[6.3,1];
% r0=[0,1];the0=interp1(r0,th0,r);
% cl0=[1.5*pi 1.7*pi];
% %[1.2/degtorad(14) 1.1/degtorad(12)];
% cla=interp1(r0,cl0,r);
% 

% CFD 2-D data used in journal 

the0=[ 1.73, 2.5, 3.19, 3.8, 4.33, 4.77, 5.14, 5.42, 5.62, 5.73, 5.77, 5.72, 5.59, 5.38, 5.08, 4.71, 4.25, 3.71, 3.09, 2.38];
cla=[ 5.21, 5.19, 5.17, 5.16, 5.15, 5.15, 5.15, 5.15, 5.15, 5.16, 5.17, 5.18, 5.2, 5.22, 5.24, 5.26, 5.29, 5.32, 5.36, 5.39];

th=degtorad(the0+[56.7,56.7,56.7,56.7, 51.3, 47.0, 43.5, 40.7, 38.3, 36.3, 34.5, 32.8, 31.3, 29.8, 28.3, 26.9, 25.5, 24.0, 22.4, 20.8]);

% %th=degtorad([56.7,56.7,56.7,56.7, 51.3, 47.0, 43.5, 40.7, 38.3, 36.3, 34.5, 32.8, 31.3, 29.8, 28.3, 26.9, 25.5, 24.0, 22.4, 20.8]);
% % c=cm*[0, 0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.989, 0.994, 0.988, 0.961, 0.894, 0.732, 0.513];

%asitave chords....
c=[0
0
0
0.064815
  0.064815
  0.064815
  0.064815
  0.064815
  0.064815
  0.064815
  0.064815
  0.064815
  0.064815
  0.064137
  0.064427
  0.064008
  0.062297
  0.057939
  0.047481
  0.033279]';

% new chord for testing -- didn't use in the paper
% c=[0
%    0
%    0
%     0.0520
%     0.0585
%     0.0604
%     0.0607
%     0.0634
%     0.0649
%     0.0655
%     0.0654
%     0.0648
%     0.0636
%     0.0636
%     0.0638
%     0.0634
%     0.0617
%     0.0574
%     0.0470
%     0.0330]';
%% BEMT to be used to find the  / bound vortex on the disk

% Prandtl's tip-loss function
r(end)=0.99;r(1)=0.01;
for i=1:length(r)
    Fold=1;dF=1;
    while abs(dF)>.0001
        
        lam(i)=sqrt((sig*cla(i)/(16*Fold)-lamc*0.5)^2 + sig.*cla(i)*th(i)*r(i)*0.125/Fold)- (sig*cla(i)/(16*Fold)-lamc*0.5);
        phi(i)=atan(lam(i)/r(i));
        froot=0.5*nb*r(i)/((1-r(i))*phi(i));
        ftip=0.5*nb*(1-r(i))/(r(i)*phi(i));
     %   ftip=1000;
     %   froot=1000;
        
        %% based on pp 142 and 144 of Leishman -- does not work
        %  Fnew=(2/pi)*acos(exp(-froot*ftip))
        
        %% My propossal  
        Fnew=(2/pi)*acos(exp(-froot))*(2/pi)*acos(exp(-ftip));
        
        %%
        dF=abs(Fnew-Fold);
        Fold=Fnew;         
    end
   
    F(i)=Fnew;
end


% figure
% plot(r(4:end),F(4:end),'LineWidth',2);ylabel('Prandtl F function','FontSize',18);xlabel('r/R','FontSize',18);
% figure
%  plot(r(4:end),(th(4:end)- phi(4:end))*90/pi,'LineWidth',2);ylabel('\alpha_{eff}','FontSize',18);xlabel('r/R','FontSize',18);
figure;
plot(r(4:end),phi(4:end)*90/pi,'LineWidth',2);ylabel('Inflow angle (\phi)','FontSize',18);xlabel('r/R','FontSize',18);
cl=cla.*(th-phi);

% kutta-Joukowski theorem -- note we have to use the actual R which is r*R
% G=cl.*r*R.*c*om*0.5;

cl_asitav2=cl.*((r*R*om).^2+uinf^2).*c/340^2;
cl_asitav=cl_asitav2';
G=cl.*c.*sqrt((r*R*om).^2+uinf^2)*0.5;

vinf=lam*Vt;

% figure
% plot(r(4:end),cl(4:end));ylabel('C_l');xlabel('r/R');
% figure
% plot(r(4:end),G(4:end));ylabel('Circulation (G)');xlabel('r/R');

%plot(r,[G(1:end-1),0]);ylabel('G')
%plot(r,F)
% note that G is assumed to be zero at the tip, also the dG/dr is assumed
% to be zero at the root

gx=[0 nb*diff([G(1:end-1),0])/(2*pi*R)];
gt=[0 nb*om1/(2*pi)*diff([G(1:end-1),0])/uinf];
%gd=[G(1:end-1),0]./(2*pi);
gd=G;

%%
%figure(2)

for i=1:length(r)
    
    [X,Y,Z] = cylinder(i/length(r)*R*ones(sqrt(num2)),num1-1);
    % surf(X,Y,Z);
    x(:,:,i)=b*Z+xp;
    y(:,:,i)=Y+yp;
    z(:,:,i)=X+zp;
%     
%     surf(x(:,:,i),y(:,:,i),z(:,:,i));hold on
%     xlabel('x');ylabel('y');zlabel('z');
    size(X);
    
end                                                                                                                                 

end

