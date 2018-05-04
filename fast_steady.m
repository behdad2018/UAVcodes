% simplified unsteady 3D wing lifting line method by B. Davoudi
% Aerospace Engineering Department, University of Michigan 7/27/2016
function [G,Am_wing,A,a,a_d,w_ind_drag]=fast_steady(x,y,z,xcol,zcol,n,dl_x,dly,Nxw,Nxt,Nyw,Nyt,u,alpha,Lam,dih,b)
% genrating the wake for a steady flight
% finding the steady solution with the panel method

% influence matrix is the summation of the wing and wake contributions:

Awing=zeros(2*Nxw*Nyw+2*Nxt*Nyt,2*Nxw*Nyw+2*Nxt*Nyt);
Awing_drag=zeros(2*Nxw*Nyw+2*Nxt*Nyt,2*Nxw*Nyw+2*Nxt*Nyt);
Awake=zeros(2*Nxw*Nyw+2*Nxt*Nyt,2*Nxw*Nyw+2*Nxt*Nyt);

% this is only the contribution of the wakes on each panel, velocity in the three direction 

a1=zeros(2*Nxw*Nyw+2*Nxt*Nyt,2*Nxw*Nyw+2*Nxt*Nyt,3);
a1_d=zeros(2*Nxw*Nyw+2*Nxt*Nyt,2*Nxw*Nyw+2*Nxt*Nyt,3);
%for i=1:2*Nx*Ny

for i=1:2*Nxw*Nyw+2*Nxt*Nyt  % total of the panels in the wing and wake 
    
    % once the panels on the wing are considered, tail panels are the next
    % starting from 2*Nxw*Nyw+1 to 2*Nxw*Nyw+2*Nxt*Nyt
    
    % for mIN wing
    index1col=ceil(i/(2*Nyw));
    index2col=i-(2*Nyw)*(index1col-1);
    
    % once we start consdiering the tail panels
    if index1col>Nxw
        
        index1col=Nxw+ceil((i-2*Nxw*Nyw)/(2*Nyt));        % goes from Nxw+1 to Nxt
        index2col=i-2*Nxw*Nyw-(2*Nyt)*(index1col-Nxw-1);  % goes from 1 to Nyt
        
    end
    
    % collocation points for both wing and tail
    xcol1=xcol(index1col,index2col);
    ycol1=y(index1col,index2col);
    zcol1=zcol(index1col,index2col);
    %% for finding matrix A only panels on the wing and tail
    for j=1:2*Nxw*Nyw+2*Nxt*Nyt
        
        index1=ceil(j/(2*Nyw));
        index2=j-(2*Nyw)*(index1-1);
        index3=j-(2*Nyw)*(index1-1);
        
        if index1>Nxw
            index1=Nxw+ceil((j-2*Nxw*Nyw)/(2*Nyt));
            index2=j-2*Nxw*Nyw-(2*Nyt)*(index1-Nxw-1);
            index3=index2+2*Nyw;
        end
        x1=x(index1,index2);
        y1=y(index1,index2);
        z1=z(index1,index2);
        % y=1/2*x+3/2 provides 1 and 2 indices to be used for dly
        [~,Awing(i,j)]=vortexring(n(:,i),dl_x(index3),dly(0.5*sign(j-2*Nxw*Nyw-0.01)+1.5),alpha(index1,index2),Lam,dih,xcol1,ycol1,zcol1,x1,y1,z1,1,0);       % Effect of ring vortices on the wing
        % for drag calculation
        [~,Awing_drag(i,j)]=vortexring_drag(n(:,i),dl_x(index3),dly(0.5*sign(j-2*Nxw*Nyw-0.01)+1.5),alpha(index1,index2),Lam,dih,xcol1,ycol1,zcol1,x1,y1,z1,1,0);       % Effect of ring vortices on the wing
    end                
    %% wake effect to be included in the A matrix
    % Note: we have two wakes one from the wing and one from the tail. Also
    % we think of the wake as if it were extended of the last panel
    % wing:
    for j=1:2*Nyw
        
        % first row of the wake, index1=Nxw
        xw1=x(Nxw,j)+dl_x(j)*cos(alpha(Nxw,j));
        yw1=y(Nxw,j);
        zw1=z(Nxw,j)-dl_x(j)*sin(alpha(Nxw,j))*cos(dih);
        
        [a1(i,2*(Nxw-1)*Nyw+j,:),Awake(i,2*(Nxw-1)*Nyw+j)]=vortexring(n(:,i),20*b,dly(1),0,Lam,dih,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,0);
        [a1_d(i,2*(Nxw-1)*Nyw+j,:),~]=vortexring_drag(n(:,i),20*b,dly(1),0,Lam,dih,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,0);
    end
    % tail:
    for j=1:2*Nyt
        
        % first row of the wake index1=Nxt
        xw1=x(Nxt+Nxw,j)+dl_x(2*Nyw+j)*cos(alpha(Nxt+Nxw,j));
        yw1=y(Nxt+Nxw,j);
        zw1=z(Nxt+Nxw,j)-dl_x(2*Nyw+j)*sin(alpha(Nxt+Nxw,j))*cos(dih);
        [a1(i,2*Nxw*Nyw+2*(Nxt-1)*Nyt+j,:),Awake(i,2*Nxw*Nyw+2*(Nxt-1)*Nyt+j)]=vortexring(n(:,i),20*b,dly(2),0,Lam,dih,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,0);
        [a1_d(i,2*Nxw*Nyw+2*(Nxt-1)*Nyt+j,:),~]=vortexring_drag(n(:,i),20*b,dly(2),0,Lam,dih,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,0);
   
    end
end

A=Awing+Awake;

% the inverse matrix A that will be used in the unsteady solution

Am_wing=inv(Awing);

U=repmat(u',1,2*Nxw*Nyw+2*Nxt*Nyt);
%temp
% U(1,Nyw)=0;U(1,Nyw+1)=0;
% U(1,3*Nyw)=0;U(1,3*Nyw+1)=0;
% U(1,5*Nyw)=0;U(1,5*Nyw+1)=0;
% U(1,7*Nyw)=0;U(1,7*Nyw+1)=0;
%
%% propeller % momentum thoery -- not a good thing
%
% slope=2*pi;
% R=b/10;
% cprop=b/70;
% sig=4*(cprop)/(pi*R);
% om=3000*2*pi/60*b/20;
% th=0.2;
% v1=1/16*slope*sig*om*R*(sqrt(1+24*th/(slope*sig))-1);
%
% for i=1:Nxw
%     U(1,Nyw+(i-1)*2*Nyw)=u(1)+2*v1;
%     U(1,Nyw+1+(i-1)*2*Nyw)=u(1)+2*v1;
%     % tale
% end
%
% for i=1:Nxt
%     U(1,2*Nxw*Nyw+Nyt+(i-1)*2*Nyt)=u(1)+v1;
%     U(1,2*Nxw*Nyw+Nyt+1+(i-1)*2*Nyt)=u(1)+v1;
% end
%%
RHS=dot(-U,n)';


%RHS=dot(-repmat(u',1,2*Nxw*Nyw+2*Nxt*Nyt),n)';

G=A\RHS;    % vorticity on the wing
size(G);

% these are the induced velocty due to the wake on the wing panels
% note that a1(:,:,i) is a big influence matric that only has non zero
% values where the wake is talking to a panel on thw wing... and when
% panels are talking to each other it has zero values. When mutiply this 
% by G, we get the actual indiced velicty done by the wake on a panel

a(1,:)=a1(:,:,1)*G;
a(2,:)=a1(:,:,2)*G;
a(3,:)=a1(:,:,3)*G;

% its not clear that, in the wake, ww should eliminate the span-wise vortex
% or not --- here to include it:

 a_d(1,:)=a1_d(:,:,1)*G;
 a_d(2,:)=a1_d(:,:,2)*G;
 a_d(3,:)=a1_d(:,:,3)*G;

w_ind_drag=Awing_drag*G;

end



