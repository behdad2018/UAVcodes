% figure(1)
% plot(y,G);grid on;
% xlabel('wing span');ylabel('Vorticity distribution');
% figure(2)
% plot(y,cl_cL(5));grid on;
% xlabel('wing span');ylabel('C_l/C_L');
%  legend('\Lambda = 0','\Lambda = 10','\Lambda = 20','\Lambda = 30')
%  title('Effect of sweep angles')
%  print(gcf,'-djpeg',sprintf('-r%d',400),['sweep.jpg'])
% L_method(N,tr,b,AR,aoa,Lam,dih)
clc;clear all;
% L_method(N,tr,b,AR,aoa,Lam,dih)
[y,cL,cl_cL,G]=L_method(6,1,1,8,7,0,0);

% for i=1:16
% aoa=i-1
% [y,cL(i),cl_cL]=L_method(40,1,1.8,6,aoa,0,0);
% 
% end
% 
% plot([0:15],cL)
% grid on;

% xlabel('\alpha');ylabel('c_L')
% 
% print(gcf,'-djpeg',sprintf('-r%d',400),['cL_for_Lawrens_wing.jpg'])
