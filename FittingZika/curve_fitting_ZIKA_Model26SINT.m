function curve_fitting_ZIKA_Model26SINT

clear all
close all
clc


load ZikaLocal2.txt
tdata = ZikaLocal2(:,1);
ydata = ZikaLocal2(:,2);


tpdata = [ 3 30 32 43 44 64];
ypdata = [ 1 2 3 4 5 6]';

tforward = [tdata(1):0.1:tdata(end)+2]';
tpmeasure = [21 291 311 421 431 631];
format long

Lambda_v = 1/10;
mu_v = 1/10;

mu = 1/(79*365);
Lambda = 19380000*mu;
omega = 0.0;
xi = 0.0000007;


%    k = 1.0e+03 *[0.331799685664064   0.053624425966729   0.000100000000004   0.000041864211775...
%                       3.298223121567662   0.000000000004620   0.000020000000003];
 
   k =1.0e+03 *[2.551547085324936   0.093007119888124...
               0.000100000008075...
               0.385140893575393...
               0.000020000001438]
         
initial_cond = [0.99 0.01 19380000 3 0 1 387600 1 1];




function dy = model_zikav1(y,k)

dy = zeros(8,1);



beta = k(1);
beta_v = k(2);
gamma = k(3);
beta_vp = k(4);
gammap = k(5);


Sv = y(1);
Iv = y(2);
S = y(3);
I = y(4);
R = y(5);
C = y(6);

Sp = y(7);
Ip = y(8);
Cp = y(9);

dy(1) = Lambda_v - beta* (I+ Ip)* Sv /(S+I+R+Sp +Ip) - mu_v*Sv;
dy(2) = beta*(I+Ip)*Sv/(S+I+R+Sp +Ip) - mu_v*Iv;
dy(3) = Lambda - beta_v* Iv*S/(S+I+R+Sp +Ip) - (mu+xi)*S + omega*R + (1/280)*Sp;
dy(7) = xi*S - beta_vp * Iv * Sp/(S+I+R+Sp +Ip)  - mu*Sp - (1/280)*Sp;
dy(4) = beta_v*Iv*S/(S+I+R+Sp +Ip)  - (mu + gamma)*I;
dy(8) = beta_vp * Iv * Sp/(S+I+R+Sp +Ip)  -(mu+gammap) *Ip;
dy(5) = gamma*I +gammap * Ip - (mu + omega)*R;
dy(6) = beta_v* Iv*S/(S+I+R+Sp+Ip) +beta_vp * Iv * Sp/(S+I+R+Sp+Ip) ;
dy(9) = beta_vp * Iv * Sp/(S+I+R+Sp+Ip);
end

function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
 [t,yp] = ode45(@(t,y)model_zikav1(y,k),tforward,initial_cond);
 CI = y(:,6); %cumulative incidences
 CIP = yp(tpmeasure(:),9); 
 %Incidences = k(2)*y(:,2)*y(:,3); 
 error_in_data = sum((CI - ydata).^2) + sum((CIP - ypdata).^2) ;
 
 

end


lb = [0        0       1/10     0            1/50 ];
%ub = [1000000  100000  0.5    1000000  1/7 ];

 for j=1:2
% 
% k = fmincon(@err_in_data,k,[],[],[],[],lb,ub,[],...
%     optimoptions('fmincon','Display','iter'))
% end
[k,fval] =  fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter','MaxIter',2000,'MaxFunEvals',2000)) 
 end
S0 = mu/(mu+xi);
Sp0 = xi/(mu+xi);
M = Lambda_v*mu/(mu_v*Lambda);

% r1 = k(4)*S0/(mu+k(3));
% r2 = k(6)*Sp0/(mu+k(7)) * (k(4)*S0/(mu+k(3)) + k(2)*S0*k(1)*M/(mu_v*(mu+k(3))));
% r3 = k(2)*S0*k(1)*M/(mu_v*(mu+k(3)));
% r4 = k(5)*Sp0*k(1)*M/(mu_v*(mu+k(7)));
% Reproduction_Number = r1+r2+r3+r4

 figure(1)

 [t y] = ode45(@(t,y)model_zikav1(y,k),tforward,initial_cond);
 plot(tdata, ydata, 'o', tforward, y(:,6), '-r')
 
  figure(2)
  plot(tpdata, ypdata, 'o', tforward, y(:,9), '-r')
% 
%  tprojection = 0:0.1:2000;
% %  [t y] = ode45(@(t,y)model_zikav1(y,k),tprojection,initial_cond);
% %  plot(tdata, ydata, 'o', tprojection, y(:,6), '-r')
%  
%  
% 
%    [t y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
% residuals = (ydata - y(:,6))./y(:,6);
% mn = mean(residuals)*ones(1,length(ydata));
%   plot(tdata, residuals, 'or', tdata, mn,'-b')
%   
%  figure(3)
% 
%  [t y] = ode45(@(t,y)model_zikav1(y,k),tprojection,initial_cond);
%  plot(tprojection, y(:,4), '-r')
% % % 

%The results

k =1.0e+03 *[2.617723435993325   0.091377380457317   0.000100000000001   0.399007853530409...
             0.000020000000003];

 end