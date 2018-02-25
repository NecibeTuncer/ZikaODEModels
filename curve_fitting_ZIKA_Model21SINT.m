function curve_fitting_ZIKA_Model21SINT

clear all
close all
clc


% load ZikaFlorida.txt
% tdata = ZikaFlorida(:,1);
% ydata = ZikaFlorida(:,2);

% load ZikaLocal.txt
% tdata = ZikaLocal(:,1);
% ydata = ZikaLocal(:,2);

load ZikaLocal2.txt
tdata = ZikaLocal2(:,1);
ydata = ZikaLocal2(:,2);

tforward = [tdata(1):0.1:tdata(end)+2]';

format long

Lambda_v = 1/10;
mu_v = 1/10;

mu = 1/(79*365);
Lambda = 19380000*mu;
omega = 0.0;
% Lambda = 20000000*mu;



% 
% k = 1.0e+03 * [0.000000028227958   0.000000005328385   
%                0.000026657369222   8.327830807908896] %res 310, best fit so far
           
% k = 1.0e+03 * [0.000000028227958     0.000000005328385   
%                0.00026657369222    0.00008327830807908896]
% 
% k = [0.000015472837978   0.000005456979318...
%      0.14                0.0013298604935864];

% k = [0.000047753182369   0.000004812724812   0.060000126691383   0.010000826532716]; %res 3.520741420749052e+02
%k  =[ 0.000029750989454   0.000008481267878   0.010359716858378   0.081662056817416]; % res 2.991445365967680e+02

% k=[0.000145767388383   0.000005180243571   0.100000000000010   0.185071294814277]; %res 3.71e+02  
%k = [0.000047753182369   0.000004812724812   0.060000126691383];
% k = [0.145767388383   1.5180243571   0.100000000000010];
k= 1.0e+03 *[2.824984510512934   0.100392926290024   0.000100000000016]; %SI 
% initial_cond = [10^6 10^2 19380000 3 0 1];

initial_cond = [0.99 0.01 19380000 3 0 1];




function dy = model_zikav1(y,k)

dy = zeros(6,1);



beta = k(1);
beta_v = k(2);
gamma = k(3);


Sv = y(1);
Iv = y(2);
S = y(3);
I = y(4);
R = y(5);
C = y(6);

dy(1) = Lambda_v - beta* I* Sv/(S+I+R) - mu_v*Sv;
dy(2) = beta*I*Sv/(S+I+R) - mu_v*Iv;
dy(3) = Lambda - beta_v* Iv*S/(S+I+R) - mu*S + omega*R;
dy(4) = beta_v*Iv*S/(S+I+R) - (mu + gamma)*I;
dy(5) = gamma*I - (mu + omega)*R;
dy(6) = beta_v* Iv*S/(S+I+R);
end

function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
 
 CI = y(:,6); %cumulative incidences
  
 %Incidences = k(2)*y(:,2)*y(:,3); 
 error_in_data = sum((CI - ydata).^2) ;
 
 

end


lb = [0      0     1/10  ];
ub = [10000  10000   0.5 ];

%for j=1:2
%[k,fval] =  fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter','MaxIter',2000,'MaxFunEvals',2000)) 
% k = fmincon(@err_in_data,k,[],[],[],[],lb,ub,[],...
%     optimoptions('fmincon','Display','iter'))
%end


r1 = (k(1)*k(2)*Lambda_v*mu);
r2 = (mu_v^2 * (mu + k(3)) *Lambda);
Reproduction_Number = r1/r2

 figure(1)

 [t y] = ode45(@(t,y)model_zikav1(y,k),tforward,initial_cond);
 plot(tdata, ydata, 'o', tforward, y(:,6), '-r')
 
 figure(2)

 tprojection = 0:0.1:2000;
%  [t y] = ode45(@(t,y)model_zikav1(y,k),tprojection,initial_cond);
%  plot(tdata, ydata, 'o', tprojection, y(:,6), '-r')
 
 

   [t y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
residuals = (ydata - y(:,6))./y(:,6);
mn = mean(residuals)*ones(1,length(ydata));
  plot(tdata, residuals, 'or', tdata, mn,'-b')
  
 figure(3)

 [t y] = ode45(@(t,y)model_zikav1(y,k),tprojection,initial_cond);
 plot(tprojection, y(:,4), '-r')
% % 
 end