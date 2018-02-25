function  total_ARE = MonteCarlo_ZikaEpiModel2

clear all
close all
clc

numiter = 1000; 

conv = zeros(6,numiter);

X = zeros(4,numiter);

Lambda_v = 1/10;
mu_v = 1/10;

mu = 1/(79*365);
Lambda = 19380000*mu;
omega = 0.0;
phi=0.2;
q=1;


tdata = [1 3 9 14 15 18 21 22 23 24 25 30 31 32 36 37 39 42 43 44 50 56 57 58 59 60 63 64 66 67 71 72 73];

tforward = (1:0.1:73);

tmeasure = [1 21 81 131 141 171 201 211 221 231 241 291  301 311 351 361 381 411 421 431 491 551 561 571 581 591 ...
            621 631 651 661 701 711 721];


true_params =1.0e+02 *[5.475785680175260   5.193474972033285   0.001000000002151   0.001000000000086];
   


initial_cond = [0.99 0.01 19380000 3 0 1  10];

[~,y_trp] = ode45(@(t,y)(model_zikav1(y,true_params)),tforward,initial_cond);





function dy = model_zikav1(y,k)

dy = zeros(7,1);

beta = k(1);
beta_v = k(2);
gamma = k(3);
gammaA = k(4);


Sv = y(1);
Iv = y(2);
S = y(3);
I = y(4);
R = y(5);
C = y(6);

A = y(7);

dy(1) = Lambda_v - beta* (I+q*A)* Sv/(S+I+A+R) - mu_v*Sv;
dy(2) = beta*(I+q*A)*Sv/(S+I+A+R) - mu_v*Iv;
dy(3) = Lambda - beta_v* Iv*S/(S+I+A+R)  - mu*S + omega*R;
dy(7) = (1-phi)*beta_v*Iv*S/(S+I+A+R)  - (mu + gammaA)*A;
dy(4) = phi* beta_v*Iv*S/(S+I+A+R)  - (mu + gamma)*I;
dy(5) = gamma*I + gammaA * A - (mu + omega)*R;
dy(6) = phi* beta_v* Iv*S/(S+I+A+R);
end

function error_in_data = err_in_data(k)


 [t,y] = ode45(@(t,y)model_zikav1(y,k),tdata,initial_cond);
 
 CI = y(:,6)'; %cumulative incidences
  
 error_in_data = sum((CI - ydata).^2) ;
 
 

end


noiselevel = [0, 0.01, 0.05, 0.1, 0.2, 0.3];
total_ARE =  zeros(length(noiselevel), length(true_params));

for noisei = 1:6
    
rng default
noiselev = noiselevel(noisei)


for i= 1:numiter
    i
    ytrp = y_trp(tmeasure(:),6);
  ydata =  (noiselev*(ytrp').*randn(1,size(tmeasure,2))) + ytrp';


k = true_params;

lb = [0      0      0.1   0.1]; 
ub = [10000  10000  0.5   0.5];
 
[k,fval] =  fminsearchbnd(@err_in_data,k,lb,ub,optimset('Display','iter','MaxIter',2000,'MaxFunEvals',2000)) 
                           
X(:,i) = k';
conv(noisei,i) = exitflag;
 

plot(tdata,ydata,'o',tforward,y_trp(:,6),'LineWidth',1.5)
title('Cumulative New Cases Generated Data')

end



arescore = zeros(1,length(true_params));
    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = arescore;



end
end