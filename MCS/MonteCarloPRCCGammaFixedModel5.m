clear all
close all
clc

global tdata ydata tpdata ypdata tforward tpmeasure initial_cond Lambda_v mu_v Lambda mu omega gamma gammap xi

format long

numiter = 2000;%1000;
 
Lambda_v = 1/10;
mu_v = 1/10;

mu = 1/(79*365);
Lambda = 19380000*mu;
omega = 0.0;
xi = 0.0000007;



load ZikaLocal2.txt
tdata = ZikaLocal2(:,1);
ydata = ZikaLocal2(:,2)';

tpdata = [ 3 30 32 43 44 64];
ypdata = [ 1 2 3 4 5 6];

tforward = (1:0.1:73);
tpmeasure = [21 291 311 421 431 631];

initial_cond = [0.99 0.01 19380000 3 0 1 387600 1 1];


lb = [0  0  0 0 0];



% true_params = 1.0e+03 *[2.551547085324936   0.093007119888124...
%                         0.000000000000033...
%                         0.385140893575393   0.000000000219522];

true_params = 1.0e+03 *[2.617716723300200   0.091377567889342...
                        0.000000000006291   0.398991459973893   0.000000024919599];
                    
true_reproduction_number = 1.74;                     

Reproduction_Number =  zeros(1, numiter);
betavfitted =  zeros(1, numiter);
betafitted =  zeros(1, numiter);
betadfitted =  zeros(1, numiter);
betavpfitted =  zeros(1, numiter);
betadpfitted =  zeros(1, numiter);

k = true_params;

 for i= 1:numiter
    i
    
    gamma = normrnd(0.1,0.1);
    gammap = normrnd(0.02,0.005);
    
if gamma > 0.07 
  
k =  fminsearchbnd(@err_in_data_model5,k,lb,[]);%,optimset('Display','iter','MaxIter',200,'MaxFunEvals',200)) 

 
betafitted(i) = k(1);
betavfitted(i) = k(2);
betadfitted(i) = k(3);
betavpfitted(i) = k(4);
betadpfitted(i) = k(5);

S0 = mu/(mu+xi);
Sp0 = xi/(mu+xi);
M = Lambda_v*mu/(mu_v*Lambda);

r1 = k(3)*S0/(mu+gamma);
r2 = k(5)*Sp0/(mu+gammap) * (k(3)*S0/(mu+gamma) + k(2)*S0*k(1)*M/(mu_v*(mu+gamma)));
r3 = k(2)*S0*k(1)*M/(mu_v*(mu+gamma));
r4 = k(4)*Sp0*k(1)*M/(mu_v*(mu+gammap));
Reproduction_Number(i) = r1+r2+r3+r4;



  
end
 end
 
 
 
 betafitted = betafitted(betafitted ~=0);
 betavfitted = betavfitted(betavfitted ~=0);
 betadfitted = betadfitted(betadfitted ~=0);
 betavpfitted = betavpfitted(betavpfitted ~=0);
 betadpfitted = betadpfitted(betadpfitted ~=0);
 
 Reproduction_Number = Reproduction_Number(Reproduction_Number~=0);
%  
%  histogram(Reproduction_Number)
%  hold on 
%  plot(true_reproduction_number, 0,'pr','MarkerSize',20,'MarkerFaceColor','r')
%  figure
%  histogram(betafitted)
%  hold on 
%  plot(true_params(1), 0,'pr','MarkerSize',20,'MarkerFaceColor','r')
%  figure
%  histogram(betavfitted)
%  hold on 
%  plot(true_params(2), 0,'pr','MarkerSize',20,'MarkerFaceColor','r')
%  
%  
 nl = length(betafitted);
 ares_R0 = 100*sum(abs(true_reproduction_number - Reproduction_Number)/abs(true_reproduction_number))/nl;
 ares_beta = 100*sum(abs(true_params(1) - betafitted)/abs(true_params(1)))/nl;
 ares_beta_v = 100*sum(abs(true_params(2) - betavfitted)/abs(true_params(2)))/nl;
 ares_beta_d = 100*sum(abs(true_params(3) - betadfitted)/abs(true_params(3)))/nl;
 ares_beta_vp = 100*sum(abs(true_params(4) - betavpfitted)/abs(true_params(4)))/nl;
 ares_beta_dp = 100*sum(abs(true_params(5) - betadpfitted)/abs(true_params(5)))/nl;
 
 save('Model5Results')