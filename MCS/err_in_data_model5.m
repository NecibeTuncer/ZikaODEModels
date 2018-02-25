function error_in_data5 = err_in_data_model5(k)

global tdata initial_cond  ydata  tpmeasure tforward ypdata


                 

[~,y] = ode45(@(t,y)model5_zika(y,k),tdata,initial_cond);
[~,yp] = ode45(@(t,y)model5_zika(y,k),tforward,initial_cond);
 CI = y(:,6)'; %cumulative incidences
 CIP = yp(tpmeasure(:),9)';

 error_in_data5 = sum((CI - ydata).^2) + sum((CIP - ypdata).^2) ;
 
 
end