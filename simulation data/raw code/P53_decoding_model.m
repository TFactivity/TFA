function [finalResults]=P53_decoding_model(alpha_0,beta_0,gamma_0,n_0,kd_0,...
    A_basal_0,A_max_0,T_0,sigma_0,...
    tau_0,timelimit_0,timelag_record_0,cellnum_0,...
    sample_name,result_folder)

% This function is used for simulating unspliced RNA (u) and spliced RNA
% (s) of targets of oscillating TF (e.g., P53)
% Suppose that initial condition is the steady-state solution
% ODE:
    % du/dt=alpha*TF^n/(TF^n+kd^n)-beta*u
    % ds/dt=beta*u-gamma*s
    % TF(t)=A_max*(1-cos(2*pi*t/T))/2+A_basal*ги(1+cos(2*pi*t/T))/2)+N(0,sigma^2)
    % u0=(alpha*TF^n/(TF^n+kd^n))/beta and TF=A_basal
    % s0=(alpha*TF^n/(TF^n+kd^n))/gamma and TF=A_basal

mkdir(result_folder);


alpha=alpha_0;
beta=beta_0;
gamma=gamma_0; 
n=n_0;
kd=kd_0;

A_basal=A_basal_0;
A_max=A_max_0;
T=T_0;
sigma=sigma_0;

u0=floor((alpha*A_basal^n/(A_basal^n+kd^n))/beta);
s0=floor((alpha*A_basal^n/(A_basal^n+kd^n))/gamma);

tau=tau_0;
timelimit=timelimit_0;
timelag_record=timelag_record_0;

cellnum=cellnum_0;

tdur = 0:tau:timelimit;

for cellindex=1:cellnum
    
    %simulating TF dynamics
    TF_dynamics=A_max*(1-cos(2*pi*tdur/T))/2+A_basal*(1-(1-cos(2*pi*tdur/T))/2)+normrnd(0,sigma,1,length(tdur));
         
    x=zeros(1,4);
    x(1,1)=TF_dynamics(1);%initial value of TF dynamics
    x(1,2)=u0;%initial value of unspliced RNA
    x(1,3)=s0;%initial value of spliced RNA
    x(1,4)=0;%initial value of time
    
    counter=1;
    t=0;
    tcounter=0;
    
    tindex=1;
    %stochastic simulation of u and s based on Tau-leap method
    while t < timelimit
        
        if counter==1
            counter2=2;
        elseif counter==2
            counter2=1;
        end
        
        a(1)=alpha*(x(counter,1)^n)/(x(counter,1)^n+kd^n);
        a(2)=beta*x(counter,2);        
        
        a(3)=beta*x(counter,2);
        a(4)=gamma*x(counter,3);
               
        delta_x=[];
        for i=1:length(a)
            delta_x(i)=random('Poisson',a(i)*tau);
        end
        
        x(counter2,:)=x(counter,:);
        
        counter=counter2;
        x(counter2,1)=TF_dynamics(tindex);
        x(counter2,2)=max(x(counter,2)+delta_x(1)-delta_x(2),0);
        x(counter2,3)=max(x(counter,3)+delta_x(3)-delta_x(4),0);
        
        t=t+tau;
        x(counter2,4)=t;
        
        tindex=tindex+1;
        
        if (t-timelag_record*tcounter)>=0 %These are to reduce RAM consumption and speed up simulation.
            tcounter=tcounter+1;
            result(tcounter,:)=x(counter2,:); 
        end
        
    end
    
    time=result(:,4);
    TF=result(:,1);
    unspliced=result(:,2);
    spliced=result(:,3);
    
    finalResults{cellindex}.time=time;
    finalResults{cellindex}.TF=TF;
    finalResults{cellindex}.unspliced=unspliced;
    finalResults{cellindex}.spliced=spliced;
 
    finalResults{cellindex}.alpha=alpha;
    finalResults{cellindex}.beta=beta;
    finalResults{cellindex}.gamma=gamma;
    finalResults{cellindex}.n=n;
    finalResults{cellindex}.kd=kd;
    
    finalResults{cellindex}.A_basal=A_basal;
    finalResults{cellindex}.A_max=A_max;
    finalResults{cellindex}.T=T;
    finalResults{cellindex}.sigma=sigma;
    
    finalResults{cellindex}.u0=u0;
    finalResults{cellindex}.s0=s0;
    
    finalResults{cellindex}.tau=tau;
    finalResults{cellindex}.timelimit=timelimit;
    finalResults{cellindex}.cellnum=cellnum;    
    finalResults{cellindex}.tdur=tdur;
    
    save([result_folder,'\',sample_name,'_finalResults','.mat'],'finalResults');
end



