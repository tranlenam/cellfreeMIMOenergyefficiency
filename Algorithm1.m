clc
clear variables
warning('off','all')
rng('default')

K=4; % number of terminals
M=20; % number of APs
N=2; % number of antennas/AP
B=20; % bandwidth in Mhz

tau_c=200; % coherence time (in symbols)
tau_p=20; % length of pilot sequences (in symbols)
D=1; %in kilometer.
[U,~,~]=svd(randn(tau_p,tau_p));% U includes tau_p orthogonal sequences 

Hb = 15; % Base station height in m
Hm = 1.65; % Mobile height in m
f = 1900; % Frequency in MHz

% path loss parameters
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

power_f=N*1; %downlink power: 1W

noise_figure = 9; % noise figure
noise_p = 10^((-203.975+10*log10(B*10^6)+noise_figure)/10); %noise power

rho_d = power_f/noise_p; % nomalized tx power
rho_p= 0.2/noise_p; % nomalized pilot power

sigma_shd=8; % standard deviation with shadowing, in dB

d0=0.01;% km
d1=0.05;% km

channelparams.nAPs = M;
channelparams.nUsers = K;
channelparams.pathloss = L;
channelparams.dim = D;
channelparams.shadowdev = sigma_shd;
channelparams.refdist0 = d0;
channelparams.refdist1 = d1;


%Power consumption parameters:
myalpha=(1/0.4)*ones(M,1);
P_fix=0;
P_tc=0.2*ones(M,1);
P_bt=0.25*10^(-3)*ones(M,1);
P_0=0.825*ones(M,1);
P_fix_bar=P_fix + N*sum(P_tc) + sum(P_0);

% Generate large-scale fading matrix
mybeta=getslowfading(channelparams);

% Pilot Asignment: (random choice)
pilotseq=zeros(tau_p,K); % pilot sequences, the length of each sequence is tau_p
if tau_p<K
    pilotseq(:,1:tau_p)=U;
    for iUser=(tau_p+1):K
        pilotseq(:,iUser)=U(:,randi([1,tau_p]));
    end
else
    pilotseq=U(:,1:K);
end

% Create gamma matrix defined in (5)
den=zeros(M,K);
for iAP=1:M
    for iUser=1:K
        den(iAP,iUser)=norm((mybeta(iAP,:).^(1/2)).*(pilotseq(:,iUser)'*pilotseq))^2;
    end
end

mygamma=tau_p*rho_p*(mybeta.^2)./(tau_p*rho_p*den + 1);


maxIteration = 30;
EE_max_sub=zeros(maxIteration,1);
EE_max=zeros(maxIteration,1);
epsi = 0.1;

RateQoS = (tau_c/(tau_c-tau_p))*1*ones(K,1); % minimum spectral efficiency
[c_n,u_n,t_n] = generateinitialpoint(M,K,N,mygamma,mybeta,rho_d,pilotseq,RateQoS);
if isnan(c_n) % problem is infeasible, stop. If required, reduce RateQoS to make it feasible
    return
else
    cdot_n = c_n;
    %% define optimization variables
    cdot = sdpvar(M,K,'full') ;
    tdot = sdpvar(K,1);
    udot = sdpvar(K,1);
    mytheta = sdpvar;
    
    obj= sum(tdot); % objective to be maximized; (36a); B is omitted
    opts=sdpsettings('solver','mosek','verbose',0,'dualize',0); % set internal solver to mosek
    %% main SCA loop
    for iIter=1:maxIteration
        F=[]; % reset the constraints to emmpty set
        F = [F,tdot(:)>=0];
        F = [F,cone([mytheta/sqrt(N)*ones(1,M);(sqrt(mygamma).*cdot)'])]; %(36b) 
        F = [F,cdot(:)>=0]; %  (36c)
        Gammaan_temp = sqrt(rho_d*noise_p*N*(repmat(myalpha,1,K).*mygamma));
        F = [F,cone([sqrt(P_fix_bar)*mytheta;Gammaan_temp(:).*cdot(:);0.5*(mytheta-1)],0.5*(mytheta+1))]; % (36e)
        F = [F,cone([(udot+mytheta*(log(u_n)+1)-log(2)*tdot)';(udot-mytheta*(log(u_n)+1)+log(2)*tdot)';...
            2*sqrt(u_n)'*mytheta])]; %(36f)
        
        for iUser=1:K
            F = [F,cone([(1/sqrt(2^( RateQoS(iUser)) - 1))*cdot(:,iUser)'*(sqrt(rho_d)*mygamma(:,iUser));...
                interferencevector(M,N,K,cdot,sqrt(rho_d)*mygamma,sqrt(rho_d)*mybeta,pilotseq,iUser);mytheta/N])]; %(36d)
            
            approx = approxfunction(M,N,K,mygamma,mybeta,pilotseq,rho_d,cdot,udot,cdot_n,u_n,mytheta,iUser);
            F = [F,cone([2*[sqrt(rho_d)*N*interferencevector(M,N,K,cdot,mygamma,mybeta,pilotseq,iUser);mytheta]; ...
                mytheta - approx],...
                mytheta + approx)];% 36(g)
            
        end
        diagnotics = optimize(F,-obj,opts);
        if (diagnotics.problem==0)
            u_n = value(udot/mytheta);
            cdot_n = value(cdot/mytheta);
        else
            disp('potential numerical issue, disregard the result, and move on to the next run')
            break
        end
        
        %% Energy efficiency
        EE_max_sub(iIter) = B*(1-tau_p/tau_c)*sum(double(tdot)); % (36a) 
        EE_max(iIter)=  1/(1/EE_max_sub(iIter) + sum(P_bt)); % the sequence of energy efficiency objective
        
    end
    %plotting
    plot(1:length(EE_max),EE_max,'r')
    xlabel('Iteration count')
    ylabel('Energy Efficiency')
end

