clc
clear variables
warning('off','all')
rng(1)

K=10; % number of terminals
Marray = 40:20:60; % number of APs
nRuns =length(Marray);
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


nLargescaleruns=2;%large scale fading loop

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


MMselectaverage=zeros(1,nRuns);
maxIteration = 50;
errtol = 0.01;

for iM=1:nRuns
    M = Marray(iM);
    channelparams.nAPs = M;
    channelparams.nUsers = K;
    channelparams.pathloss = L;
    channelparams.dim = D;
    channelparams.shadowdev = sigma_shd;
    channelparams.refdist0 = d0;
    channelparams.refdist1 = d1;
    
    MMselect=M*K*ones(1,nLargescaleruns);
    
    %Power consumption parameters:
    myalpha=(1/0.4)*ones(M,1);
    P_fix=0;
    P_tc=0.2*ones(M,1);
    P_bt=0.25*10^(-3)*ones(M,1);
    P_0=0.825*ones(M,1);
    P_fix_bar=P_fix + N*sum(P_tc) + sum(P_0);
    
    for l=1:nLargescaleruns
        % Large-scale fading matrix
        mybeta=getslowfading(channelparams);
        
        
        % Create Gamma matrix
        den=zeros(M,K);
        for m=1:M
            for k=1:K
                den(m,k)=norm( (mybeta(m,:).^(1/2)).*(pilotseq(:,k)'*pilotseq))^2;
            end
        end
        
        mygamma=tau_p*rho_p*(mybeta.^2)./(tau_p*rho_p*den + 1);
        
        RateQoS=(tau_c/(tau_c-tau_p))*1*ones(K,1); % QoS requirement
        %Find the intial power control matrix
        [c_n,u_n,t_n] = generateinitialpoint(M,K,N,mygamma,mybeta,rho_d,pilotseq,RateQoS);
        if(~isnan(c_n)) % problem is feasible
            % Run Algorithm 1
            cdot = sdpvar(M,K,'full') ;
            tdot = sdpvar(K,1);
            udot = sdpvar(K,1);
            mytheta = sdpvar;
            cdot_n = c_n;
            obj= sum(tdot); % objective to be maximized; (36a); B is omitted
            opts=sdpsettings('solver','mosek','verbose',0,'dualize',0); % set internal solver to mosek
            convergence = 0;
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
                        interferencevectorvectorised(M,N,K,cdot,sqrt(rho_d)*mygamma,sqrt(rho_d)*mybeta,pilotseq,iUser);mytheta/N])]; %(36d)
                    
                    approx = approxfunctionvectorised(M,N,K,mygamma,mybeta,pilotseq,rho_d,cdot,udot,cdot_n,u_n,mytheta,iUser);
                    F = [F,cone([2*[sqrt(rho_d)*N*interferencevectorvectorised(M,N,K,cdot,mygamma,mybeta,pilotseq,iUser);mytheta]; ...
                        mytheta - approx],...
                        mytheta + approx)];% 36(g)
                    
                end
                diagnotics = optimize(F,-obj,opts);
                if (diagnotics.problem==0)
                    relincease=norm(u_n-value(udot/mytheta))/norm(u_n); % relative increase
                    u_n = value(udot/mytheta);
                    cdot_n = value(cdot/mytheta);
                else
                    disp('potential numerical issue, disregard the result, and move on to the next run')
                    break
                end
                if(relincease<errtol)
                    convergence = 1; % convergence is reached
                    break
                end
            end
            
            if(convergence)
                mygammaselect=mygamma;
                %thresh=0.1*1/M;
                c=value(cdot/mytheta);
                A=(c).*mygammaselect;
                PerMatrix=zeros(M,K);
                for k=1:K
                    for m=1:M
                        PerMatrix(m,k)=A(m,k)/sum(A(:,k));
                    end
                end
                for k=1:K
                    sumpercent = min(PerMatrix(:,k));
                    m_index=1;
                    [sortx,sortxindex]=sort(PerMatrix(:,k));
                    while sumpercent<0.05
                        mygammaselect(sortxindex(m_index),k)=0;
                        c(sortxindex(m_index),k)=0;
                        m_index=m_index+1;
                        sumpercent = sumpercent + sortx(m_index);
                        MMselect(l)=MMselect(l)-1;
                    end
                end
            end
        else
            disp('problem is not feasible')
            MMselect(l)=M*K;
        end
    end
    MMselectaverage(iM) = mean(MMselect)/K;
end
% 

plot(Marray,MMselectaverage,'r')

