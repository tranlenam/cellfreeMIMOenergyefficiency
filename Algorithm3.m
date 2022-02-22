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
            mygammaselect=mygamma;
            PerMatrix=zeros(M,K);
            for k=1:K
                for m=1:M
                    PerMatrix(m,k)=mygamma(m,k)/sum(mygamma(:,k));
                end
            end
            for k=1:K
                sumpercent = min(PerMatrix(:,k));
                m_index=1;
                [sortx,sortxindex]=sort(PerMatrix(:,k));
                while sumpercent<0.05
                    mygammaselect(sortxindex(m_index),k)=0;
                    m_index=m_index+1;
                    sumpercent = sumpercent + sortx(m_index);
                    MMselect(l)=MMselect(l)-1;
                end
            end
            % Run Algorithm 1 to find optimal power allocation, omit here
            % to simplicity
 
        else
            disp('problem is not feasible')
            MMselect(l)=M*K;
        end
    end
    MMselectaverage(iM) = mean(MMselect)/K;
end
%

plot(Marray,MMselectaverage,'r')

