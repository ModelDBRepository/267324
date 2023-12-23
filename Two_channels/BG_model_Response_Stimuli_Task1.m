function [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,k_tap_vett,Uchi,ChI,t] = BG_model_Response_Stimuli_Task1(S,Wgc,Wgs,Wnc,Wns,STN_ON,T_ON,Dop_tonic,noiseC)

%% function used to simulate the task
tau = 15;   %basal time constant
tauL = 5*tau;   %time constant of lateral inhibition
%tauS = tau/5;

dt = 0.1;   %step
t = (0:dt:800)';   %time [ms] %800
D = length(t);   %number of samples


% %--------------------------------------------------------------------------
% % Initialisation structures
% %--------------------------------------------------------------------------

%C: cortex
Nc = 2;   %neurons in the cortex
C = zeros(Nc,D);   %activity of neurons
Uc = zeros(Nc,D);   %input to the sigmoid
Ul = zeros(Nc,D);   %input from the lateral ihibition
%Go
Go = zeros(Nc,D);
Ugo = zeros(Nc,D);
%NoGo
NoGo = zeros(Nc,D);
Unogo = zeros(Nc,D);
%Gpe: globus pallidus pars externa
Gpe = zeros(Nc,D);
Ugpe = zeros(Nc,D);
%Gpi: globus pallidus pars interna
Gpi = zeros(Nc,D);
Ugpi = zeros(Nc,D);
%T: thalamus
T = zeros(Nc,D);
Ut = zeros(Nc,D);
%STN: sub thalamic nucleus
STN = zeros(1,D);
Ustn = zeros(1,D);
%E: energy (as an index of cortical conflict)
E = zeros(1,D);
%ChI: cholinergic interneurons
ChI = zeros(1,D);
Uchi = zeros(1,D);

%input from DA+ACh to Go and NoGo
IGo_DA_Ach = zeros(Nc,D);
INoGo_DA_Ach = zeros(Nc,D);


%% initialization of synapses

%weights from stimulus to cortex
Wcs = 1*ones(Nc,Nc);   %0.2 extradiag, 1.1 su diag
%lateral inhibition
L0 = 1.2;
L = -L0*ones(Nc,Nc)+diag(L0*ones(Nc,1));   %tolto autoanello 1.2
%weights from thalamus to cortex
Wct = 4*diag(ones(Nc,1));  %diagonal

% %%%%%%%%%%%% trained synapses are passed as external parameters %%%%%%%%%%%%%%%%%%%%%
% 

% Wgc = 
% Wgs =
% Wnc = 
% Wns =
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%weights from No-Go a Gpe (inhibitory))
Wen = -2*diag(ones(Nc,1));   %diagonal
Wen = Wen+10/100*Wen;   %diagonal

%weights from Gpe to Gpi (inhibitory)
Wie = -3*diag(ones(Nc,1));   %diagonal
%weights from  Go to Gpi (inhibitory)
Wig = -36*diag(ones(Nc,1));   %diagonal

%weights from cortex to thalamus (excitatory)
Wtc = 3*diag(ones(Nc,1));   %diagonal
%weights from Gpi to thalamus (inhibitory)
Wti = -3*diag(ones(Nc,1));   %diagonal


%%%%%%%%%%%%%%%%%%% weights from - to STN %%%%%%%%%%%%%%%%%%%
%weight from energy to STN (excitatory)
Ke = 7;   % now it is a parameter 7
%weight from  Gpe to STN (inhibitory)
Kgpe = -1;

%weight from  STN to Gpe (excitatory)
Westn = 1;

%weight from  STN to Gpi (excitatory)
Wistn = 30;   %14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% weight from ChI to Go-NoGo %%%%%%%%%%%%%%%%%%%
%weight from ChI to Go (inhibitory)
wgchi = -1;

%weight from ChI to NoGo (excitatory)
wnchi = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% guadagni %%%%%%%%%%%%%%%%%%%
%gain from DA to Go (excitation)
Ugo_trigger = 1.074;

alpha   =0.75 ; 
beta    =-1; 
gamma   =-0.5 ;
%% calcolo attività dinamica delle strutture: dinamica (passabasso) + sigmoide

%parameters of the sigmoid
a = 4;
U0 = 1.0;

%tonic activity of the Gpe
Igpe = 1.0;
%tonic activity of the Gpi
Igpi = 3;

%tonic activity of the Chi
Ichi = 1.00;

% tonic dopamine
DA = Dop_tonic;



%% initial conditions


Ugpe(:,1) = Igpe;
Ugpi(:,1) = Igpi;
Uchi(1) = Ichi+gamma*DA;


Ns = 4;

t_alto = 0;
t_reset = 0;
vincitore = 0;
T_reset = 70;
T_alto = 45;
k_tap_vett = [];
caso_neuron = 0.0;
%%
% S0(1) = 1;
% S0(2:Nc) = 0;
% noise1 =   0.2*rand(Nc,1);
for k = 1:D
    
   
   
%     if (t(k) > t_alto) && (vincitore == 1)   || ( (t(k) > t_reset + 500) && (vincitore == 0) ) %I have no response for 200 ms
%     S0(1:2) = S0(1:2)*(-1) +1;   % I change the sign of the inputs, but only after the time t_alto
%     S = S0';
%     vincitore = 0;
%     t_reset = t(k) + T_reset;
%     end
    
%     if t(k) < t_reset
%         C(:,k) = 0;
%         Uc(:,k) = 0;
%     else
        C(:,k) = 1./(1+exp(-a*(Uc(:,k)-U0)));
%     end
    
    if (vincitore == 0) && ( max(C(:,k)) > 0.9)
        t_alto = t(k) + T_alto;
        k_tap_vett = [k_tap_vett t(k)];
        vincitore = 1;
    end
    
    
    Go(:,k) = 1./(1+exp(-a*(Ugo(:,k)-U0)));
    NoGo(:,k) = 1./(1+exp(-a*(Unogo(:,k)-U0)));
    Gpe(:,k) = 1./(1+exp(-a*(Ugpe(:,k)-U0)));
    Gpi(:,k) = 1./(1+exp(-a*(Ugpi(:,k)-U0)));
    ChI(k) = 1./(1+exp(-a*(Uchi(k)-U0)));
    
    if T_ON==1
        T(:,k) = 1./(1+exp(-a*(Ut(:,k)-U0)));
        %otherwise it remains at zero 
    elseif T_ON~=0
        disp('Wrong value for T_ON')
        return
    end
    
    if STN_ON==1
        STN(k) = 1./(1+exp(-a*(Ustn(k)-U0)));
        %otherwise it remains at zero 
    elseif STN_ON~=0
        disp('Wrong value for STN_ON')
        return
    end
    
    
   
    %Energy in the cortex (as an indicator of conflict)
    for i = 1:Nc
        for j = i:Nc
            E(k) = E(k)+C(i,k)*C(j,k);
        end
    end
    E(k) = E(k)-(sum(C(:,k).^2));
    
    %%%% differential equations solved with the Euler method
    Ul(:,k+1) = Ul(:,k)+dt/tauL*(-Ul(:,k)+L*C(:,k));
    Uc(:,k+1) = Uc(:,k)+dt/tau*(-Uc(:,k)+Wcs*S+Ul(:,k)+Wct*T(:,k)+noiseC);
    
    IGo_DA_Ach(:,k) = alpha*DA*(Go(:,k)-0.35)+wgchi*ChI(k); %+0.18; 
    INoGo_DA_Ach(:,k) = (beta*DA+wnchi*ChI(k))*ones(Nc,1);
    
    Ugo(:,k+1) = Ugo(:,k)+dt/tau*(-Ugo(:,k)+Wgs*S+Wgc*C(:,k)+IGo_DA_Ach(:,k));     %versione originaria
    Unogo(:,k+1) = Unogo(:,k)+dt/tau*(-Unogo(:,k)+Wns*S+Wnc*C(:,k)+INoGo_DA_Ach(:,k));
    Ugpe(:,k+1) = Ugpe(:,k)+dt/tau*(-Ugpe(:,k)+Wen*NoGo(:,k)+Westn*STN(k)+Igpe);
    Ugpi(:,k+1) = Ugpi(:,k)+dt/tau*(-Ugpi(:,k)+Wig*Go(:,k)+Wie*Gpe(:,k)+Wistn*STN(k)+Igpi);
    Uchi(k+1) = Uchi(k)+dt/tau*(-Uchi(k)+Ichi+gamma*DA);   %+gammaDop_tonic
    
    
    if T_ON==1
        Ut(:,k+1) = Ut(:,k)+dt/tau*(-Ut(:,k)+Wti*Gpi(:,k)+Wtc*C(:,k));
         %otherwise it remains at zero 
    elseif T_ON~=0
        disp('Wrong value for STN_ON')
        return
    end
    
    if STN_ON==1
        Ustn(k+1) = Ustn(k)+dt/tau*(-Ustn(k)+Ke*E(k)+Kgpe*sum(Gpe(:,k)));
         %otherwise it remains at zero 
    elseif STN_ON~=0
        disp('Wrong value for STN_ON')
        return
    end
    %%%%
    
end


%% computation of the tapping frequency

% if length(k_tap_vett)< 3
%     ft = 0;
% else
% Tt = k_tap_vett(end) - k_tap_vett(end-2);
% ft = 1/Tt*1000;
% end
% frequenza = ft*60





