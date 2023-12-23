                                                                                                                                                                                                                                                                
function [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_new(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,gain_drop_dop)
%% function used to simulate the network during training
% BG_model_function_Ach                         returns dynamical behaviour of Basal Ganglia structures and cortex
% S                                             stimulus. Column vector of 4 elements Ei, 0<=Ei<=1 for i=1,2,3,4
% Wgc,Wgs,Wnc,Wns                               sets corresponding synaptic weights
% C_CORRECT_WINNER                              sets the desired correct response, if possible, depending on the stimulus
% Uc,Ugo,Unogo,Ugpe,Ugpi,Ut,Ustn                returns the input to the sigmoidal function within time of the corresponding brain structures
% C,Go,NoGo,Gpe,Gpi,T,STN,E                     returns activity within time of the corresponding brain structures and energy in the cortex
% IGo_DA_Ach,INoGo_DA_Ach                       returns the input due to Dopa and Ach to Go and NoGo units
% t                                             returns time
% Wgc_post,Wgs_post,Wnc_post,Wns_post           returns corresponding synaptic weights after Hebbian learning
% r                                             returns +1 for reward, -1 for punishment, NaN for no feedback
% r_story                                       retyrn the historical record of rewards and punishments
% k_reward                                      returns position of feedback, NaN for no feedback
% ChI                                           returns activity within time of the cholinegic interneuron

% It is different from the other functions since it use the positive part for all the post-synaptic activities 
% both for Wgs and Wns, and  Wgc and Wnc

%% tempi

tau = 15;       %basal time constant
tauL = 5*tau;   %time constant of lateral inhibition
%tauS = tau/5;

dt = 0.1;   %step
t = (0:dt:800)';   %time [ms]
D = length(t);   %number of samples


%% initialization of the structure

%C: cortex
Nc = 2;   %neurons in the cortex   Nc = 4
C = zeros(Nc,D);   % activity of neurons
Uc = zeros(Nc,D);   %input to the sigmoidal function
Ul = zeros(Nc,D);   %contribution from the lateral inhibition
%Go: striatum, Go
Go = zeros(Nc,D);
Ugo = zeros(Nc,D);
%NoGo: striatum, No-Go
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
%STN: sub-thalamic nucleus
STN = zeros(1,D);
Ustn = zeros(1,D);
%E: energy (it is an index of conflict in the cortex)
E = zeros(1,D);
%ChI: cholinergic interneuron
ChI = zeros(1,D);
Uchi = zeros(1,D);

%input DA+ACh to Go and NoGo
IGo_DA_Ach = zeros(Nc,D);
INoGo_DA_Ach = zeros(Nc,D);


%% initialization of synapses

%weights from stimulus to cortex
%Wcs = 1*ones(4,4);   %0.2 extradiag, 1.1 su diag
Wcs = 1*ones(Nc,Nc); % + 0.2*diag(ones(Nc,1));
%lateral inhibition
L0 = 1.2;   %normale 1.2
L = -L0*ones(Nc,Nc)+diag(L0*ones(Nc,1));   %tolto autoanello
%weights from thalamus to  cortex
Wct = 4*diag(ones(Nc,1));  %diagonale

%%%%%%%%%%%% the following weights are passed as input parameters since subject to learning %%%%%%%%%%%%%%%%%%%%%

% %weights from cortex to Go: Wgc (diagonal)
% %weights from stimulus to Go: Wgs 
%
% %weights from cortex to NoGo: Wnc (diagonal)
% %weights from stimulus to NoGo: Wns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%weights from No-Go to Gpe (inhibitory)
Wen = -2.2*diag(ones(Nc,1));   %diagonal

%weights from Gpe to Gpi (inhibitory)
Wie = -3*diag(ones(Nc,1));   %diagonal
%weights from Go to Gpi (inhibitory)
Wig = -36*diag(ones(Nc,1));   %diagonal

%weights from cortex to thalamus (excitatory)
Wtc = 3*diag(ones(Nc,1));   %diagonale
%weights from Gpi to thalamus (inibitoriy)
Wti = -3*diag(ones(Nc,1));   %diagonal

%%%%%%%%%%%%%%%%%%% weights from-to STN %%%%%%%%%%%%%%%%%%%
%weight from energy to STN (excitatory)
Ke = 7;   % normale 7
%weight from Gpe to STN (inibitory)
Kgpe = -1;

%weight from STN to Gpe (sinapsi excitatory)
Westn = 1;

%weight from STN to Gpi (excitatory)
Wistn = 30;   %14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% weight from ChI to Go-NoGo %%%%%%%%%%%%%%%%%%%
%weight from ChI to Go(inhibitory)
wgchi = -1;

%weight from ChI to NoGo (excitatory)
wnchi = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% gains %%%%%%%%%%%%%%%%%%%
%gain from DA to Go (excitation)
Ugo_trigger = 1.074;
alpha = 0.75;  %0.75/2

%gain from DA to No-Go (inhibition)
beta = -1.0;%-2.0

%gain from DA a Chi (inibition)
gamma = -0.5;   %-0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% computation of the dynamics: low-pass filter  + sigmoid

%parameters of the sigmoid
a = 4;
U0 = 1.0;


%tonic activity Gpe
Igpe = 1.0;
%tonic activity Gpi
Igpi = 3;%3

%tonic activity ChI
Ichi = 1.00;   

%reward o punishment
%r = NaN   no reward o initialization
%r = 1   reward
%r = -1   punishment
r = NaN;
sw = [];
%time of reward o punishment
%k_reward = NaN   no reward o initialization
k_reward = NaN;


%time pattern of phasic dopamine
latency = 100;  %ms, physiological values 50-110 ms (Schultz 1998) %100
klatency = ceil(latency/dt);
duration = 50;   %ms, physiological values <200 ms (Schultz 1998)
kduration = ceil(duration/dt);


%% initial conditions


Ugpe(:,1) = Igpe;
Ugpi(:,1) = Igpi;
Uchi(1) = Ichi+gamma*Dop_tonic;

Ns = 4;

r_story = [];

% noise1 =   0.2*randn(Nc,1);%0.2

gain_ADHD = 1;

for k = 1:D
    
    %PRODUCTION OF PHASIC DOPAMINE
    
    if isnan(r);   %nothing
        dop_phasic = 0;
    else   %k_reward initialization
                if k>=k_reward+klatency && k<=k_reward+klatency+kduration   %time interval of dopamine change
                    if r>=0   %reward
                        delta_Dop = gain_ADHD*Dop_tonic*r;
                    elseif r==-1   %punishment
                        delta_Dop = -gain_drop_dop*Dop_tonic;
                    end
                    dop_phasic = delta_Dop;
                else   %time interval without dopamine change
                    dop_phasic = 0;
                end
    end
% 
    DA = Dop_tonic+dop_phasic;
%     
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% sigmoids
    C(:,k) = 1./(1+exp(-a*(Uc(:,k)-U0)));
    Go(:,k) = 1./(1+exp(-a*(Ugo(:,k)-U0)));
    NoGo(:,k) = 1./(1+exp(-a*(Unogo(:,k)-U0)));
    Gpe(:,k) = 1./(1+exp(-a*(Ugpe(:,k)-U0)));
    Gpi(:,k) = 1./(1+exp(-a*(Ugpi(:,k)-U0)));
    ChI(k) = 1./(1+exp(-a*(Uchi(k)-U0)));
    T(:,k) = 1./(1+exp(-a*(Ut(:,k)-U0)));
    STN(k) = 1./(1+exp(-a*(Ustn(k)-U0)));
    %%%%
    
        %REWARD, PUNISHMENT OR NOTHING
        %I look for a winner; in case, I must produce reward or punishment
        if length(Correct_winner)==0
            Correct_winner = 0;
        end
        if isnan(Correct_winner)==0; 
            if isnan(r);   %I still have no result
                C_winner_list = find(C(:,k)>=0.9);
                lost = isempty(C_winner_list);
                if lost==0   %there is al least one winner
                    n_C_winner = length(C_winner_list);
                    if n_C_winner == 1   %there is just one winner, I give reward or punishment
                        k_reward = k;   %the time of victory
                        if C_winner_list == Correct_winner
                            r = 1;   %reward
                            sw = 1;
                            gain = 0.02;
                        elseif C_winner_list == Small_winner
                            r = 1;
                            sw = 2;
                            gain = 0.02;
                        else
                            r = -1;   %punishment
                            sw = 0;
                            gain = 0.02;
                        end
                        r_story = [r_story; r];
                    end
                end
            end
        end
    
    
    %energy in the cortex (index of a conflict)
    for i = 1:Nc
        for j = i:Nc
            E(k) = E(k)+C(i,k)*C(j,k);
        end
    end
    E(k) = E(k)-(sum(C(:,k).^2));
    
    %%%% differential equations with Euler method
    Ul(:,k+1) = Ul(:,k)+dt/tauL*(-Ul(:,k)+L*C(:,k));
    Uc(:,k+1) = Uc(:,k)+dt/tau*(-Uc(:,k)+Wcs*S+Ul(:,k)+Wct*T(:,k)+noiseC);
    
    IGo_DA_Ach(:,k) = alpha*DA*(Go(:,k)-0.35)+wgchi*ChI(k);
    INoGo_DA_Ach(:,k) = (beta*DA+wnchi*ChI(k))*ones(Nc,1);
    
    Ugo(:,k+1) = Ugo(:,k)+dt/tau*(-Ugo(:,k)+Wgs*S+Wgc*C(:,k)+IGo_DA_Ach(:,k));     
    Unogo(:,k+1) = Unogo(:,k)+dt/tau*(-Unogo(:,k)+Wns*S+Wnc*C(:,k)+INoGo_DA_Ach(:,k));
    Ugpe(:,k+1) = Ugpe(:,k)+dt/tau*(-Ugpe(:,k)+Wen*NoGo(:,k)+Westn*STN(k)+Igpe);
    Ugpi(:,k+1) = Ugpi(:,k)+dt/tau*(-Ugpi(:,k)+Wig*Go(:,k)+Wie*Gpe(:,k)+Wistn*STN(k)+Igpi);
    Uchi(k+1) = Uchi(k)+dt/tau*(-Uchi(k)+Ichi+gamma*DA);
    Ut(:,k+1) = Ut(:,k)+dt/tau*(-Ut(:,k)+Wti*Gpi(:,k)+Wtc*C(:,k));
    Ustn(k+1) = Ustn(k)+dt/tau*(-Ustn(k)+Ke*E(k)+Kgpe*sum(Gpe(:,k)));
    %%%%
    
end

%%  HEBB RULE
%inizialization
S_th = 0.5*ones(length(S),1);
C_th = 0.5*ones(Nc,1);
Go_th = 0.5*ones(Nc,1);%0.5 0.4
NoGo_th = 0.5*ones(Nc,1);%0.5
% gain = 0.02;


delta_Wgc = zeros(Nc,Nc);
delta_Wgs = zeros(Nc,Nc);

delta_Wnc = zeros(Nc,Nc);
delta_Wns = zeros(Nc,Nc);


if isnan(r)==0   % we had reward o punishment
    
    if ((k_reward+klatency)>D)==0 &&((k_reward+klatency+kduration)>D)==0   %I do not train the net if still varying
%% cambio
if r>=0%r==1   %reward
    %maximum of the winner Go and minimum of the loosing Go
    %minimum of all NoGO
    val_Go(1) = min(Go(1,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(2) = min(Go(2,(k_reward+klatency):(k_reward+klatency+kduration)));
%     val_Go(3) = min(Go(3,(k_reward+klatency):(k_reward+klatency+kduration)));
%     val_Go(4) = min(Go(4,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(C_winner_list) = max(Go(C_winner_list,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(1) = min(NoGo(1,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(2) = min(NoGo(2,(k_reward+klatency):(k_reward+klatency+kduration)));
%     val_NoGo(3) = min(NoGo(3,(k_reward+klatency):(k_reward+klatency+kduration)));
%     val_NoGo(4) = min(NoGo(4,(k_reward+klatency):(k_reward+klatency+kduration)));
elseif r<0   %punishment
    %minimum of all GO
    %maximum of all NoGo
    val_Go(1) = min(Go(1,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(2) = min(Go(2,(k_reward+klatency):(k_reward+klatency+kduration)));
%     val_Go(3) = min(Go(3,(k_reward+klatency):(k_reward+klatency+kduration)));
%     val_Go(4) = min(Go(4,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(1) = max(NoGo(1,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(2) = max(NoGo(2,(k_reward+klatency):(k_reward+klatency+kduration)));
%     val_NoGo(3) = max(NoGo(3,(k_reward+klatency):(k_reward+klatency+kduration)));
%     val_NoGo(4) = max(NoGo(4,(k_reward+klatency):(k_reward+klatency+kduration)));
end


        
%%
        for l=1:Nc
            delta_Wgc(l,l) = gain*(C(l,k_reward)-C_th(l)).*max(0,val_Go(l)-Go_th(l));
            delta_Wnc(l,l) = gain*(C(l,k_reward)-C_th(l)).*max(0,(val_NoGo(l)-NoGo_th(l)));
        end
        
        for row=1:Nc
            for col=1:Nc
                delta_Wgs(row,col) = gain*max(0,(val_Go(row)-Go_th(row)))*(S(col)-S_th(col));
                delta_Wns(row,col) = gain*max(0,(val_NoGo(row)-NoGo_th(row)))*(S(col)-S_th(col));
            end
        end
    end
    
end

Wgc_post = Wgc+delta_Wgc;
Wgs_post = Wgs+delta_Wgs;
Wnc_post = Wnc+delta_Wnc;
Wns_post = Wns+delta_Wns;

%% CONTROL ON SYNAPSES

%synapse saturation
delta_opt = 0.2;

sinapsi_ecc_max = 1.2;
max_ecc_go = sinapsi_ecc_max-2*delta_opt;
max_ecc_nogo = sinapsi_ecc_max+delta_opt;

Wgc_opt_max = max_ecc_go*diag(ones(Nc,1));   %diagonal
Wgc_opt_min = zeros(Nc,Nc);
Wnc_opt_max = max_ecc_nogo*diag(ones(Nc,1));   %diagonal
Wnc_opt_min = zeros(Nc,Nc);
Wgs_opt_max = max_ecc_go*ones(Nc,Nc);   %full
Wgs_opt_min = zeros(Nc,Nc);
Wns_opt_max = max_ecc_nogo*ones(Nc,Nc);   %full
Wns_opt_min = zeros(Nc,Nc);

[rowgc,colgc,valgc] = find(Wgc_post>Wgc_opt_max);
if isempty(rowgc)==0
    for i=1:length(rowgc)
        Wgc_post(rowgc(i),colgc(i)) = Wgc_opt_max(rowgc(i),colgc(i));
    end
end
[rowgc,colgc,valgc] = find(Wgc_post<Wgc_opt_min);
if isempty(rowgc)==0
    for i=1:length(rowgc)
        Wgc_post(rowgc(i),colgc(i)) = Wgc_opt_min(rowgc(i),colgc(i));
    end
end

[rownc,colnc,valnc] = find(Wnc_post>Wnc_opt_max);
if isempty(rownc)==0
    for i=1:length(rownc)
        Wnc_post(rownc(i),colnc(i)) = Wnc_opt_max(rownc(i),colnc(i));
    end
end
[rownc,colnc,valgnc] = find(Wnc_post<Wnc_opt_min);
if isempty(rownc)==0
    for i=1:length(rownc)
        Wnc_post(rownc(i),colnc(i)) = Wnc_opt_min(rownc(i),colnc(i));
    end
end

[rowgs,colgs,valgs] = find(Wgs_post>Wgs_opt_max);
if isempty(rowgs)==0
    for i=1:length(rowgs)
        Wgs_post(rowgs(i),colgs(i)) = Wgs_opt_max(rowgs(i),colgs(i));
    end
end
[rowgs,colgs,valgs] = find(Wgs_post<Wgs_opt_min);
if isempty(rowgs)==0
    for i=1:length(rowgs)
        Wgs_post(rowgs(i),colgs(i)) = Wgs_opt_min(rowgs(i),colgs(i));
    end
end

[rowns,colns,valns] = find(Wns_post>Wns_opt_max);
if isempty(rowns)==0
    for i=1:length(rowns)
        Wns_post(rowns(i),colns(i)) = Wns_opt_max(rowns(i),colns(i));
    end
end
[rowns,colns,valns] = find(Wns_post<Wns_opt_min);
if isempty(rowns)==0
    for i=1:length(rowns)
        Wns_post(rowns(i),colns(i)) = Wns_opt_min(rowns(i),colns(i));
    end
end
