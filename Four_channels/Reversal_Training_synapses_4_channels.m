%==========================================================================
%Title: Phasic dopamine changes and Hebbian mechanisms during probabilistic
%reversal learning in striatal circuits: a computational study
%
%Description:Script for training synapses: Reversal Learning.
%
% To be run after Basal_Training_synapses_4_channels.m The script loads
% simulations result obtained from Basal_Training_synapses_4_channels.m
%
% Each epoch consists of the four possibile stimuli. During reversal
% learning, the probabilities were inverted.
%
% Reversal learning was attempted at different moments: after 100 epochs
% and after 150 epochs of previous learning. To be indicated in line 123.
%
% The user chooses the Hebb rule by entering the corresponding value.
%
% The results are saved in a .mat file which indicates:
%    The type of learning (Rev = Reversal Learning); The name of the rule
%    used; The number of the learning procedure, with ii procedure index.
% (ex: Reversal Learning, Post Post Rule , ii = 1
% 'Rev_W_tot_post_post_1.mat')
%
% The results can be shown by running the program
% 'plot_training_synapses_4_ch.m'
%
% Mauro Ursino, Miriam Schirru 
% Jan. 2022
%==========================================================================

clear all
close all
clc

oldpath = path;
selpath = uigetdir;
path(selpath,oldpath)

        P_high = 0.8;  % probability of reward for the best choice
        P_small = 0.3; % probability of reward for the small choice
        P_low = 0.1;

answer = input('scegli regola (1 = post-post, 2 = post-pre, 3 = pre-pre, 4 = total, 5 = Oja, 6 = adapt)');

%% basal stimuli

% four stimuli may be applied to the network
Ns = 4;
Nc = 4;
N_training = 10;

S_high = 1.0;
S_small = 0.3;

%S1: stimulus 1
S1(1) = S_high;   
S1(2) = S_small;   
S1(3) = 0.1;
S1(4) = 0.1;
S1 = S1';

Correct_winner_1 = 2;
Small_winner_1 = 1;
Low_winner_1 = 3;
Low_winner_2 = 4;

%S2: stimulus 2
S2(1) = S_small;   
S2(2) = S_high;   
S2(3) = 0.1;
S2(4) = 0.1;
S2 = S2';

Correct_winner_2 = 1;
Small_winner_2 = 2;
Low_winner_1 = 4;
Low_winner_2 = 3;

%S3: stimulus 3
S3(1) = 0.1;   
S3(2) = 0.1;   
S3(3) = S_high;
S3(4) = S_small;
S3 = S3';

Correct_winner_3 = 4;
Small_winner_3 = 3;
Low_winner_1 = 1;
Low_winner_2 = 2;

%S4:stimulus 4
S4(1) = 0.1;   
S4(2) = 0.1;   
S4(3) = S_small;
S4(4) = S_high;
S4 = S4';

Correct_winner_4 = 3;
Small_winner_4 = 4;
Low_winner_1 = 2;
Low_winner_2 = 1;


%%
for ii = 1:N_training
    
    switch answer
        case 1
            name = strcat('W_tot_post_post_',num2str(ii));
        case 2
            name = strcat('W_tot_post_pre_',num2str(ii));
        case 3
            name = strcat('W_tot_pre_pre_',num2str(ii));
        case 4
            name = strcat('W_tot_total_',num2str(ii));
        case 5
            name = strcat('W_tot_oja_',num2str(ii));
        case 6
            name = strcat('W_tot_post_adapt_',num2str(ii));
    end
    
load(name)
Epoca = 150; %100
Wgc = squeeze(Wgc_epocs(:,:,4*Epoca));

Wgs = squeeze(Wgs_epocs(:,:,4*Epoca));
Wnc = squeeze(Wnc_epocs(:,:,4*Epoca));
Wns = squeeze(Wns_epocs(:,:,4*Epoca));

N_epoche = 400;  
Dop_tonic = 1.0; % value of the dopaminergic input used during training, default 1.0
gain_drop_dop = 0.5; %1.0 %1.5

Wgc_epocs = zeros(Nc,Nc,4*N_epoche);
Wgs_epocs= zeros(Nc,Nc,4*N_epoche);
Wnc_epocs = zeros(Nc,Nc,4*N_epoche);
Wns_epocs = zeros(Nc,Nc,4*N_epoche);

Wgc_epocs(:,:,1) = Wgc;
Wgs_epocs(:,:,1) = Wgs;
Wnc_epocs(:,:,1) = Wnc;
Wns_epocs(:,:,1) = Wns;

vett_reward = zeros(4*N_epoche,1);
vett_punishment = zeros(4*N_epoche,1);
vett_small_winn = zeros(4*N_epoche,1);
vett_low_winn_1 = zeros(Ns*N_epoche,1);
vett_low_winn_2 = zeros(Ns*N_epoche,1);
vett_no_risposta = zeros(4*N_epoche,1);

S_vett = zeros(Nc,4*N_epoche);
S1_vett = zeros(Nc,4*N_epoche);
S2_vett = zeros(Nc,4*N_epoche);
S3_vett= zeros(Nc,4*N_epoche);
S4_vett= zeros(Nc,4*N_epoche);

exitC_tot = zeros(Nc,Ns*N_epoche);

rng(15)%12
noise1 = zeros(Nc,Ns*N_epoche);
noise1(1:Ns,1:Ns*N_epoche) = 0.20*randn(Ns,Ns*N_epoche);% rand
%%
for i = 1:N_epoche
    
    action = randperm (Ns);
    
    for j = 1:Ns
        Wgc = squeeze(Wgc_epocs(:,:,j+4*(i-1)));
        Wgs = squeeze(Wgs_epocs(:,:,j+4*(i-1)));
        Wnc = squeeze(Wnc_epocs(:,:,j+4*(i-1)));
        Wns = squeeze(Wns_epocs(:,:,j+4*(i-1)));
%       resto = rem(i,4);% try with switch
        noise = 0.0*randn(4,1);
        noiseC = noise1(:,j+4*(i-1));
        num_rand = rand(1,1);
        Correct_winner = [];
        Small_winner = [];
        Low_winner_1 = [];
        Low_winner_2 = [];

        % if the sum is less than 1, the remaining probability id divided
        % between the remaining chices
        
        switch action (j)
            
            case 1
                S = S1;
                S(1) = S(1)+noise(1);
                S(2) = S(2)+noise(2);
                S(3) = S(3)+noise(3);
                S(4) = S(4)+noise(4);
                
                if num_rand < P_high
                    Correct_winner = Correct_winner_1;
                end
                if num_rand < P_small
                    Correct_winner = Correct_winner_1; 
                    Small_winner = Small_winner_1;
                end
                if num_rand < P_low
                    Correct_winner = Correct_winner_1;
                    Small_winner = Small_winner_1;
                    Low_winner_1 = 3;
                    Low_winner_2 = 4;
                end
                
            case 2
                S = S2;
                S(1) = S(1)+noise(1);
                S(2) = S(2)+noise(2);
                S(3) = S(3)+noise(3);
                S(4) = S(4)+noise(4);
                
                if num_rand < P_high
                    Correct_winner = Correct_winner_2;
                end
                if num_rand < P_small
                    Correct_winner = Correct_winner_2;
                    Small_winner = Small_winner_2;
                end
                if num_rand < P_low
                    Correct_winner = Correct_winner_2;
                    Small_winner = Small_winner_2;
                    Low_winner_1 = 4;
                    Low_winner_2 = 3;
                end
                
            case 3
                S = S3;
                S(1) = S(1)+noise(1);
                S(2) = S(2)+noise(2);
                S(3) = S(3)+noise(3);
                S(4) = S(4)+noise(4);
                
                if num_rand < P_high
                    Correct_winner = Correct_winner_3;
                end
                if num_rand < P_small
                    Correct_winner = Correct_winner_3;
                    Small_winner = Small_winner_3;
                end
                if num_rand < P_low
                    Correct_winner = Correct_winner_3;
                    Small_winner = Small_winner_3;
                    Low_winner_1 = 1;
                    Low_winner_2 = 2;
                end
                
            case 4  
                S = S4;
                S(1) = S(1)+noise(1);
                S(2) = S(2)+noise(2);
                S(3) = S(3)+noise(3);
                S(4) = S(4)+noise(4);
                
                if num_rand < P_high
                    Correct_winner = Correct_winner_4;
                end
                if num_rand < P_small
                    Correct_winner = Correct_winner_4; 
                    Small_winner = Small_winner_4;
                end
                if num_rand < P_low
                    Correct_winner = Correct_winner_4;
                    Small_winner = Small_winner_4;
                    Low_winner_1 = 2;
                    Low_winner_2 = 1;
                end
    end

    
    S(find(S>1)) = 1;
    S(find(S<0)) = 0;
%     S(3:4) = 0;
    
    % Call to the function which simulates the basal ganglia response
    switch answer
                        case 1
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_post_post(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Low_winner_1,Low_winner_2,Dop_tonic,noiseC,gain_drop_dop);
                case 2
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_post_pre(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Low_winner_1,Low_winner_2,Dop_tonic,noiseC,gain_drop_dop);
                case 3
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_pre_pre(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Low_winner_1,Low_winner_2,Dop_tonic,noiseC,gain_drop_dop);
                case 4
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_total(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Low_winner_1,Low_winner_2,Dop_tonic,noiseC,gain_drop_dop);
                case 5
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_Oja(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Low_winner_1,Low_winner_2,Dop_tonic,noiseC,gain_drop_dop);
                case 6
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_post_adapt(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Low_winner_1,Low_winner_2,Dop_tonic,noiseC,gain_drop_dop);

    end
    exitC(1,j) = C(1,end);%C(1,end); %exit cortex ch 1
    exitC(2,j) = C(2,end);
    exitC(3,j) = C(3,end);
    exitC(4,j) = C(4,end);
    
    [i r] 

    
    if r==1 & sw==1
        vett_reward(j+4*(i-1)) = 1;
    elseif r==-1 & sw==0
        vett_punishment(j+4*(i-1)) = 1;
    elseif r==1 & sw==2
        vett_small_winn(j+4*(i-1)) = 1;
    elseif r==1 & sw==3
        vett_low_winn_1(j+4*(i-1)) = 1;
    elseif r==1 & sw==4
        vett_low_winn_2(j+4*(i-1)) = 1;
    else
        vett_no_risposta(j+4*(i-1)) = 1;
    end
    
    S_vett(1,j+4*(i-1)) = S(1);
    S_vett(2,j+4*(i-1)) = S(2);
    S_vett(3,j+4*(i-1)) = S(3);
    S_vett(4,j+4*(i-1)) = S(4);
    
    Wgc_epocs(:,:,j+4*(i-1)+1) = Wgc_post;
    Wgs_epocs(:,:,j+4*(i-1)+1) = Wgs_post;
    Wnc_epocs(:,:,j+4*(i-1)+1) = Wnc_post;
    Wns_epocs(:,:,j+4*(i-1)+1) = Wns_post;
    
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    exitC_tot (:,Ns*i-(Ns-1):Ns*i) = exitC;
    
    clear Wgc Wgs Wnc Wns
    
end 
 %choices
    v = zeros(1,Ns*N_epoche);
    
    idx_corr_act = find((~(sum (exitC_tot>=0.9)> 1)) & (~(sum(exitC_tot < 0.9)> Ns-1))  );
    idx_wrong_act = find((sum (exitC_tot>=0.9)> 1) | (sum(exitC_tot < 0.9)> Ns-1));
    v (idx_wrong_act) = NaN;
    
    [row,col] = find(exitC_tot(:,idx_corr_act) > 0.9);
    v (idx_corr_act) = row;
    for jj = 1 : length(v)
        if isnan (v(jj))
            choices (jj)= NaN;
        else
            choices(jj) = S_vett(v(jj),jj);
        end
    end
    
%rewards
big_reward_tot = sum(vett_reward)
%punishments
punishment_tot = sum(vett_punishment)
%no answers
no_answer_tot = sum(vett_no_risposta)
%small reward
small_reward_tot = sum(vett_small_winn)
%low reward
low_reward_tot_1 = sum(vett_low_winn_1)
low_reward_tot_2 = sum(vett_low_winn_2)


% save the date to the file named W_tot_new. This name can be subsequently changed 
save (strcat('Rev_',name), 'Wgc_epocs', 'Wgs_epocs', 'Wnc_epocs', 'Wns_epocs', 'vett_reward', 'vett_punishment', 'vett_no_risposta', 'vett_small_winn','vett_low_winn_1','vett_low_winn_2','S_vett', 'Dop_tonic', 'N_epoche','choices')

end

