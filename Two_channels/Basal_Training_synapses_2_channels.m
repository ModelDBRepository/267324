%==========================================================================
%Title: Phasic dopamine changes and Hebbian mechanisms during probabilistic
%reversal learning in striatal circuits: a computational study
% 
%Description: Script for training synapses during the Basal Learning phase.
%
% 10 independent learning procedures (each with 200 epochs starting from
% the naïve network and different noise realizations, with different seeds
% of noise) are performed.
%
% Each epoch consists of the two stimuli, S1 = [1 0.3] and S2 = [0.3 1].
%
% In response to the first stimulus, S1, the first action is randomly
% rewarded with a probability P_high as high as 0.8 (otherwise the action
% is punished), while the second action is randomly rewarded with a
% probability P_small = 0.3 as low as 0.3 (otherwise punished). When a
% second stimulus S2 is given the reward probabilities are P_small for the
% first choice and P_high for the second one.
% 
% The user chooses the Hebb rule by entering the corresponding value.
% 
% The results are saved in a .mat file which indicates:
%    The name of the rule used; The number of the learning procedure, with
%    ii procedure index.
% (ex: Post Post Rule , ii = 1 'W_tot_post_post_1.mat')
%
% The results can be plotted by the program 'plot_training_synapses_2_ch.m'
%
% Mauro Ursino, Miriam Schirru 
% Jan. 2022
%==========================================================================

clear all
close all
clc
 %%
% %--------------------------------------------------------------------------
% % Initialisation 
% %--------------------------------------------------------------------------

P_high = 0.8;  
P_small = 0.3; 

answer = input('choose the rule (1 = post-post, 2 = post-pre, 3 = pre-pre, 4 = total, 5 = Oja, 6 = post-adapt) ');
%% basal stimuli
% Ns stimuli may be applied to the network
Ns = 2;
Nc = 2;
S_high = 1.0;
S_small = 0.3;

%S1: stimulus 1
S1(1) = S_high;
S1(2) = S_small;
% S1(3) = 0.;
% S1(4) = 0.;
S1 = S1';

Correct_winner_1 = 1;
Small_winner_1 = 2;

%S2: stimulus 2
S2(1) = S_small;
S2(2) = S_high;
% S2(3) = 0.;
% S2(4) = 0.;
S2 = S2';

Correct_winner_2 = 2;
Small_winner_2 = 1;

%S3: stimulus 3
S3(1) = 0.1;
S3(2) = 0.1;
S3(3) = S_high;
S3(4) = S_small;
S3 = S3';

Correct_winner_3 = 3;
Small_winner_3 = 4;

%S4:stimulus 4
S4(1) = 0.1;
S4(2) = 0.1;
S4(3) = S_small;
S4(4) = S_high;
S4 = S4';

Correct_winner_4 = 4;
Small_winner_4 = 3;

N_training = 10;
N_epochs = 200;
Dop_tonic = 1.0; % value of the dopaminergic input used during training, default 1.0
gain_drop_dop = 0.5;

%----------------------------------------------------------------------
% Training synapses
%----------------------------------------------------------------------
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

    rng(10+ii)
    noise1=zeros(Nc,Ns*N_epochs);
    noise1(1:Ns,1:Ns*N_epochs) =  0.20*randn(Ns,Ns*N_epochs);

    par = 0;
    %% initial value of synapses before learning
    % default value = 0.5 before learning
    % weights from cortex to GO
    Wgc = 0.44*diag(ones(Nc,1));%0.44
    %     Wgc(3,3) = 0;
    %     Wgc(4,4) = 0;

    %  weights from stimuli to GO
    Wgs = 0.44*(ones(Nc,Nc));%0.44
    %     Wgs(1,2) = 0.5;
    %     Wgs(2,1) = 0.5;

    %  weights from cortex to NOGO
    Wnc = 0.82*diag(ones(Nc,1));%0.82
    %     Wnc(3,3) = 0;
    %     Wnc(4,4) = 0;

    %  weights from stimuli to NOGO
    Wns = 0.82*(ones(Nc,Nc));%0.82
    %     Wns(1,2) = 0.5;
    %     Wns(2,1) = 0.5;

    Wgc_epocs = zeros(Nc,Nc,Ns*N_epochs);
    Wgs_epocs= zeros(Nc,Nc,Ns*N_epochs);
    Wnc_epocs = zeros(Nc,Nc,Ns*N_epochs);
    Wns_epocs = zeros(Nc,Nc,Ns*N_epochs);

    Wgc_epocs(:,:,1) = Wgc;
    Wgs_epocs(:,:,1) = Wgs;
    Wnc_epocs(:,:,1) = Wnc;
    Wns_epocs(:,:,1) = Wns;

    vett_reward = zeros(Ns*N_epochs,1);
    vett_punishment = zeros(Ns*N_epochs,1);
    vett_small_winn = zeros(Ns*N_epochs,1);
    vett_no_risposta = zeros(Ns*N_epochs,1);

    S_vett = zeros(Nc,Ns*N_epochs);
    S1_vett = zeros(Nc,Ns*N_epochs);
    S2_vett = zeros(Nc,Ns*N_epochs);

    exitC_tot = zeros(Nc,Ns*N_epochs);
    %%
    for i = 1:N_epochs

        action = randperm(Ns);

        for j = 1:Ns
            Wgc = squeeze(Wgc_epocs(:,:,j+Ns*(i-1)));
            Wgs = squeeze(Wgs_epocs(:,:,j+Ns*(i-1)));
            Wnc = squeeze(Wnc_epocs(:,:,j+Ns*(i-1)));
            Wns = squeeze(Wns_epocs(:,:,j+Ns*(i-1)));
            %       resto = rem(i,4);% try with switch
            noise = 0.0*randn(Nc,1);
            noiseC = noise1(:,j+Ns*(i-1));
            num_rand = rand(1,1);
            Correct_winner = [];
            Small_winner = [];


            switch action (j)

                case 1
                    S = S1;
                    S(1) = S(1)+noise(1);
                    S(2) = S(2)+noise(2);

                    if num_rand < P_high
                        Correct_winner = Correct_winner_1;
                    end
                    if num_rand < P_small
                        Correct_winner = Correct_winner_1;
                        Small_winner = Small_winner_1;
                    end

                case 2
                    S = S2;
                    S(1) = S(1)+noise(1);
                    S(2) = S(2)+noise(2);

                    if num_rand < P_high
                        Correct_winner = Correct_winner_2;
                    end
                    if num_rand < P_small
                        Correct_winner = Correct_winner_2;
                        Small_winner = Small_winner_2;
                    end
            end


            S(find(S>1)) = 1;
            S(find(S<0)) = 0;
            %     S(3:4) = 0;

            switch answer
                % Call to the function which simulates the basal ganglia response
                case 1
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_post_post(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,gain_drop_dop);
                case 2
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_post_pre(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,gain_drop_dop);
                case 3
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_pre_pre(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,gain_drop_dop);
                case 4
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_total(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,gain_drop_dop);
                case 5
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_Oja(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,gain_drop_dop);
                case 6
                    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_post_adapt(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,gain_drop_dop);


            end
            exitC(1,j) = C(1,end); %exit cortex ch 1
            exitC(2,j) = C(2,end);

            [i r]

            if r==1 & sw==1
                vett_reward(j+Ns*(i-1)) = 1;
            elseif r==-1 & sw==0
                vett_punishment(j+Ns*(i-1)) = 1;
            elseif r==1 & sw==2
                vett_small_winn(j+Ns*(i-1)) = 1;
            else
                vett_no_risposta(j+Ns*(i-1)) = 1;
            end

            S_vett(1,j+Ns*(i-1)) = S(1);
            S_vett(2,j+Ns*(i-1)) = S(2);
            %
            Wgc_epocs(:,:,j+Ns*(i-1)+1) = Wgc_post;
            Wgs_epocs(:,:,j+Ns*(i-1)+1) = Wgs_post;
            Wnc_epocs(:,:,j+Ns*(i-1)+1) = Wnc_post;
            Wns_epocs(:,:,j+Ns*(i-1)+1) = Wns_post;

        end

        exitC_tot (:,Ns*i-(Ns-1):Ns*i) = exitC;

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clear Wgc Wgs Wnc Wns

    end

    %choices
    v = zeros(1,Ns*N_epochs);
    idx_corr_act = find((~(sum (exitC_tot>=0.9)> 1)) & (~(sum(exitC_tot < 0.9)>1))  );
    idx_wrong_act = find((sum (exitC_tot>=0.9)> 1) | (sum(exitC_tot < 0.9)>1));
    v (idx_wrong_act) = NaN;

    [row,col] = find(exitC_tot(:,idx_corr_act) >= 0.9);
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

    % save the date to the file named W_tot_new. This name can be subsequently changed
    save (name, 'Wgc_epocs', 'Wgs_epocs', 'Wnc_epocs', 'Wns_epocs', 'vett_reward', 'vett_punishment', 'vett_no_risposta', 'vett_small_winn', 'S_vett', 'Dop_tonic', 'N_epochs','exitC_tot','choices')

end


