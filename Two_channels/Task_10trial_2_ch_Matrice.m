%==========================================================================
%Title: Phasic dopamine changes and Hebbian mechanisms during probabilistic
%reversal learning in striatal circuits: a computational study
%
%Description: Script to perform the testing procedure: 2 choices
%experiment.
% To be run after 'Basal_Training_synapses_2_channels.m' or  Reversal_Training_synapses_2_channels.m’.
% The script calls the function used to simulate the task named
% 'BG_model_Response_Stimuli_Task1'.
% 
% 100 trials are executed, by using the inputs S1 and S2 50 times each.
% During the test, random Gaussian noise was applied not only to cortical
% neurons but also to the input stimuli (simulating a noisy environment).
%
% This testing procedure has been performed starting from the network
% obtained after different epochs of learning indicated by vectors in lines
% 46-47.
%
% The results are shown at the end of the simulation in a table organized
% as follow: 
% 1st column: epoch 
% 2nd column: mean 1st Action 
% 3rd column: standard deviation 1st Action 
% 4th column: mean 2nd Action 
% 5th column: standard deviation 2nd Action
%
% Mauro Ursino, Miriam Schirru Jan. 2022
%==========================================================================
clc
clear
close all

oldpath = path;
selpath = uigetdir;
path(selpath,oldpath)
answer = input('choose the rule (1 = post-post, 2 = post-pre, 3 = pre-pre, 4 = total, 5 = Oja, 6 = adapt )');
answer1 = input ('choose the training (1 = basal, 2 = reversal)');

% %--------------------------------------------------------------------------
% % Initialisation 
% %--------------------------------------------------------------------------
Nprove = 10;
Ns = 2;
Nc = 2;

% Epoche = [1 50 75 100 150 200 250 300 350 400];  % case of reversal
Epoche =  [1 25 50 75 100 150 200 ];   % case basal
L_epoche = length(Epoche);

Matrice_risultati = zeros(L_epoche,5);

for indice_epoca = 1 : L_epoche

    Epoca = Epoche(indice_epoca);

    Array_n_succ_tot = zeros(Nc,Nprove);
    Array_n_small_succ_tot = zeros(Nc,Nprove);
    S1 = zeros(Nc,1);
    S2 = zeros(Nc,1);
    % S3 = zeros(Nc,1);
    % S4 = zeros(Nc,1);

    for ii = 1:Nprove

        switch answer
            case 1
                if answer1 == 1
                    name = strcat('W_tot_post_post_',num2str(ii));
                elseif answer1 == 2
                    name = strcat('Rev_W_tot_post_post_',num2str(ii));
                end
            case 2
                if answer1 == 1
                    name = strcat('W_tot_post_pre_',num2str(ii));
                elseif answer1 == 2
                    name = strcat('Rev_W_tot_post_pre_',num2str(ii));
                end
            case 3
                if answer1 == 1
                    name = strcat('W_tot_pre_pre_',num2str(ii));
                elseif answer1 == 2
                    name = strcat('Rev_W_tot_pre_pre_',num2str(ii));
                end
            case 4
                if answer1 == 1
                    name = strcat('W_tot_total_',num2str(ii));
                elseif answer1 == 2
                    name = strcat('Rev_W_tot_total_',num2str(ii));
                end
            case 5
                if answer1 == 1
                    name = strcat('W_tot_oja_',num2str(ii));
                elseif answer1 == 2
                    name = strcat('Rev_W_tot_oja_',num2str(ii));
                end
            case 6
                if answer1 == 1
                    name = strcat('W_tot_post_adapt_',num2str(ii));
                elseif answer1 == 2
                    name = strcat('Rev_W_tot_post_adapt_',num2str(ii));
                end
        end

        load(name)

        % basal stimuli
        caso = 0;

        STN_ON = 1;
        T_ON = 1;

        %trained synapses
        Wgc = squeeze(Wgc_epocs(:,:,2*Epoca));
        Wgs = squeeze(Wgs_epocs(:,:,2*Epoca));
        Wnc = squeeze(Wnc_epocs(:,:,2*Epoca));
        Wns = squeeze(Wns_epocs(:,:,2*Epoca));

        if answer1 == 1
            S_high = 1.0;
            S_small = 0.3;
        elseif answer1 == 2
            S_high = 0.3;
            S_small = 1.0;
        end

        %S1: stimulus 1
        S1(1) = S_high;
        S1(2) = S_small;
        % S1(3) = 0.1;
        % S1(4) = 0.1;
        % % S1 = S1';

        Correct_winner_1 = 1;
        Small_winner_1 = 2;

        %S2: stimulus 2
        S2(1) = S_small;
        S2(2) = S_high;
        % S2(3) = 0.1;
        % S2(4) = 0.1;
        % S2 = S2';

        Correct_winner_2 = 2;
        Small_winner_2 = 1;

        %S3: stimulus 3
        S3(1) = 0.1;
        S3(2) = 0.1;
        S3(3) = S_high;
        S3(4) = S_small;
        % S3 = S3';

        Correct_winner_3 = 3;
        Small_winner_3 = 4;

        %S4:stimulus 4
        S4(1) = 0.1;
        S4(2) = 0.1;
        S4(3) = S_small;
        S4(4) = S_high;
        % S4 = S4';

        Correct_winner_4 = 4;
        Small_winner_4 = 3;


        % task parameters

        trials = 100/Ns;
        Dop_tonic = 1.0; % value of the dopaminergic input used during training, default 1.2

        S_vett = [];
        exitC = [];
        Winner = [];

        S_vett_tot = zeros(Nc,Ns*trials);
        exitC_tot = zeros(Nc,Ns*trials);
        Winner_tot = zeros(1,Ns*trials);
        n_err_act_tot = zeros(1,Ns*trials);

        rng(21)
        noise1=zeros(Nc,Ns*trials);
        noise1(1:Ns,1:Ns*trials) =  0.15*randn(Ns,Ns*trials);
        % noise1 =   0.20*randn(Nc,4*trials);% noise to the cortex

        rng(31)
        noise2=zeros(Nc,Ns*trials);
        noise2(1:Ns,1:Ns*trials) = 0.2*randn(Ns,Ns*trials);% noise to S
        %%
        for i = 1:trials

            action = randperm (Ns);

            for j = 1:Ns
                noiseS = noise2(:,j+Ns*(i-1));
                noiseC = noise1(:,j+Ns*(i-1));
                %
                switch action(j)

                    case 1
                        S = S1;
                        S(1) = S(1)+noiseS(1);
                        S(2) = S(2)+noiseS(2);
                        %                 S(3) = S(3)+noiseS(3);
                        %                 S(4) = S(4)+noiseS(4);
                        Correct_winner = Correct_winner_1;
                        Minor_winner = Small_winner_1;

                    case 2
                        S = S2;
                        S(1) = S(1)+noiseS(1);
                        S(2) = S(2)+noiseS(2);
                        %                 S(3) = S(3)+noiseS(3);
                        %                 S(4) = S(4)+noiseS(4);
                        Correct_winner = Correct_winner_2;
                        Minor_winner = Small_winner_2;

                    case 3
                        S = S3;
                        S(1) = S(1)+noiseS(1);
                        S(2) = S(2)+noiseS(2);
                        S(3) = S(3)+noiseS(3);
                        S(4) = S(4)+noiseS(4);
                        Correct_winner = Correct_winner_3;
                        Minor_winner = Small_winner_3;

                    case 4
                        S = S4;
                        S(1) = S(1)+noiseS(1);
                        S(2) = S(2)+noiseS(2);
                        S(3) = S(3)+noiseS(3);
                        S(4) = S(4)+noiseS(4);
                        Correct_winner = Correct_winner_4;
                        Minor_winner = Small_winner_4;
                end


                %%
                S(find(S>1)) = 1;
                S(find(S<0)) = 0;
                %     S(3:4) = 0;


                % Call to the function which simulates the basal ganglia response
                [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,k_tap_vett,Uchi,ChI,t] = BG_model_Response_Stimuli_Task1(S,Wgc,Wgs,Wnc,Wns,STN_ON,T_ON,Dop_tonic,noiseC);


                S_vett(1,j) = S(1);
                S_vett(2,j) = S(2);
                %     S_vett(3,j) = S(3);
                %     S_vett(4,j) = S(4);

                Winner(j) = Correct_winner;
                Small_winner(j) = Minor_winner;

                exitC(1,j) = C(1,end); %exit cortex ch 1
                exitC(2,j) = C(2,end);

                if max(exitC(:,j)) < 0.9 | sum (exitC(:,j)>0.9)>1 % number of no act or multiple act
                    n_err_act(j) = 1;
                    exitC(:,j) = zeros(Nc,1); %multiple act or no act
                else
                    n_err_act(j)= NaN;
                end

            end

            S_vett_tot(:,Ns*i-(Ns-1):Ns*i) = S_vett;
            exitC_tot (:,Ns*i-(Ns-1):Ns*i) = exitC;
            Winner_tot(1,Ns*i-(Ns-1):Ns*i) = Winner;
            Small_winner_tot(1,Ns*i-(Ns-1):Ns*i) = Small_winner;

            err = isempty(n_err_act);
            if err ==0
                n_err_act_tot(1,Ns*i-(Ns-1):Ns*i) = n_err_act; %no act or mult act
            end

        end
        %%
        n_succ_tot = zeros(Ns,1);
        n_small_succ_tot = zeros(Ns,1);

        idx_act_C1 = find(exitC_tot(1,:)>=0.9);
        lostC1 = isempty(idx_act_C1);

        idx_act_C2 = find(exitC_tot(2,:)>=0.9);
        lostC2 = isempty(idx_act_C2);

        if lostC1==0  
            n_succ_tot(1) = length(intersect(find(Winner_tot==1),idx_act_C1));
            n_small_succ_tot(1) = length(intersect(find(Small_winner_tot==1),idx_act_C1));
            %                 n_err_tot(1) = length(setdiff(find(Winner_tot==1),idx_act_C1));
        end

        if lostC2==0   
            n_succ_tot(2) = length(intersect(find(Winner_tot==2),idx_act_C2));
            n_small_succ_tot(2) = length(intersect(find(Small_winner_tot==2),idx_act_C2));
            %             n_err_tot(2) = length(setdiff(find(Winner_tot==2),idx_act_C2));
        end


        n_succ = sum(n_succ_tot);
        n_small_succ = sum(n_small_succ_tot);
        n_succ_tot;
        n_small_succ_tot;
        Array_n_succ_tot(:,ii) = n_succ_tot;
        Array_n_small_succ_tot(:,ii) = n_small_succ_tot;
    end

    Array_n_succ = sum(Array_n_succ_tot);
    Array_n_small_succ = sum(Array_n_small_succ_tot);
    Array_error= 100 - Array_n_succ - Array_n_small_succ;

    mu_succ = mean(Array_n_succ);
    std_succ = std(Array_n_succ);
    mu_small_succ = mean(Array_n_small_succ);
    std_small_succ = std(Array_n_small_succ);
    Matrice_risultati(indice_epoca,:) = [Epoca mu_succ std_succ mu_small_succ std_small_succ];
end

% disp(Matrice_risultati)
TEST=array2table(Matrice_risultati);
TEST.Properties.VariableNames={'Epoch','Mean 1st Action','Std 1st Action','Mean 2nd Action','Std 2nd Action'};
TEST

