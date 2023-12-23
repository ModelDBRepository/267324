%==========================================================================
%Description: The script performs the plot of all synapses subject to
% training concerning one learning procedure.
%
% To be run after ‘Basal_Training_synapses_2_channels.m’ or
% ‘Reversal_Training_synapses_2_channels.m’.
%
% The script loads simulations results that come from .mat files saved
% during training. The user can choose the simulation in line 23.
% 
% Mauro Ursino, Miriam Schirru Jan. 2022
%==========================================================================
clear all
close all
clc
% %--------------------------------------------------------------------------
% % Simulations plot 
% %--------------------------------------------------------------------------
oldpath = path;
selpath = uigetdir;
path(selpath,oldpath)

load W_tot_post_adapt_1

N_epoche = N_epochs*2;
t = 1:1:N_epoche;

width = 1.5;
font = 18;
font1 = 12;

t = t/2;

%plot Wgc
figure
subplot(2,2,1)
plot(t,squeeze(Wgc_epocs(1,1,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,2)
plot(t,squeeze(Wgc_epocs(1,2,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,3)
plot(t,squeeze(Wgc_epocs(2,1,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,4)
plot(t,squeeze(Wgc_epocs(2,2,1:N_epoche)),'linewidth',width)
sgtitle('W^G^C','fontsize',font)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)

%plot Wgs
figure
subplot(2,2,1)
plot(t,squeeze(Wgs_epocs(1,1,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,2)
plot(t,squeeze(Wgs_epocs(1,2,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,3)
plot(t,squeeze(Wgs_epocs(2,1,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,4)
plot(t,squeeze(Wgs_epocs(2,2,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
sgtitle('W^G^S','fontsize',font)

%plot Wnc
figure
subplot(2,2,1)
plot(t,squeeze(Wnc_epocs(1,1,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,2)
plot(t,squeeze(Wnc_epocs(1,2,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,3)
plot(t,squeeze(Wnc_epocs(2,1,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,4)
plot(t,squeeze(Wnc_epocs(2,2,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
sgtitle('W^N^C','fontsize',font)

%plot Wns
figure
subplot(2,2,1)
plot(t,squeeze(Wns_epocs(1,1,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,2)
plot(t,squeeze(Wns_epocs(1,2,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,3)
plot(t,squeeze(Wns_epocs(2,1,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
subplot(2,2,4)
plot(t,squeeze(Wns_epocs(2,2,1:N_epoche)),'linewidth',width)
xlabel('Epoch')
ylabel('synapse strength')
set(gca,'fontsize',font1)
sgtitle('W^N^S', 'fontsize',font)


% 
% reward_tot = sum(vett_reward)
% punishment_tot = sum(vett_punishment)
% no_answer_tot = sum(vett_no_risposta)
% small_rew_tot = sum(vett_small_winn)
% 
% 
%  for kk = 1:N_epoche
% tot_vett_no_risposta(kk) = sum(vett_no_risposta(1:kk));
% tot_vett_punishment(kk) = sum(vett_punishment(1:kk));
% tot_vett_reward(kk) = sum(vett_reward(1:kk));
% tot_vett_small_winn(kk) = sum (vett_small_winn(1:kk));
%  end
%  index = (1:N_epoche);
% figure
% plot(index,tot_vett_no_risposta,'y',index,tot_vett_punishment,'r',index,tot_vett_reward,'g','linewidth',width)
% hold on
% plot(index,tot_vett_small_winn,'b','linewidth',width)
% xlabel('Number of epochs','fontsize',font)
% ylabel('cumulative distribution','fontsize',font)
% title('yellow: no response; red: punishment; green: reward','fontsize',font)
% set(gca,'fontsize',font)

    
    