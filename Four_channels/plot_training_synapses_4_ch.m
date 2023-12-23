%==========================================================================
%Description: The script performs the plot of all synapses subject to
% training concerning one learning procedure.
%
% To be run after ‘Basal_Training_synapses_4_channels.m’ or
% ‘Reversal_Training_synapses_4_channels.m’.
%
% The script loads simulations results that come from .mat files saved
% during training. The user can choose the simulation in line 23.
% 
% Mauro Ursino, Miriam Schirru Jan. 2022
%==========================================================================
clc
clear all
close all
% %--------------------------------------------------------------------------
% % Simulations plot 
% %--------------------------------------------------------------------------

oldpath = path;
selpath = uigetdir;
path(selpath,oldpath)
load W_tot_post_post_1.mat

width = 2;
font = 18;
N_epoche = 4*N_epoche;
t = 1:N_epoche;
%plot Wgc

figure (1)
subplot(4,4,1)
plot(t,squeeze(Wgc_epocs(1,1,1:N_epoche)),'linewidth',width)
subplot(4,4,2)
plot(t,squeeze(Wgc_epocs(1,2,1:N_epoche)),'linewidth',width)
subplot(4,4,3)
plot(t,squeeze(Wgc_epocs(1,3,1:N_epoche)),'linewidth',width)
subplot(4,4,4)
plot(t,squeeze(Wgc_epocs(1,4,1:N_epoche)),'linewidth',width)

subplot(4,4,5)
plot(t,squeeze(Wgc_epocs(2,1,1:N_epoche)),'linewidth',width)
subplot(4,4,6)
plot(t,squeeze(Wgc_epocs(2,2,1:N_epoche)),'linewidth',width)
subplot(4,4,7)
plot(t,squeeze(Wgc_epocs(2,3,1:N_epoche)),'linewidth',width)
subplot(4,4,8)
plot(t,squeeze(Wgc_epocs(2,4,1:N_epoche)),'linewidth',width)

subplot(4,4,9)
plot(t,squeeze(Wgc_epocs(3,1,1:N_epoche)),'linewidth',width)
subplot(4,4,10)
plot(t,squeeze(Wgc_epocs(3,2,1:N_epoche)),'linewidth',width)
subplot(4,4,11)
plot(t,squeeze(Wgc_epocs(3,3,1:N_epoche)),'linewidth',width)
subplot(4,4,12)
plot(t,squeeze(Wgc_epocs(3,4,1:N_epoche)),'linewidth',width)

subplot(4,4,13)
plot(t,squeeze(Wgc_epocs(4,1,1:N_epoche)),'linewidth',width)
subplot(4,4,14)
plot(t,squeeze(Wgc_epocs(4,2,1:N_epoche)),'linewidth',width)
subplot(4,4,15)
plot(t,squeeze(Wgc_epocs(4,3,1:N_epoche)),'linewidth',width)
subplot(4,4,16)
plot(t,squeeze(Wgc_epocs(4,4,1:N_epoche)),'linewidth',width)

sgtitle('W_G_C')

%%
%plot Wgs
figure (2)
subplot(4,4,1)
plot(t,squeeze(Wgs_epocs(1,1,1:N_epoche)),'linewidth',width)
subplot(4,4,2)
plot(t,squeeze(Wgs_epocs(1,2,1:N_epoche)),'linewidth',width)
subplot(4,4,3)
plot(t,squeeze(Wgs_epocs(1,3,1:N_epoche)),'linewidth',width)
subplot(4,4,4)
plot(t,squeeze(Wgs_epocs(1,4,1:N_epoche)),'linewidth',width)

subplot(4,4,5)
plot(t,squeeze(Wgs_epocs(2,1,1:N_epoche)),'linewidth',width)
subplot(4,4,6)
plot(t,squeeze(Wgs_epocs(2,2,1:N_epoche)),'linewidth',width)
subplot(4,4,7)
plot(t,squeeze(Wgs_epocs(2,3,1:N_epoche)),'linewidth',width)
subplot(4,4,8)
plot(t,squeeze(Wgs_epocs(2,4,1:N_epoche)),'linewidth',width)

subplot(4,4,9)
plot(t,squeeze(Wgs_epocs(3,1,1:N_epoche)),'linewidth',width)
subplot(4,4,10)
plot(t,squeeze(Wgs_epocs(3,2,1:N_epoche)),'linewidth',width)
subplot(4,4,11)
plot(t,squeeze(Wgs_epocs(3,3,1:N_epoche)),'linewidth',width)
subplot(4,4,12)
plot(t,squeeze(Wgs_epocs(3,4,1:N_epoche)),'linewidth',width)

subplot(4,4,13)
plot(t,squeeze(Wgs_epocs(4,1,1:N_epoche)),'linewidth',width)
subplot(4,4,14)
plot(t,squeeze(Wgs_epocs(4,2,1:N_epoche)),'linewidth',width)
subplot(4,4,15)
plot(t,squeeze(Wgs_epocs(4,3,1:N_epoche)),'linewidth',width)
subplot(4,4,16)
plot(t,squeeze(Wgs_epocs(4,4,1:N_epoche)),'linewidth',width)

sgtitle('W_G_S')
%%
%plot Wnc
figure (3)
subplot(4,4,1)
plot(t,squeeze(Wnc_epocs(1,1,1:N_epoche)),'linewidth',width)
subplot(4,4,2)
plot(t,squeeze(Wnc_epocs(1,2,1:N_epoche)),'linewidth',width)
subplot(4,4,3)
plot(t,squeeze(Wnc_epocs(1,3,1:N_epoche)),'linewidth',width)
subplot(4,4,4)
plot(t,squeeze(Wnc_epocs(1,4,1:N_epoche)),'linewidth',width)

subplot(4,4,5)
plot(t,squeeze(Wnc_epocs(2,1,1:N_epoche)),'linewidth',width)
subplot(4,4,6)
plot(t,squeeze(Wnc_epocs(2,2,1:N_epoche)),'linewidth',width)
subplot(4,4,7)
plot(t,squeeze(Wnc_epocs(2,3,1:N_epoche)),'linewidth',width)
subplot(4,4,8)
plot(t,squeeze(Wnc_epocs(2,4,1:N_epoche)),'linewidth',width)

subplot(4,4,9)
plot(t,squeeze(Wnc_epocs(3,1,1:N_epoche)),'linewidth',width)
subplot(4,4,10)
plot(t,squeeze(Wnc_epocs(3,2,1:N_epoche)),'linewidth',width)
subplot(4,4,11)
plot(t,squeeze(Wnc_epocs(3,3,1:N_epoche)),'linewidth',width)
subplot(4,4,12)
plot(t,squeeze(Wnc_epocs(3,4,1:N_epoche)),'linewidth',width)

subplot(4,4,13)
plot(t,squeeze(Wnc_epocs(4,1,1:N_epoche)),'linewidth',width)
subplot(4,4,14)
plot(t,squeeze(Wnc_epocs(4,2,1:N_epoche)),'linewidth',width)
subplot(4,4,15)
plot(t,squeeze(Wnc_epocs(4,3,1:N_epoche)),'linewidth',width)
subplot(4,4,16)
plot(t,squeeze(Wnc_epocs(4,4,1:N_epoche)),'linewidth',width)

sgtitle('W_N_C')
%%
%plot Wns
figure (4)
subplot(4,4,1)
plot(t,squeeze(Wns_epocs(1,1,1:N_epoche)),'linewidth',width)
subplot(4,4,2)
plot(t,squeeze(Wns_epocs(1,2,1:N_epoche)),'linewidth',width)
subplot(4,4,3)
plot(t,squeeze(Wns_epocs(1,3,1:N_epoche)),'linewidth',width)
subplot(4,4,4)
plot(t,squeeze(Wns_epocs(1,4,1:N_epoche)),'linewidth',width)

subplot(4,4,5)
plot(t,squeeze(Wns_epocs(2,1,1:N_epoche)),'linewidth',width)
subplot(4,4,6)
plot(t,squeeze(Wns_epocs(2,2,1:N_epoche)),'linewidth',width)
subplot(4,4,7)
plot(t,squeeze(Wns_epocs(2,3,1:N_epoche)),'linewidth',width)
subplot(4,4,8)
plot(t,squeeze(Wns_epocs(2,4,1:N_epoche)),'linewidth',width)

subplot(4,4,9)
plot(t,squeeze(Wns_epocs(3,1,1:N_epoche)),'linewidth',width)
subplot(4,4,10)
plot(t,squeeze(Wns_epocs(3,2,1:N_epoche)),'linewidth',width)
subplot(4,4,11)
plot(t,squeeze(Wns_epocs(3,3,1:N_epoche)),'linewidth',width)
subplot(4,4,12)
plot(t,squeeze(Wns_epocs(3,4,1:N_epoche)),'linewidth',width)

subplot(4,4,13)
plot(t,squeeze(Wns_epocs(4,1,1:N_epoche)),'linewidth',width)
subplot(4,4,14)
plot(t,squeeze(Wns_epocs(4,2,1:N_epoche)),'linewidth',width)
subplot(4,4,15)
plot(t,squeeze(Wns_epocs(4,3,1:N_epoche)),'linewidth',width)
subplot(4,4,16)
plot(t,squeeze(Wns_epocs(4,4,1:N_epoche)),'linewidth',width)

sgtitle('W_N_S')


reward_tot = sum(vett_reward)
punishment_tot = sum(vett_punishment)
no_answer_tot = sum(vett_no_risposta)
small_rew_tot = sum(vett_small_winn)


    
%  % I compute the cumulative sum of all responses
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
% % title('Hebb rule post-pre','fontsize',font)
% legend('no response','punishment','reward','location','northwest')
% set(gca,'fontsize',font)

