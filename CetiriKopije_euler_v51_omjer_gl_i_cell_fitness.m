clear;
close all
clc;
tic

% Parameters
chr = 20;
idxL = zeros(1, 4);
idxR = zeros(1, 4);
idx = [0, 0, 0];
kor = 1/(2*log(2));
p0 = (1.7e-3)*kor;%5e-3;%2e-5
pa = 0e-3; %širi
b0 = log(2)/(24/24);
b1 = 0.2*b0;
% b2 = 0.2;
a2 = 0.0;
nmax = 1e100;

% Diskretisation
T = 17; %500
dt = 0.1; %0.01

% Initialisation
Nx = chr;
Nt = round(T/dt);%1000;
omjer_gl = nan(4,4);
loss_param = nan(4,4);
gain_param = nan(4,4);
fitness2N = nan(4,4);
fitness2Nplus1 = nan(4,4);
fitness2Nminus1 = nan(4,4);


i1 = 1;
i2 = 1;
for b2 = (0.2:0.05:0.55)*b0
    for a1 = (0:0.02:0.2)

    CetiriKopije_euler_v52_for_forloops
    
    omjer_gl(i1,i2) = gain(end)/loss(end);
    loss_param(i1,i2) = loss_norm(end);
    gain_param(i1, i2) = gain_norm(end);
    
%     fitness2N(i1, i2) = (1 - ...
%         ((0 + 2*chr + 3*0 + 4*0) * p(0, chr, 0, 0)) - ...
%         (alpha(0, chr, 0, 0))) * B(0, chr, 0, 0);

%     fitness2Nplus1(i1, i2) = (1 - ...
%         ((0 + 2*(chr-1) + 3*1 + 4*0) * p(0, (chr-1), 1, 0)) - ...
%         (alpha(0, (chr-1), 1, 0))) * B(0, (chr-1), 1, 0);

    fitness2Nminus1(i1, i2) = (1 - ...
        (2*(1*1 + 2*(chr-1) + 3*0 + 4*0) * p(1, (chr-1), 0, 0)) - ...
        (2*alpha(1, (chr-1), 0, 0))) * B(1, (chr-1), 0, 0);
    
    i2 = i2+1;
    end
    i2=1;
    i1 = i1+1;
end


toc
% 
% figure()
% imagesc((0:0.1:0.6) - a2, (0.2:0.1:0.8) - b1, omjer_gl(:,:))
% colorbar
% colormap jet
% xlabel('apoptosis of monosomies')
% ylabel('proliferation of monosomies')
% set(gca,'FontSize',20)
% 
% figure()
% plot(0:0.1:0.6, omjer_gl(1,:), 'LineWidth',2, 'DisplayName','Prolif1N = 0.8')
% hold on
% % plot(0:0.1:0.6, omjer_gl(2,:), 'LineWidth',2)
% % hold on
% % plot(0:0.1:0.6, omjer_gl(3,:), 'LineWidth',2)
% % hold on
% plot(0:0.1:0.6, omjer_gl(4,:), 'LineWidth',2, 'DisplayName','Prolif1N = 0.5')
% hold on
% % plot(0:0.1:0.4, omjer_gl(5,:), 'LineWidth',2)
% % hold on
% % plot(0:0.1:0.4, omjer_gl(6,:), 'LineWidth',2)
% % hold on
% plot(0:0.1:0.6, omjer_gl(7,:), 'LineWidth',2, 'DisplayName','Prolif1N = 0.2')
% hold on
% legend('Location','southeast')
% xlabel('apoptosis of monosomies')
% ylabel('gain/loss')
% set(gca,'FontSize',20)
% 
% 
% figure()
% plot(0.2:0.1:0.8, omjer_gl(:,1), 'LineWidth',2)
% hold on
% plot(0.2:0.1:0.8, omjer_gl(:,2), 'LineWidth',2)
% hold on
% plot(0.2:0.1:0.8, omjer_gl(:,3), 'LineWidth',2)
% hold on
% plot(0.2:0.1:0.8, omjer_gl(:,4), 'LineWidth',2)
% hold on
% plot(0.2:0.1:0.8, omjer_gl(:,5), 'LineWidth',2)
% hold on
% plot(0.2:0.1:0.8, omjer_gl(:,6), 'LineWidth',2)
% hold on
% plot(0.2:0.1:0.8, omjer_gl(:,7), 'LineWidth',2)
% 
% xlabel('proliferation of monosomies')
% ylabel('gain/loss')
% set(gca,'FontSize',20)
% 
% figure()
% plot(0.2:0.1:0.8, fitness2Nminus1(:,1),'LineWidth', 2)
% hold on
% % plot(omjer_gl(:,2),fitness2Nminus1(:,2), 'o','LineWidth', 2)
% % hold on
% plot(0.2:0.1:0.8, fitness2Nminus1(:,3),'LineWidth', 2)
% hold on
% % plot(omjer_gl(:,4),fitness2Nminus1(:,4), 'o','LineWidth', 2)
% % hold on
% plot(0.2:0.1:0.8, fitness2Nminus1(:,5),'LineWidth', 2, 'DisplayName','a1=0.4')
% hold on
% % plot(omjer_gl(:,6),fitness2Nminus1(:,6), 'o','LineWidth', 2)
% % hold on
% % plot(omjer_gl(:,7),fitness2Nminus1(:,7), 'o','LineWidth', 2)
% legend('Location','northeast')
% % xlim([0.8 6.5])
% ylabel('cell fitness')
% xlabel('Prolif')
% set(gca,'FontSize',20)
% 
% figure()
% plot(loss_param(:,1), fitness2Nminus1(:,1),'LineWidth', 2)
% hold on
% % plot(omjer_gl(:,2),fitness2Nminus1(:,2), 'o','LineWidth', 2)
% % hold on
% plot(loss_param(:,3), fitness2Nminus1(:,3),'LineWidth', 2)
% hold on
% % plot(omjer_gl(:,4),fitness2Nminus1(:,4), 'o','LineWidth', 2)
% % hold on
% plot(loss_param(:,5), fitness2Nminus1(:,5),'LineWidth', 2, 'DisplayName','a1=0.4')
% hold on
% % plot(omjer_gl(:,6),fitness2Nminus1(:,6), 'o','LineWidth', 2)
% % hold on
% % plot(omjer_gl(:,7),fitness2Nminus1(:,7), 'o','LineWidth', 2)
% % legend('Location','northeast')
% % % xlim([0.8 6.5])
% % ylabel('cell fitness')
% % xlabel('loss(%)')
% % set(gca,'FontSize',20)
% 
% % figure()
% plot(gain_param(:,1), fitness2Nminus1(:,1),'LineWidth', 2, 'DisplayName','a1=0')
% hold on
% % plot(omjer_gl(:,2),fitness2Nminus1(:,2), 'o','LineWidth', 2)
% % hold on
% plot(gain_param(:,3), fitness2Nminus1(:,3),'LineWidth', 2, 'DisplayName','a1=0.2')
% hold on
% % plot(omjer_gl(:,4),fitness2Nminus1(:,4), 'o','LineWidth', 2)
% % hold on
% plot(gain_param(:,5), fitness2Nminus1(:,5),'LineWidth', 2, 'DisplayName','a1=0.4')
% hold on
% % plot(omjer_gl(:,6),fitness2Nminus1(:,6), 'o','LineWidth', 2)
% % hold on
% % plot(omjer_gl(:,7),fitness2Nminus1(:,7), 'o','LineWidth', 2)
% % legend('Location','northeast')
% % xlim([0.8 6.5])
% ylabel('cell fitness')
% xlabel('gain (%)')
% set(gca,'FontSize',20)
% 
% 
% figure()
% plot(omjer_gl(:,1),fitness2Nminus1(:,1), 'LineWidth', 2, 'DisplayName','a1=0')
% hold on
% % plot(omjer_gl(:,2),fitness2Nminus1(:,2), 'o','LineWidth', 2)
% % hold on
% plot(omjer_gl(:,3),fitness2Nminus1(:,3),'LineWidth', 2, 'DisplayName','a1=0.2')
% hold on
% % plot(omjer_gl(:,4),fitness2Nminus1(:,4), 'o','LineWidth', 2)
% % hold on
% plot(omjer_gl(:,5),fitness2Nminus1(:,5),'LineWidth', 2, 'DisplayName','a1=0.4')
% hold on
% % plot(omjer_gl(:,6),fitness2Nminus1(:,6), 'o','LineWidth', 2)
% % hold on
% % plot(omjer_gl(:,7),fitness2Nminus1(:,7), 'o','LineWidth', 2)
% legend('Location','northeast')
% xlim([0.8 6.5])
% xlabel('gain/loss')
% ylabel('cell fitness')
% set(gca,'FontSize',20)
% 
% figure()
% imagesc((0.2:0.1:0.8), 0:0.05:0.6 ,fitness2Nminus1(:,:))
% colorbar
% colormap winter
% ylabel('apoptosis of monosomies')
% xlabel('proliferation of monosomies')
% set(gca,'FontSize',20)

figure()
contourf((0:0.02:0.2), (b0-((0.2:0.05:0.55)*b0))/b0, omjer_gl(:,:), 11, 'ShowText','on')
colorbar
colormap('winter');
% ylim([0.48 0.8])
% xlim([0 0.2])
pbaspect([1 1 1])
ylabel('proliferation of monosomies')
xlabel('apoptosis of monosomies')
set(gca,'FontSize',20)


figure()
contourf((0:0.02:0.2), (b0-((0.2:0.05:0.55)*b0))/b0, fitness2Nminus1(:,:),'ShowText','on')
colorbar
colormap('winter');
% ylim([0.48 0.8])
pbaspect([1 1 1])
ylabel('proliferation of monosomies')
xlabel('apoptosis of monosomies')
set(gca,'FontSize',20)
