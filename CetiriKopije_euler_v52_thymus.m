clear;
% close all
clc;
tic

% Parameters

chr = 19;
idxL = zeros(1, 4);
idxR = zeros(1, 4);
idx = [0, 0, 0];
kor = 1/(2*log(2));
p0 = (1.7e-3)*kor;%5e-3;%2e-5
pa = 0e-3; %širi
b0 = log(2)/(24/24);
b1 = 0.2*b0;
b2 = 0.2*b0;%0.8
a1 = 0.197;%0.4;
a2 = 0.0;
nmax = 1e100;



% Diskretisation
T = 17; %500
dt = 0.1; %0.01

% Initialisation
Nx = chr;
Nt = round(T/dt);%1000;
nold(1:Nx+3, 1:Nx+3, 1:Nx+3,1) = 0.0;
nnew(1:Nx+3, 1:Nx+3, 1:Nx+3,1) = 0.0;
n(1:Nx+3, 1:Nx+3, 1:Nx+3, 1:Nt+1) = 0.0;
fitness(1:Nx+3, 1:Nx+3, 1:Nx+3) = 0;
time = linspace(0, T, Nt+1);
t = 0; m = 1; 

chrom = zeros((4*chr)-(chr-1), Nt);
% chrom = zeros((4*copies)+1, Nt);
ntot_all = zeros((4*chr)-(chr-1),1);
ntot = zeros((4*chr)-(chr-1), Nt);
ntot2 = zeros((4*chr)-(chr-1), Nt);

% Initial condition & Boundary condition %DIPLOIDNI PU!!
% nold(:, :, :, 1) = 0;
% nold(chr+2, 2, 2, 1) = 1;
nold(:, :, :, 1) = 0;
nold(2, 2, 2, 1) = 1e0;
apoptoza = zeros(1, Nt+1);
apoptoza_test = zeros(1, Nt+1);
prolif = zeros(1, Nt+1);
gain = zeros(1, Nt+1);
gain_norm = zeros(1, Nt+1);
loss = zeros(1, Nt+1);
loss_norm = zeros(1, Nt+1);
br_st_aneu = zeros(1, Nt+1);


%%% Funcije parametri

p = @(x1, x2, x3, x4) p0 + pa * (x2 ~= chr);                  

pwgd = @(x1, x2, x3, x4) 0.0;% + pwda * (3.164 * (1 - (x1/chr)) * (1 - ((x2)/chr)) * ...
                      %(1 - (x3/chr)) * (1 - (x4/chr)));

            
B = @(x1, x2, x3, x4) b0 - b1 * ((x3 > 0 || x4 > 0) && x1 == 0 && x2 ~= chr) - b2 * (x1 > 0 && x2 ~= chr);

alpha = @(x1, x2, x3, x4) (x1>0)*a1;

                  

while m <= Nt+1  %Time loop
   
     if sum(nold(:,:,:),'all') <= nmax
        saturation = (1 - (sum(nold(:,:,:),'all')/nmax));
    else
        saturation = 0;
    end
    
    for i = 2 : Nx+2
        idx(1) = i; x(1) = i-2;
        for j = 2 : ((Nx+4)-i)
            idx(3) = j; x(3) = j-2;
            for k = 2 : ((Nx+6)-i-j)
                idx(4) = k; x(4) = k-2;
                x(2) = chr - x(1) - x(3) - x(4); 
                idx(2) = chr + (4 * 2) - idx(1) - idx(3) -idx(4);
             
                x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
                
                n0gd = nold(idx(1), idx(3), idx(4))* B(x1, x2, x3, x4)* pwgd(x1, x2, x3, x4);
                n0a = nold(idx(1), idx(3), idx(4)) * alpha(x1, x2, x3, x4) * B(x1, x2, x3, x4);
                n0 = nold(idx(1), idx(3), idx(4))* B(x1, x2, x3, x4);
                n0p = nold(idx(1), idx(3), idx(4)) * p(x1, x2, x3, x4) * B(x1, x2, x3, x4);
                n1 = nold(idx(1)+1, idx(3), idx(4)) * p(x1 + 1, x2 - 1, x3, x4) * B(x1 + 1, x2 - 1, x3, x4);
                n4 = nold(idx(1), idx(3)-1, idx(4)+1) * p(x1, x2, x3 - 1, x4 + 1) * B(x1, x2, x3 - 1, x4 + 1);
%                 gainTrm = gain2E2(x, idx, nold) * saturation;

                apoptoza(m) = apoptoza(m) + n0a;
                apoptoza_test(m) = apoptoza_test(m) + nold(idx(1), idx(3), idx(4));
                prolif(m) = prolif(m) + n0;

                if m == Nt+1
                    apoptoza = apoptoza + n0a;
                    apoptoza_test = apoptoza_test + nold(idx(1), idx(3), idx(4));
                    prolif = prolif + n0; 
                    
                    fitness(idx(1), idx(3), idx(4)) = (1 - ...
                        ((x1 + 2*x2 + 3*x3 + 4*x4) * p(x1, x2, x3, x4)) - ...
                        (alpha(x1, x2, x3, x4))) * B(x1, x2, x3, x4);
                end
                
                
                gainTrm = 0;
                for ncopy = 2:3
                    iL=idx; iR=idx;
                    iL(ncopy) = idx(ncopy) + 1; iL(ncopy-1) = iL(ncopy-1) - 1;
                    iR(ncopy) = idx(ncopy) + 1; iR(ncopy+1) = iR(ncopy+1) - 1;
                    
                    gainTrm = gainTrm + ncopy * (x(ncopy) + 1) * (nold(iL(1),iL(3),iL(4)) * p(iL(1)-2,iL(2)-2,iL(3)-2,iL(4)-2) * B(iL(1)-2,iL(2)-2,iL(3)-2,iL(4)-2)...
                                + nold(iR(1),iR(3),iR(4)) * p(iR(1)-2,iR(2)-2,iR(3)-2,iR(4)-2) * B(iR(1)-2,iR(2)-2,iR(3)-2,iR(4)-2));
                end

                if x1 == 0 && x3 ==0
                    ngd = nold(idx(2), 2, 2) * pwgd(x2, (chr-x2), 0, 0) * B(x2, (chr-x2), 0, 0) * saturation;%ngd = n(idx(2), idx(4), 2, m);
                else
                    ngd = 0;
                end
                
                chrom_tot = 1 * x(1) + 2 * x(2) + 3 * x(3) + 4 * x(4);

                
                F = ((1  * n0) - 2 * n0gd  - 2 * (x1 + 2*x2 + 3*x3 + 4*x4) * n0p + ...
                    + (x1+1) * n1 + gainTrm + 4 * (x4 + 1) * n4 + ngd) * saturation - 2 * n0a;
                
                loss(m) = loss(m) + nold(idx(1), idx(3), idx(4))*x1;
                gain(m) = gain(m) + nold(idx(1), idx(3), idx(4))*(x3 + x4);

                nnew(idx(1), idx(3), idx(4)) = nold(idx(1), idx(3), idx(4)) + dt * F;
                ntot_all(chrom_tot-(chr-1)) = ntot_all(chrom_tot-(chr-1)) + nold(idx(1), idx(3), idx(4));
%                 ntot2(chrom_tot-22, m) = ntot2(chrom_tot-22, m) + n(idx(1), idx(2), idx(3), m+1);
                
%                 if m == Nt+1
%                     loss = loss + nnew(idx(1), idx(3), idx(4))*x1;
%                     gain = gain + nnew(idx(1), idx(3), idx(4))*(x3 + x4);
%                 end

                    
               
            end
        end
    end
    
%     gain_norm(m) = 100*(gain(m)/((sum(nold(:,:,:),'all')-nold(2,2,2))*chr));
%     
%     loss_norm(m) = 100*(loss(m)/((sum(nold(:,:,:),'all')-nold(2,2,2))*chr));

    gain_norm(m) = 100*(gain(m)/((sum(nold(:,:,:),'all'))*chr));
    
    loss_norm(m) = 100*(loss(m)/((sum(nold(:,:,:),'all'))*chr));
    
    br_st_aneu(m) = sum(nold(:,:,:),'all')-nold(2,2,2);
    
    t = t + dt;
    
    if t >= time(m)
        n(:, :, :, m) = nold(:, :, :);
        ntot(:, m) = ntot_all(:);
        m = m + 1;
    end
    
    nold = nnew;
    ntot_all(:) = 0;   
     
end

% omjer = apoptoza / prolif;

toc
 
attime = T;
maxnb = max(n(:,:,:,time == attime),[],'all');


figure()
bar(1:(chr),100*ones(1,(chr))*(gain(end)/((sum(nold(:,:,:),'all')-nold(2,2,2))*chr)), 'r')
hold on 
bar(1:(chr),-100*(loss(end)/(chr*(sum(nold(:,:,:),'all')-nold(2,2,2))))*ones(1,(chr)), 'g')
ylabel('Chr Gains or Losses')
xlabel('Chromosomes 1-19')
xlim([0 20])
xticks(1:1:19)
xtickangle(45)
ylim([-10 20])
yticks(-10:10:20)
daspect([0.3 1 1])
%  daspect([0.15 1.5 1])
% pbaspect([2 0.7 1])
set(gca,'FontSize',20)

% exp = 9.15/2.73
theory = gain(end)/loss(end)
theory2 = gain_norm(end)/loss_norm(end)

gain_exp = load('ThymusDN4_gain.txt');
gain_exp = gain_exp * (20/152); %(1pixel = 20/152%)
loss_exp = load('ThymusDN4_loss.txt');
loss_exp = loss_exp * (20/152); %(1pixel = 20/152%)
exp_GL = sum(gain_exp)/sum(loss_exp)
m_gain = sum(gain_exp)/19
m_loss = sum(loss_exp)/19


figure()
bar(1:(19),gain_exp, 'r')
hold on 
bar(1:(19),-loss_exp, 'g')
ylabel('Chr. Gains or Losses (%)')
xlabel('Chromosomes 1-19')
xlim([0 20])
xticks(1:1:19)
xtickangle(45)
ylim([-20 20])
yticks(-20:10:20)
% daspect([0.3 1 1])
daspect([0.41 1.5 1])
% pbaspect([2 0.7 1])
set(gca,'FontSize',20)

figure()
bar(1:19, m_gain * ones(1,chr), 'r')
hold on
bar(1:19, -m_loss * ones(1,chr), 'g')
xlim([0 20])
xticks(1:1:19)
xtickangle(45)
ylim([-20 20])
yticks(-20:10:20)
% daspect([0.3 1 1])
daspect([0.41 1.5 1])
% pbaspect([2 0.7 1])
set(gca,'FontSize',20)

figure()
plot(time, gain_norm, 'r')
hold on 
plot(time, loss_norm, 'g')
hold off

% 
% figure()
% 
% % M = load('ThymicLymphomas.txt'); M = M * 10/44; % *10/44
% M = load('Tymus.txt');  M = M * 100/33; % *(100/33)
% 
% % M = load('Intestine.txt'); M = M * 1000/30; % *1000/30
% 
% % Tymus
% bar(38:45,M/sum(M), 'b')
% hold on
% % bar(2*chr,(n(2,2,2, Nt+1)*max(M))/(sum(M)*max(ntot(:, Nt+1))), 'k')
% % hold on
% xlim([chr 4*chr])
% xticks([chr 2*chr 3*chr 4*chr])
% set(gca,'FontSize',18)
% 
% % plot(chr:(4*chr),((ntot(:, Nt+1)*max(M))/(sum(M)*max(ntot(:, Nt+1)))), 'r', 'LineWidth',2)
% plot(chr:(4*chr),ntot(:, Nt+1)/sum(ntot(:, Nt+1)), 'r', 'LineWidth',2)
% hold on
% % bar(chr:(4*chr),((ntot(:, Nt+1)*max(M))/max(ntot(:, Nt+1))), 'r')
% % hold on
% xlabel('Number of chromosomes')
% ylabel('Cell ratio')
% xticks([chr 2*chr 3*chr 4*chr])
% set(gca,'FontSize',18)
% hold on

% figure()
% plot(time, br_st_aneu, 'LineWidth',2)
% xlabel('Time')
% ylabel('Broj aneu stanica')
% set(gca,'FontSize',18)

%%% T-Lymphomas
% bar(40:47,M)
% xlim([chr 4*chr])
% xticks([chr 2*chr 3*chr 4*chr])
% set(gca,'FontSize',15)



% figure(7)
% M = load('Tymus.txt');
% bar(38:45,M*(100/33))
% xlim([chr 4*chr])
% xticks([chr 2*chr 3*chr 4*chr])
% set(gca,'FontSize',15)
% 
% figure(8)
% M = load('Intestine.txt');
% bar(37:44,M*(1000/30))
% xlim([chr 4*chr])
% xticks([chr 2*chr 3*chr 4*chr])
% set(gca,'FontSize',15)
% figure(5)
% bar(gain/(chr))
% hold on 
% bar(-loss/(chr))
% ylabel('% Chr. Gains/Losses')
% xlabel('Chromosomes 1-20')
% %xticks(1:1:20)
% set(gca,'FontSize',15)

% 
% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'Intestine_June_cellRatio.gif';
% for n = 1:50
%     % Draw plot for y = x.^n
%     plot(chr:(4*chr),((ntot(:,round(n*(Nt+1)/50))*max(M))/(sum(M)*max(ntot(:, round(n*(Nt+1)/50))))), 'r', 'LineWidth',2)
%     xlabel('Number of chromosomes')
%     ylabel('Cell ratio')
%     xticks([chr 2*chr 3*chr 4*chr])
%     set(gcf,'color','white');
%     set(gca,'FontSize',18)
%     drawnow 
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.05); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.05); 
%       end 
% end

% function b = B(x1, x2, x3, x4)
%     b = 1;
%     if (x3 > 0 || x4 > 0)
%        b=0.8;
%     end
%     if (x1 > 0)
%        b=0.8;
%     end
% end


% B = @(x1, x2, x3, x4) b0 - 0.2 * ((x1 + 2*x2 + 3*x3 + 4*x4)<(2*chr)) - 0.2 * ((x1 + 2*x2 + 3*x3 + 4*x4)>(2*chr));

