clear;
tic
% RATE Euler
%%%%%%%%%%%%% 
%%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%beta([1,18,0,0],3)

chr=19;

kor = 1/(2*log(2));
%kor = 1;
%p_mis = 1.38*0.05*xlsread('pmis.xlsx','D3:D25')';
%p_mis= 0.05/15*(-9/44*(1:23)+6.704545); %%perc. of aneupl 
%p_mis= 0.05/100*(-9/44*(1:23)+6.704545); %%perc. of aneupl 
%p_mis=0.002*ones(23,1);


% p0 = 0.001;
% p15 = 0.001;
% p14 = 0.003;


% p0  = @(j14) 0.001 + (j14>2)*0.002;
% p15 = @(j14) 0.001 + (j14>2)*0.002;
% p14 = @(j14) 0.001 + (j14>2)*0.002;

p0  = @(j14) (0.0017 - (j14>2)*0.0016)*kor;
p15 = @(j14) (0.0017 - (j14>2)*0.0016)*kor;
p14 = @(j14) (0.0017 - (j14>2)*0.0016)*kor;

% p0  = @(j14) 0.005 - (j14>2)*0.0025;
% p15 = @(j14) 0.005 - (j14>2)*0.0025;
% p14 = @(j14) 0.005 - (j14>2)*0.0025;

% pp = 5;
% p0  = @(j14) (j14<=2)*0.01 + (j14>2)*0.001*pp;
% p15 = @(j14) (j14<=2)*0.01 + (j14>2)*0.001*pp;
% p14 = @(j14) (j14<=2)*0.01 + (j14>2)*0.001*pp;

% p0  = @(j14) 0.001 + (j14<2)*0.01;
% p15 = @(j14) 0.001 + (j14<2)*0.01;
% p14 = @(j14) 0.001 + (j14<2)*0.01;



% p0 = mean(p_mis(2:22));
% p15 = p_mis(1);
% p14 = p_mis(23);
%apop = @(x,j15,j14) (x(1)*(x(1)>0)+(j15<2)+(j14<2))*0.06; %+ (x(3)>0)*0.13;
apop = @(x,j15,j14) ((x(1)>0)||(j15<2)||(j14<2))*0.2*kor;
%apop = @(x,j15,j14) 0;

T = 50;
dt = 0.1;
points = round(T/dt);
%nmax = 1e100;   % nmax za saturaciju
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0;
time = linspace(0, T, points+1);   %trenutci t u kojima se spremju vrijednosti n i ntot
Nx=chr-2;

%%% Initialisation
%n = zeros(chr+1,chr+1,chr+1,6,6, points); 
nnew = zeros(chr+1,chr+1,chr+1,6,6);
ntot_everyt = zeros((4*chr)-(chr-1),1);
ntot_everyt_j15 = zeros((4*chr)-(chr-1),4); 
ntot = zeros((4*chr)-(chr-1), points);
ntot_j15 = zeros((4*chr)-(chr-1),4, points);

%%%% INITIAL CONDITIONS
nold = zeros(chr+1,chr+1,chr+1,6,6);
% %%Diploidna stanica
nold(2,chr,2,3,3) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = 1;
loss=zeros(points+1,1);
gain=zeros(points+1,1);
loss_j15 = zeros(points+1,1);
gain_j15 = zeros(points+1,1);
loss_j14 = zeros(points+1,1);
gain_j14 = zeros(points+1,1);
loss_j14i15 = zeros(points+1,1);
gain_j14i15 = zeros(points+1,1);
br_st = zeros(points+1,1);
an=0;
bn=0;
while ii <= points+1
%     if sum(nold(:)) <= nmax    %sum(nold(:),'all')
%         saturation = (1 - (sum(nold(:))/nmax));   %saturation = (1 - (sum(nold(:),'all')/nmax));
%     else
%         saturation = 0;
%     end
%     
    for i1 = 2:chr
        for i2 = 2:(chr-i1+2)
            for i3 = 2:(chr-i1-i2+4)
                x = [i1-2,i2-2,i3-2,chr-i3-i2-i1+4];
                for ji15 = 2:5
                    j15=ji15 - 1;
                    for ji14 = 2:5
                        j14= ji14 - 1;
                        
                        %                     if ji15==4
                        %                         fnj1 = fn(i1,i2,i3,x,ji15-1,nold);
                        %                         fnj2=0;
                        %                     elseif ji15==1
                        %                         fnj1 = 0;
                        %                         fnj2= fn(i1,i2,i3,x,ji15+1,nold);
                        %                     else
                        %                         fnj1 = fn(i1,i2,i3,x,ji15-1,nold);
                        %                         fnj2= fn(i1,i2,i3,x,ji15+1,nold);
                        %                     end
                        
                        chrom_tot = sum((1:4).*x) + j15 + j14;
                        
                        F = (1 - 2*(sum((1:4).*x)*p0(j14)+j15*p15(j14) +j14*p14(j14))-2*apop(x,j15,j14) )*beta(x,j15,j14)*nold(i1,i2,i3,ji15,ji14) ...
                            +(( (x(1)+1)*beta([x(1)+1,x(2)-1,x(3),x(4)],j15,j14)*nold(i1+1,i2-1,i3,ji15,ji14) ...
                            + 2*(x(2)+1)* ( beta([x(1)-1,x(2)+1,x(3),x(4)],j15,j14)*nold(i1-1,i2+1,i3,ji15,ji14) + beta([x(1),x(2)+1,x(3)-1,x(4)],j15,j14)*nold(i1,i2+1,i3-1,ji15,ji14) )...
                            + 3*(x(3)+1)* ( beta([x(1),x(2)-1,x(3)+1,x(4)],j15,j14)*nold(i1,i2-1,i3+1,ji15,ji14) + beta([x(1),x(2),x(3)+1,x(4)-1],j15,j14)*nold(i1,i2,i3+1,ji15,ji14) )...
                            + 4*(x(4)+1)* beta([x(1),x(2),x(3)-1,x(4)+1],j15,j14)*nold(i1,i2,i3-1,ji15,ji14))*p0(j14)...
                            +  ((j15-1)*beta(x,j15-1,j14)*nold(i1,i2,i3,ji15-1,ji14) + (j15+1)*beta(x,j15+1,j14)*nold(i1,i2,i3,ji15+1,ji14))*p15(j14)...
                            + ((j14-1)*beta(x,j15,j14-1)*nold(i1,i2,i3,ji15,ji14-1)*p14(j14-1) + (j14+1)*beta(x,j15,j14+1)*nold(i1,i2,i3,ji15,ji14+1)*p14(j14+1) ) );
                        
                        nnew(i1,i2,i3,ji15,ji14) = nold(i1,i2,i3,ji15,ji14) + dt * F;
                        
                        ntot_everyt(chrom_tot-chr+1) = ntot_everyt(chrom_tot-chr+1) + nold(i1,i2,i3,ji15,ji14);
                        ntot_everyt_j15(chrom_tot-chr+1,j15) = ntot_everyt_j15(chrom_tot-chr+1,j15) + nold(i1,i2,i3,ji15,ji14);
                        
                        
                        loss(ii) = loss(ii) + nnew(i1,i2,i3,ji15,ji14)*x(1);
                        gain(ii) = gain(ii) + nnew(i1,i2,i3,ji15,ji14)*(x(3) + x(4));
                        loss_j15(ii) = loss_j15(ii) +nnew(i1,i2,i3,ji15,ji14)*(j15==1);
                        gain_j15(ii) = gain_j15(ii) +nnew(i1,i2,i3,ji15,ji14)*((j15==3)+(j15==4));
                        loss_j14(ii) = loss_j14(ii) +nnew(i1,i2,i3,ji15,ji14)*(j14==1);
                        gain_j14(ii) = gain_j14(ii) +nnew(i1,i2,i3,ji15,ji14)*((j14==3)+(j14==4));
                        loss_j14i15(ii) = loss_j14i15(ii) + nnew(i1,i2,i3,ji15,ji14)*((j14<2)*(j15<2));
                        gain_j14i15(ii) = gain_j14i15(ii) + nnew(i1,i2,i3,ji15,ji14)*((j14>2)*(j15>2));
                        br_st(ii) = br_st(ii) + nnew(i1,i2,i3,ji15,ji14);
                        
%                         an= an + (apop(x,ji15)*nold(i1,i2,i3,ji15));
%                         bn= bn + (beta(x,ji15)*nold(i1,i2,i3,ji15));
                        
                    end
                end
            end
        end
    end
    
    t = t + dt;
    
    if t >= time(ii)    %spremanje za plotanje
        %n(:,:,:,:,:,ii) = nold;
        ntot(:,ii) = ntot_everyt;
        ntot_j15(:,:,ii) = ntot_everyt_j15;
        ii = ii + 1;
    end
    nold = nnew;
    ntot_everyt(:) = 0;
    ntot_everyt_j15(:) = 0;
    %an/bn
    %an=0;
    %bn=0;
    
end
toc

%an/bn

% beta([0,19,0,0],2)
% beta([0,19,0,0],3)
% beta([0,19,0,0],1)
% beta([1,18,0,0],2)
% beta([0,18,1,0],2)

tren = T;
%%
% figure()
% M = load('ThymicLymphomas.txt')*(10/44) ;   %Tumor
% M = M/sum(M);
% bar(40:47,M)
% xlim([chr 4*chr])
% xticks([chr 2*chr 3*chr 4*chr])
% set(gca,'FontSize',15)
% hold on
% 
% % figure()
% % M = load('Tymus.txt')*(100/33);
% % bar(38:45,M)
% % xlim([chr 4*chr])
% % xticks([chr 2*chr 3*chr 4*chr])
% % set(gca,'FontSize',15)
% % 
% %  
% % figure()
% % M = load('Intestine.txt')*(1000/30);
% % bar(37:44,M)
% % xlim([chr 4*chr])
% % xticks([chr 2*chr 3*chr 4*chr])
% % set(gca,'FontSize',15)
% % hold on
% 
% plot(chr:(4*chr),((ntot(:,time == tren)*max(M))/max(ntot(:, time == tren)) /sum(M) ), 'r', 'LineWidth',2)
% xlabel('Chromosome number')
% ylabel('Cell fraction')
% xticks([chr 2*chr 3*chr 4*chr])
% set(gca,'FontSize',15)
% legend('Experimental data','Theory')
% 
% % bar((2*chr),((n(2,chr+1,2,2,time == tren)*max(M))/max(ntot(:, time == tren))), 'r', 'LineWidth',2)
% % xlabel('Number of chromosomes')
% % ylabel('Cell number')
% % xticks([chr 2*chr 3*chr 4*chr])
% % set(gca,'FontSize',15)
%%
% figure();
% subplot(1,2,1)
% bar(chr:(4*chr),ntot(:, time == tren))
% xlabel('Chromosome number')
% ylabel('Number of cells')
% set(gca,'FontSize',18)
% 
% subplot(1,2,2)
% plot(chr:(4*chr),ntot_j15(:,1, time == tren) , 'LineWidth',2)
% xlabel('Tot # of chrom')
% ylabel('# of cells')
% hold on
% plot(chr:(4*chr),ntot_j15(:,2, time == tren) , 'LineWidth',2)
% plot(chr:(4*chr),ntot_j15(:,3, time == tren) , 'LineWidth',2)
% plot(chr:(4*chr),ntot_j15(:,4, time == tren) , 'LineWidth',2)
% hold off
% legend('j=1','j=2','j=3','j=4')
% set(gcf, 'Position',  [20, 20, 1500, 600])
% set(gca,'FontSize',15)
%%
%%%%% G/L u T
%br_st=sum(ntot(:, time == tren));

%br_st2 = n(2,chr+1,2,2,time == tren);

ii=points+1;

figure();
hold on;
bar([1:13,16:chr],ones(1,chr-2)*(100*gain(ii)/(chr-2)/br_st(ii)),'FaceColor','r')
bar([1:13,16:chr],-ones(1,chr-2)*(100*loss(ii)/(chr-2)/br_st(ii)),'FaceColor', 'g')
bar(15,100*gain_j15(ii)/br_st(ii),'FaceColor','r')
bar(15,-100*loss_j15(ii)/br_st(ii),'FaceColor','g')
bar(14,100*gain_j14(ii)/br_st(ii),'FaceColor','r')
bar(14,-100*loss_j14(ii)/br_st(ii),'FaceColor','g')
xlabel('Chromosome')
ylabel('% Chr. Gains/losses')
legend('Gain','Loss','','')
xticks(1:chr)
set(gca,'FontSize',19)
%%
%%%%% BR stan u t
% figure();
% plot(0:dt:T,br_st, 'LineWidth',2,'Color','b');
% xlabel('Time')
% ylabel('number of cells')

%%
%%%%%  u t test log2
% figure();
% plot(0:dt*2*log(2):T*2*log(2),100*gain_j15./br_st, 'LineWidth',2,'Color','r');
% % hold on
% % plot(0:dt:T,100*gain_j14./br_st, 'LineWidth',2,'Color',[1 0 1]);
% % hold off
% xlabel('Time')
% ylabel('Gain (%) TEST')
% legend({'Gain of chr15'},'Location','northwest');
% set(gca,'FontSize',19)

%%
%%%%% g(t) of chr15 and chr14

figure();
plot(0:dt:T,100*gain_j15./br_st, 'LineWidth',2,'Color','r');
hold on
plot(0:dt:T,100*gain_j14./br_st, 'LineWidth',2,'Color',[1 0 1]);
plot(0:dt:T,100*gain_j14i15./br_st, 'LineWidth',2,'Color','g');
plot(0:dt:T,-100*loss_j14./br_st, 'LineWidth',2);
hold off
xlabel('Time')
ylabel('Gain (%)')
legend({'Gain of chr15', 'Gain of chr14', 'Gain of chr14&15', 'loss of chr14'},'Location','northwest');
set(gca,'FontSize',19)


%%
%save('/home/npavin/Lucija/Matlab/saves/apop/d1.mat') 
%load('/home/npavin/Lucija/Matlab/saves/apop/r2.mat')
%n(2,2,2,4,time==tren)

function  b = beta(x,j15,j14)
b=1;
%%%%% Ivana
    if( x(3)~=0 || j15>2 || j14>2 || x(4)~=0 )
        b =0.8;
    end
    if( x(1)~=0 || j15<2 || j14<2)
        b = 0.8;
    end
%%%%%% chr14 i 15
% if(j15>2 && j14>2)
% %if(j15>2)    
%     b=1.8;
%     %         if(x(2)==19)
%     %             b=1.7;
%     %         else
%     %             b=1.7;
%     %         end
% %     if( x(1)~=0 )   % stanice s monosomijom i viskom j
% %         b=1.1;
% %     end
% end
%%%%%% chr 14 i 15 nezavisni
% if(j14>2)
%     b=1.8;
% end
% if(j15>2)
%     b=1.8;
% end
% 

% if(j15>2)    
%     b=1.5;
%     %         if(x(2)==19)
%     %             b=1.7;
%     %         else
%     %             b=1.7;
%     %         end
% %     if( x(1)~=0 )   % stanice s monosomijom i viskom j
% %         b=1.1;
% %     end
% end

if(j15>2)
    b=1.25;
end

% if (j15>2 && j14>2)
%     b=1.6;
% end

 b=b*log(2);

end
