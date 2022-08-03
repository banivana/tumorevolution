clear;
tic
% RATE Euler

%%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%beta([1,18,0,0],3)

chr=19;
kor = 1/(2*log(2));

p0 =@(j14) (1.7e-3)*kor;
% p0 =@(j14)  ((j14<=2)*(1.7e-3)*kor) + ((j14>2)*(3e-3)*kor);   %**** stanice s 3 kopije chr14imaju povecan p
% p0  = @(j14) (0.0017 - (j14>2)*0.0016)*kor;

%apop = @(x,j15,j14) (x(1)*(x(1)>0)+(j15<2)+(j14<2))*0.06; %+ (x(3)>0)*0.13;
apop = @(x,j15,j14) ((x(1)>0)||(j15<2)||(j14<2))*0.2;
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
br_st = zeros(points+1,1);
br_aneupl_st= zeros(points+1,1);

an=0;
bn=0;
while ii <= points+1

    for i1 = 2:chr
        for i2 = 2:(chr-i1+2)
            for i3 = 2:(chr-i1-i2+4)
                x = [i1-2,i2-2,i3-2,chr-i3-i2-i1+4];
                for ji15 = 2:5
                    j15=ji15 - 1;
                    for ji14 = 2:5
                        j14= ji14 - 1;
                        
                        
                        chrom_tot = sum((1:4).*x) + j15 + j14;
                        
                        F = (1 - 2*(sum((1:4).*x)+j15 +j14)*p0(j14) -2*apop(x,j15,j14) )*beta(x,j15,j14)*nold(i1,i2,i3,ji15,ji14) ...
                            +( (x(1)+1)*beta([x(1)+1,x(2)-1,x(3),x(4)],j15,j14)*nold(i1+1,i2-1,i3,ji15,ji14) ...
                            + 2*(x(2)+1)* ( beta([x(1)-1,x(2)+1,x(3),x(4)],j15,j14)*nold(i1-1,i2+1,i3,ji15,ji14) + beta([x(1),x(2)+1,x(3)-1,x(4)],j15,j14)*nold(i1,i2+1,i3-1,ji15,ji14) )...
                            + 3*(x(3)+1)* ( beta([x(1),x(2)-1,x(3)+1,x(4)],j15,j14)*nold(i1,i2-1,i3+1,ji15,ji14) + beta([x(1),x(2),x(3)+1,x(4)-1],j15,j14)*nold(i1,i2,i3+1,ji15,ji14) )...
                            + 4*(x(4)+1)* beta([x(1),x(2),x(3)-1,x(4)+1],j15,j14)*nold(i1,i2,i3-1,ji15,ji14)...
                            + (j15-1)*beta(x,j15-1,j14)*nold(i1,i2,i3,ji15-1,ji14) + (j15+1)*beta(x,j15+1,j14)*nold(i1,i2,i3,ji15+1,ji14))*p0(j14)...
                            + (j14-1)*beta(x,j15,j14-1)*nold(i1,i2,i3,ji15,ji14-1)*p0(j14-1) + (j14+1)*beta(x,j15,j14+1)*nold(i1,i2,i3,ji15,ji14+1)*p0(j14+1);
                        
                        nnew(i1,i2,i3,ji15,ji14) = nold(i1,i2,i3,ji15,ji14) + dt * F;
                        
                        ntot_everyt(chrom_tot-chr+1) = ntot_everyt(chrom_tot-chr+1) + nold(i1,i2,i3,ji15,ji14);
                        ntot_everyt_j15(chrom_tot-chr+1,j15) = ntot_everyt_j15(chrom_tot-chr+1,j15) + nold(i1,i2,i3,ji15,ji14);
                        
                        loss(ii) = loss(ii) + nnew(i1,i2,i3,ji15,ji14)*x(1);
                        gain(ii) = gain(ii) + nnew(i1,i2,i3,ji15,ji14)*(x(3) + x(4));
                        loss_j15(ii) = loss_j15(ii) +nnew(i1,i2,i3,ji15,ji14)*(j15==1);
                        gain_j15(ii) = gain_j15(ii) +nnew(i1,i2,i3,ji15,ji14)*((j15==3)+(j15==4));
                        loss_j14(ii) = loss_j14(ii) +nnew(i1,i2,i3,ji15,ji14)*(j14==1);
                        gain_j14(ii) = gain_j14(ii) +nnew(i1,i2,i3,ji15,ji14)*((j14==3)+(j14==4));
                        br_st(ii) = br_st(ii) + nnew(i1,i2,i3,ji15,ji14);
                        
                        
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
        br_aneupl_st(ii) = br_st(ii) - nold(2,chr,2,3,3);
        ii = ii + 1;
    end
    nold = nnew;
    ntot_everyt(:) = 0;
    ntot_everyt_j15(:) = 0;
    
end
toc

tren = T;

%%
ii=points+1;

figure();
bar([1:13,16:chr],ones(1,chr-2)*(100*gain(ii)/(chr-2)/br_st(ii)),'FaceColor','r')
hold on
bar([1:13,16:chr],-ones(1,chr-2)*(100*loss(ii)/(chr-2)/br_st(ii)),'FaceColor', 'g')
hold on
bar(15,100*gain_j15(ii)/br_st(ii),'FaceColor','r')
hold on
bar(15,-100*loss_j15(ii)/br_st(ii),'FaceColor','g')
hold on
bar(14,100*gain_j14(ii)/br_st(ii),'FaceColor','r')
hold on
bar(14,-100*loss_j14(ii)/br_st(ii),'FaceColor','g')
hold on
xlabel('Chromosome')
ylabel('% Chr. Gains/losses')
% legend('Gain','Loss')
xlim([0 20])
xticks(1:1:19)
xtickangle(45)
ylim([-10 100])
yticks(0:20:100)
pbaspect([1 1 1])
% daspect([0.15 1.5 1])
set(gca,'FontSize',19)



%%
%%%%% g(t) of chr15 and chr14
% M1 = 0.0859; stderror1 = 0.0252;
% M2 = 1.1507; stderror2 = 0.2744;
% M3 = 0.9127; stderror3 = 0.1713;
% M4 = 2.554; stderror4 = 0.4723;
% M5 = 20.6625; stderror5 = 0.8775;

M1 = 0.0859; 
M2 = 1.1507; 
M3 = 0.9127; 
M4 = 2.554;
M4_chr1415 = 6.2368; semerr4_chr1415 = 2.34;
M4_others = 2.1223; semerr4_others = 0.3971;
M5 = 20.6625; semerr5_chr1415 = 1.07; %srednja vr chr14 i 15 ; 
M5_other = 4.2028; semerr5_oth = 0.7420;
MTL = 88.2; semerrTL = 4.7;
MTL_others = 8.0624; semerrTL_others = 2.5573;
% M5 = 21.75; %chr14
% M6 = 19.5; %chr15


figure();
plot(0:dt:T,100*gain_j15./br_st, 'LineWidth',2,'Color','r');
hold on
plot(0:dt:T,100*gain_j14./br_st, 'LineWidth',2,'Color',[1 0 1]);
hold on
plot(0:dt:T,100*gain./((chr-2) * br_st), 'LineWidth',2,'Color','g');
hold on
% errorbar([5, 12, 15, 17, 19],[M1,M2,M3,M4,M5], [stderror1, stderror2, stderror3, stderror4, stderror5],'o','LineWidth',1.5, 'Color','k')
plot([5, 12, 15, 19, 19],[M1,M2,M3,M5,M5_other],'o','LineWidth',2,'MarkerFaceColor', 'k', 'Color','k')
hold on
errorbar([39,39], [MTL, MTL_others], [semerrTL, semerrTL_others], 'o','LineWidth',2,'MarkerFaceColor', 'k', 'Color','k')
hold on
errorbar([17,17], [M4_chr1415, M4_others], [ semerr4_chr1415, semerr4_others], 'o','LineWidth',2,'MarkerFaceColor', 'k', 'Color','k')
pbaspect([1 1 1])
xlabel('Time')
ylabel('Gain (%)')
% legend({'Gain of chr15', 'Gain of chr14'},'Location','northwest');
set(gca,'FontSize',19)

%%
%%%%% BR aneuploidnih stan u t
% figure();
% plot(0:dt:T,br_aneupl_st./br_st*100, 'LineWidth',2,'Color','b');
% xlabel('Time')
% ylabel('Aneuploid cells %')
% set(gca,'FontSize',19)


%%

% %gain/loss
% 
% figure();
% hold on;
% bar([1:14,16:chr],ones(1,chr-1)*(100*gain/(chr-1)/br_st2),'FaceColor',[0.8 0.0 0.1])
% bar([1:14,16:chr],-ones(1,chr-1)*(100*loss/(chr-1)/br_st2),'FaceColor',[0 .7 .4])
% bar(15,100*gain_j/br_st2,'FaceColor',[0.8 0.0 0.1])
% bar(15,-100*loss_j/br_st2,'FaceColor',[0 .7 .4])
% xlabel('Chromosome')
% ylabel('% Chr. Gains/losses')
% legend('Gain','Loss','','')
% xticks(1:20)
% set(gca,'FontSize',15)
% 
% br_st=sum(ntot(:, time == tren));
% br_st2 = sum(nold(:));
% figure();
% hold on;
% bar(1,(gain/(chr-1)),'FaceColor',[0.8 0.0 0.1])
% bar(1,-(loss/(chr-1)),'FaceColor',[0 .7 .4])
% bar(2,gain_j,'FaceColor',[0.8 0.0 0.1])
% bar(2,-loss_j,'FaceColor',[0 .7 .4])
% xlabel('Chromosome')
% ylabel('Chr. Gains/losses')




% function a = apop(x,j15,j14)
%     a=0.0;
%     if( x(3)~=0 || j15~=2 || j14~=2 || x(4)~=0 || x(1)~=0)
%         a=0.2;
%     end
%     if(j14>2)
%        a=0.1;
%     end
%     if ( x(1)~=0 || j15<2 || j14<2) 
%        a=0.2;   %najvise 0.29
%     end
% 
% end


% function a = apop(x,j15,j14)
%     a=0.0;
%     if( x(3)~=0 || j15~=2 || j14~=2 || x(4)~=0 || x(1)~=0)
%         a=0.2;
%     end
%     if(j15>2)
%        a=0;
%     end
%     if ( x(1)~=0 || j15<2 || j14<2) 
%        a=0.2;   %najvise 0.29
%     end
% 
% end

function  b = beta(x,j15,j14)
b=1;
%%%% Ivana
    if( x(3)~=0 || j15>2 || j14>2 || x(4)~=0 )
        b =0.8;
    end
    if( x(1)~=0 || j15<2 || j14<2)
        b = 0.8;
    end
%%%%%% chr14 i 15
% if(j15>2)    
%     b=1.8;
% end
% if(j14>2)
%     b=1.2;
% end
if (j15>2 && j14>2)
    b=1.8;
end


b=b*log(2);
end
