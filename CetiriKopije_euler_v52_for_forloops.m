
nold(1:Nx+3, 1:Nx+3, 1:Nx+3,1) = 0.0;
nnew(1:Nx+3, 1:Nx+3, 1:Nx+3,1) = 0.0;
n(1:Nx+3, 1:Nx+3, 1:Nx+3, 1:Nt+1) = 0.0;
fitness(1:Nx+3, 1:Nx+3, 1:Nx+3) = 0;
time = linspace(0, T, Nt+1);
apoptoza = zeros(1, Nt+1);
apoptoza_test = zeros(1, Nt+1);
prolif = zeros(1, Nt+1);
gain = zeros(1, Nt+1);
gain_norm = zeros(1, Nt+1);
gain_norm_aneu = zeros(1, Nt+1);
loss = zeros(1, Nt+1);
loss_norm = zeros(1, Nt+1);
loss_norm_aneu = zeros(1, Nt+1);
chrom = zeros((4*chr)-(chr-1), Nt);
ntot_all = zeros((4*chr)-(chr-1),1);
ntot = zeros((4*chr)-(chr-1), Nt);
ntot2 = zeros((4*chr)-(chr-1), Nt);

t = 0; m = 1; 

% Initial condition & Boundary condition %DIPLOIDNI PU!!
nold(:, :, :, 1) = 0;
nold(2, 2, 2, 1) = 1e0;


%%% Funcije parametri

p = @(x1, x2, x3, x4) p0 + pa * (x2 ~= chr);                  

pwgd = @(x1, x2, x3, x4) 0.0;

%  B = @(x1, x2, x3, x4) b0 - b1 * (1 - (x2/chr)) * ((x3 > 0 || x4 > 0) && x1 == 0 && x2 ~= chr) - b2 * (1 - (x2/chr)) * (x1 > 0 && x2 ~= chr);                     
B = @(x1, x2, x3, x4) b0 - b1 * ((x3 > 0 || x4 > 0) && x1 == 0 && x2 ~= chr) - b2 * (x1 > 0 && x2 ~= chr);
                  
alpha = @(x1, x2, x3, x4) (x1>0)*a1;
% alpha = @(x1, x2, x3, x4) (x1>0)*0.25;
% alpha = @(x1, x2, x3, x4) 0.05 + a1 * (1 - (x2/chr)) * (x1>0 && (x3 == 0 && x4 == 0) && x2 ~= chr) + a2 * (1 - (x2/chr)) * ((x3 > 0 || x4 > 0) && x1 == 0 && x2 ~= chr);
% alpha = @(x1, x2, x3, x4) 0.05 + a1 * (x1>0 && (x3 == 0 && x4 == 0) && x2 ~= chr) + a2 * ((x3 > 0 || x4 > 0) && x1 == 0 && x2 ~= chr);

                  
%  gain = 0;
%  loss = 0;
%  apoptoza = 0;
%  apoptoza_test = 0;
%  prolif = 0;
%  
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
                    
%                     fitness(idx(1), idx(3), idx(4)) = 1 + (1 - ...
%                         ((x1 + 2*x2 + 3*x3 + 4*x4) * p(x1, x2, x3, x4)) - ...
%                         (alpha(x1, x2, x3, x4))) * log(2);

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
                
            end
        end
    end
    
    gain_norm_aneu(m) = 100*(gain(m)/((sum(nold(:,:,:),'all')-nold(2,2,2))*chr));
    
    loss_norm_aneu(m) = 100*(loss(m)/((sum(nold(:,:,:),'all')-nold(2,2,2))*chr));

    gain_norm(m) = 100*(gain(m)/((sum(nold(:,:,:),'all'))*chr));
    
    loss_norm(m) = 100*(loss(m)/((sum(nold(:,:,:),'all'))*chr));

    
    t = t + dt;
    
    if t >= time(m)
        n(:, :, :, m) = nold(:, :, :);
        ntot(:, m) = ntot_all(:);
        m = m + 1;
    end
    
    nold = nnew;
    ntot_all(:) = 0;   
     
end

