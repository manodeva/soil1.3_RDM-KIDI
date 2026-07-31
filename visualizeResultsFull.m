function visualizeResultsFull(modelSetting,iPore,t,y,para,bc_m)

nxy = para.nxplus*para.nyplus;
nt  = length(t);

% -------------------------------------------------
% concentration of substrates and bacteria 
% -------------------------------------------------
c = y(:,1:para.nc*nxy)/para.cmax;
b = y(:,para.nc*nxy+1:end)/para.bmax;

for ic=1:para.nc
    cc{ic}=c(:,(ic-1)*nxy+1:ic*nxy);
end

for ib=1:para.nb
    bb{ib}=b(:,(ib-1)*nxy+1:ib*nxy);
end

% -------------------------------------------------
% simulation grids
% -------------------------------------------------
x1 = linspace(0,para.Lx/1000,para.nxplus);

if modelSetting.dimension==2
    y1 = linspace(0,para.Ly/1000,para.nyplus);
    [X,Y] = meshgrid(y1,x1);
end

% -------------------------------------------------
% subplot sizes
% -------------------------------------------------
nConc = para.nc + para.nb;
ncol1 = ceil(sqrt(nConc));
nrow1 = ceil(nConc/ncol1);

nAij = para.nb*para.nb;
ncol2 = ceil(sqrt(nAij));
nrow2 = ceil(nAij/ncol2);

% =================================================
% TIME LOOP
% =================================================
for it=1:nt

%% ===============================================
% FIGURE 1 : Concentrations of substrates and microbes
% ================================================
figure(1); clf
set(gcf,'color','w','position',[50 50 1400 850])

pid=1;

% ---------------- substrates ----------------
for ic=1:para.nc

subplot(nrow1,ncol1,pid)

if modelSetting.dimension==1

    plot(x1,cc{ic}(it,:),'LineWidth',1.6)
    ylim([0 1])
    xlim([0 x1(end)])

else

    cfield = reshape(cc{ic}(it,:),para.nxplus,para.nyplus);

    surf(X,Y,cfield,'EdgeColor','none')
    view(2)
    cb = colorbar;
    cb.Label.String = 'c/cmax';
    caxis([0 1])
    axis tight equal

end

title(['Nutrient ',num2str(ic)])
%xlabel('Position [mm]')
%ylabel('c/cmax')
set(gca,'fontsize',11,'linewidth',1.1)

pid=pid+1;

end

% ---------------- microbes ----------------
for ib=1:para.nb

subplot(nrow1,ncol1,pid)

if modelSetting.dimension==1

    plot(x1,bb{ib}(it,:),'LineWidth',1.6)
    ylim([0 1])
    xlim([0 x1(end)])

else

    bfield = reshape(bb{ib}(it,:),para.nxplus,para.nyplus);

    surf(X,Y,bfield,'EdgeColor','none')
    view(2)
    cb = colorbar;
    cb.Label.String = 'b/bmax';
    caxis([0 1])
    axis tight equal

end

title(['Bacteria ',num2str(ib)])
% xlabel('Position [mm]')
% ylabel('b/bmax')
set(gca,'fontsize',11,'linewidth',1.1)

pid=pid+1;

end

sgtitle(['Concentrations, Time = ',num2str(t(it)/3600,'%.2f'),' h'], ...
'fontsize',16,'fontweight','bold')

drawnow

%% ===============================================
% FIGURE 2 : Pairwise Interactions Coefficients (KIDI)
% ================================================

figure(2); clf
set(gcf,'color','w','position',[80 80 1450 900])

pid=1;

for i=1:para.nb
for j=1:para.nb

subplot(nrow2,ncol2,pid)

% initialize
if modelSetting.dimension==1
    aij = zeros(1,para.nxplus);
else
    aij = zeros(para.nxplus,para.nyplus);
end

for ic=1:para.nc

    if modelSetting.dimension==1
        cfield = cc{ic}(it,:);
    else
        cfield = reshape(cc{ic}(it,:),para.nxplus,para.nyplus);
    end

    aij = aij + ...
        (para.gamma(ic,i)*para.cchar(ic,i)) ./(cfield + para.cchar(ic,i)).^2 .* ... % d_mu/ds
        (-para.kappa(ic,j)/para.gamma(ic,j).*para.bmax); % ds/dx (normalized) 

end

% ----- plot -----
if modelSetting.dimension==1

    plot(x1,aij,'LineWidth',1.6)
    xlim([0 x1(end)])
    ytickformat('%.3g')

else

    surf(X,Y,aij,'EdgeColor','none')
    view(2)
    colorbar
    axis tight equal

end

title(['a_{',num2str(i),num2str(j),'}'])
%xlabel('Position [mm]')
%ylabel('Interaction')
set(gca,'fontsize',11,'linewidth',1.1)

pid=pid+1;

end
end

sgtitle(['Pairwise Interaction Coefficients, Time = ', ...
num2str(t(it)/3600,'%.2f'),' h'], ...
'fontsize',16,'fontweight','bold')

drawnow

end

end
