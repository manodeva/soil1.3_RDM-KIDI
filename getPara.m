% this is the function used for specifying the input parameters of spatial simulations

function para=getPara(modelSetting,iPore)

if modelSetting.dimension==1
    para.Lx=10; % mm
	para.Ly=0; % mm
elseif modelSetting.dimension==2
    para.Lx=10; % mm
    para.Ly=10; % mm
else
    error('Check your model setting!')    
end

% - the number of bacteria
if modelSetting.noOfBacteria==1
    para.nb=1;
elseif modelSetting.noOfBacteria==2
    para.nb=2; 
elseif modelSetting.noOfBacteria==3
    para.nb=3; 
elseif modelSetting.noOfBacteria==4
    para.nb=4;
else
    error('Check your model setting!')
end

% - the number of chemicals (nutrients)
if modelSetting.noOfChemicals==1
    para.nc=1; 
elseif modelSetting.noOfChemicals==2
    para.nc=2; 
elseif modelSetting.noOfChemicals==3
    para.nc=3; 
elseif modelSetting.noOfChemicals==4
    para.nc=4; 
else
    error('Check your model setting!')
end

% - maximal concentrations of b and c
% para.bmax=1e12; % <---- theoretical maximal number density of cells
para.bmax=1e10; 
para.cmax=15;

% - length and grid size
para.h=100; % um
%para.h=1; % mm

para.Lx=para.Lx*1000; % mm -> um
para.nx=round(para.Lx/para.h); % the number of grids
para.nxplus=para.nx+1;        

para.Ly=para.Ly*1000; % mm -> um
para.ny=round(para.Ly/para.h); % the number of grids
para.nyplus=para.ny+1;

nxyplus=para.nxplus*para.nyplus;

% - porosity-dependent parameters 
if iPore==1
    para.epsilon=0.36; % porosity
    para.Db=2.32*ones(para.nb,1); % um^2/s
    para.l_c=4.6; % um
%     para.chi0=1e4*ones(para.nc,para.nb); % chemotatic parameter (guessed)
    if para.nb == 1
    para.chi0=5e3; % chemotatic parameter (guessed)
    else
        para.chi0 = ones(para.nc,para.nb)*5e3/4;
        para.chi0(1) = 0;
    end
elseif iPore==2
    para.epsilon=0.17; % porosity
    para.Db=0.93*ones(para.nb,1); % um^2/s
    para.l_c=3.1; % um
	para.chi0=1e4*0.8*ones(para.nc,para.nb); % chemotatic parameter (guessed)
elseif iPore==3
    para.epsilon=0.04; % porosity
    para.Db=0.42*ones(para.nb,1); % um^2/s
    para.l_c=2.4; % um
	para.chi0=1e4*0.8^2*ones(para.nc,para.nb); % chemotatic parameter (guessed)
end

if strcmp(modelSetting.chemotaxis,'off')
    para.chi0=zeros(para.nc,para.nb);
end

% - time span depending on initial cell position and speed
if strcmp(modelSetting.initialDist_Bacteria,'center') 
    if modelSetting.dimension==1
%         para.tspan=0:0.01:0.2; % hours  
        para.tspan=0:0.01:0.8; % hours  
    elseif modelSetting.dimension==2
        para.tspan=0:0.01:0.5; % hours  
    end    
elseif strcmp(modelSetting.initialDist_Bacteria,'end')        
    if iPore==1
%         para.tspan=0:0.01:0.5; % hours
        para.tspan=0:0.01:1; % hours
    elseif iPore==2
        para.tspan=0:0.01:1; % hours 
    elseif iPore==3
        para.tspan=0:0.01:1; % hours  
    end
end

% - initial distribution of chemicals (nutrients)
if modelSetting.dimension==1
    if strcmp(modelSetting.initialDist_Chemicals,'uniform')
    c0=[];
    for ic=1:para.nc
        c0=[c0;para.cmax*ones(para.nxplus,1)]; % uniform distribution
    end
    para.c0=c0;    
    else
    error('Check your model setting!')
    end
elseif modelSetting.dimension==2
    if strcmp(modelSetting.initialDist_Chemicals,'uniform')
    c0=[];
    for ic=1:para.nc
        c0_=para.cmax*ones(para.nxplus,para.nyplus);
        c0=[c0;c0_(:)];
    end
    para.c0=c0;
    elseif strcmp(modelSetting.initialDist_Chemicals,'patch')
        c0 = zeros(para.nxplus,para.nyplus,para.nc);
        np = 5; % change this value to alter no of patches 
        np_x = round(linspace(1,para.nxplus,np));
        np_y = round(linspace(1,para.nyplus,np));
        
        for ic = 1:para.nc
        for ip_x = round((np/2)-1):2:np
            for ip_y = round((np/2)-1):2:np
                c0(np_x(ip_x)-2:np_x(ip_x)+2,np_y(ip_y)-2:np_y(ip_y)+2,ic) = para.cmax;
            end
        end
        c0(np_x(round(np/2))-2:np_x(round(np/2))+2,np_x(round(np/2))-2:np_x(round(np/2))+2,ic) = para.cmax;
        end
   
        
       c0 = permute(c0,[2 1 3]);
       para.c0=c0(:);
    end
else
    error('Check your model setting!')
end

% - initial distribution of bacteria 
if modelSetting.dimension==1
    if strcmp(modelSetting.initialDist_Bacteria,'end')
        para.b0=zeros(para.nb*para.nxplus,1);
        if para.nb == 1
        para.b0(1:3)=para.bmax;
        else
            for i = 0:para.nb-1
            idx = i*para.nxplus + (1:3);   % first three grid points of each block
            para.b0(idx) = para.bmax;
            end
        end
    elseif strcmp(modelSetting.initialDist_Bacteria,'center')
        para.b0=zeros(para.nb*para.nxplus,1);
        para.nxhalf=round(para.nx/2)+1;
        if para.nb==1
            para.b0(para.nxhalf-2:1:para.nxhalf+2)=para.bmax;
        else
            for ib = 1:para.nb
                if ib == 1
                para.b0(para.nxhalf-2:1:para.nxhalf+2)=para.bmax;
                else     
                para.b0((nxyplus*(ib-1))+para.nxhalf-2:1:(nxyplus*(ib-1))+para.nxhalf+2)=para.bmax;
                end
            end
        end
        
    else
        error('Check your model setting!')
    end
    
elseif modelSetting.dimension==2
    if strcmp(modelSetting.initialDist_Bacteria,'end')
        if para.nb == 1
        b0=zeros(para.nxplus,para.nyplus);
        b0(1:para.nxplus,[para.nyplus]) = para.bmax;
        b0([para.nxplus],1:para.nyplus) = para.bmax;
        para.b0=b0(:);
        else
        b0=zeros(para.nxplus,para.nyplus,para.nb); 
        for ib = 1:para.nb
        b0(:,1:3,ib) = para.bmax;
        end
        b0 = permute(b0,[2 1 3]);
        para.b0=b0(:);
        end
    elseif strcmp(modelSetting.initialDist_Bacteria,'center')
        if para.nb == 1
        b0=zeros(para.nxplus,para.nyplus);
        para.nxhalf=round(para.nx/2)+1;
        para.nyhalf=round(para.ny/2)+1;
        b0(para.nxhalf,para.nyhalf)=para.bmax;
        para.b0=b0(:);
        else
        b0=zeros(para.nxplus,para.nyplus,para.nb);
        para.nxhalf=round(para.nx/2)+1;
        para.nyhalf=round(para.ny/2)+1;
        for ib = 1:para.nb
        b0(para.nxhalf,para.nyhalf,ib)=para.bmax;
        end
        b0 = permute(b0,[2 1 3]);
        para.b0=b0(:);
        end
   else
        error('Check your model setting!')
    end
end

% elseif modelSetting.dimension==2
%     if strcmp(modelSetting.initialDist,'center')
%         b0=zeros(para.nxplus,para.nyplus);
%         para.nxhalf=round(para.nx/2)+1;
%         para.nyhalf=round(para.ny/2)+1;
% %         b0(para.nxhalf-2:1:para.nxhalf+2,para.nyhalf-2:1:para.nyhalf+2)=para.bmax;
%         b0(para.nxhalf,para.nyhalf)=para.bmax;
% %         b0(para.nxhalf-5:1:para.nxhalf+5,para.nyhalf-2:1:para.nyhalf+2)=para.bmax;
%         para.b0=b0(:);
%     else
%         error('Check your model setting!')
%     end
% end

% - other parameters (setting-agnostic)
para.tspan=para.tspan*3600; % hours to seconds 

% Diffusivity of substrate
if strcmp(modelSetting.initialDist_Chemicals,'uniform')
para.Dc=800*ones(para.nc,1); % um^2/s
elseif strcmp(modelSetting.initialDist_Chemicals,'patch')
para.Dc=800*ones(para.nc,1); % um^2/s    
end

para.kappa=1.6e-11*ones(para.nc,para.nb); % Max uptake rate; mM/[(cell/mL)s]=uM/[(cell/L)s]
para.cchar=1*ones(para.nc,para.nb); % half saturation constant; uM 
para.cminus=1*ones(para.nc,para.nb); % lower bound of concentration sensitivity for chemotaxis; uM 
para.cplus=30*ones(para.nc,para.nb); % upper bound of concentration sensitivity for chemotaxis; uM 
para.gamma=log(2)/3600*ones(para.nc,para.nb); % max growth rate; 1/hr to 1/sec 
para.d=1*ones(para.nb,1); % bacterial length; um 
para.kd = para.gamma.*0.02; % Death rate; 1/sec
