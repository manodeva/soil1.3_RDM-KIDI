% code for simulating growth of microbial community using reaction-diffusion model
% The code can simulate upto of community of 4, substrates of upto 4, and perform 1D and 2D simulations 
% we assumed no flux boundary conditions for substrates.  
% for bacteria you can use different boundary conditions such as no flux, periodic and constant flux boundary condition
% In the simulation we have substrates are uniformly distributed or in patches while bacteria are assumed to be in centre or at the end of simulated region


function runRxnDiff
 clc
 clear
 close all

% - choose your parameter set
modelSetting.chemotaxis='on'; % 'on' or 'off'
modelSetting.visualization='on'; % 'on' or 'off'
modelSetting.video='off'; % 'on' or 'off' (only when visualization is on)
%---
modelSetting.noOfBacteria=3; % 1 to 4
modelSetting.noOfChemicals=1; % 1 to 4
modelSetting.dimension=1; % 1 or 2
modelSetting.initialDist_Bacteria='center'; % 'end' or 'center'
modelSetting.initialDist_Chemicals='uniform'; % 'uniform' or 'patch' (patch considered to be 5x5 matrix at 5 locations. Use patch for 2D simulations only) 
bc_m='no_flux'; % boundary condition for microbes; 'no_flux' or 'periodic' or 'const_flux'

%----------------------------

% - start simulation 

for iPore=1
    para=getPara(modelSetting,iPore);
    y0=[para.c0;para.b0];
    % Assumed 5% of max flux for inlet and outlet fluxes
    if(strcmp(bc_m,'const_flux'))
    for ib = 1:para.nb    
    max_flux(ib) = (para.bmax/para.h)*(para.Db(ib)+para.chi0(ib)*log(1+(para.c0(1)/para.cminus(1))/1+(para.c0(1)/para.cplus(1))));
    end
    para.bflux_in = max_flux.*ones(para.nc,para.nb).*0.05; 
    para.bflux_out = max_flux.*ones(para.nc,para.nb).*0.05;
    end
    odeOptions = odeset('NonNegative', 1:length(y0));
    if modelSetting.dimension==1
        [t,y] = ode23(@myOde,para.tspan,y0,odeOptions,para,bc_m);
    elseif modelSetting.dimension==2
        [t,y] = ode23(@myOde,para.tspan,y0,odeOptions,para,bc_m);
    end
    
    % - visualize c and b profiles 
    if strcmp(modelSetting.visualization,'on')
        visualizeResultsFull(modelSetting,iPore,t,y,para,bc_m)
    end

end

