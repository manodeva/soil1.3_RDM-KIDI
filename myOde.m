% This is the function where the PDE of reaction-diffusion model is
% discretized using finite difference

function [dy,correction]=myOde(t,y,para,bc_m)

%-- model complexity
nb=para.nb; % # of bacteria
nc=para.nc; % # of chemicals (nutrients)

%-- geometric parameters
h=para.h; % grid size
nxplus=para.nxplus;
nyplus=para.nyplus;
nxyplus=nxplus*nyplus;
epsilon=para.epsilon;
l_c=para.l_c;

%-- intrinsic parameters 
Dc=para.Dc;
kappa=para.kappa;
cchar=para.cchar;
cminus=para.cminus;
cplus=para.cplus;
gamma=para.gamma;
chi0=para.chi0;
Db=para.Db;
d=para.d;
kd=para.kd;

if(strcmp(bc_m,'const_flux'))
bflux_in = para.bflux_in;
bflux_out = para.bflux_out;
end
%-- initialize output variables
for ic=1:nc 
    dcc{ic}=zeros(nxplus,nyplus);
end
for ib=1:nb 
    dbb{ib}=zeros(nxplus,nyplus);
end

%-- reallocate input variables
c=y(1:nc*nxyplus);
b=y(nc*nxyplus+1:end);
for ic=1:nc
    cc{ic}=reshape(c((ic-1)*nxyplus+1:ic*nxyplus),nxplus,nyplus);
end
for ib=1:nb
    bb{ib}=reshape(b((ib-1)*nxyplus+1:ib*nxyplus),nxplus,nyplus);
end

% FORMULATE EQUATIONS =====================================================
%-- multiple chemicals
%---- mass balances for cc{ic} 
for ic=1:nc 
    
    for i=1:nxplus
        for j=1:nyplus
            switch i
                case 1 
                    c_left = cc{ic}(1,j);
                    c_right = cc{ic}(i+1,j);
                case nxplus
                    c_left = cc{ic}(i-1,j);
                    c_right = cc{ic}(nxplus,j);
                otherwise
                    c_left = cc{ic}(i-1,j);
                    c_right = cc{ic}(i+1,j);
            end
            
            cdiff = Dc(ic)*(c_right-2*cc{ic}(i,j)+c_left)/h^2;

            % add additional terms along y-axis if nd > 1
            if nyplus>1
                switch j
                    case 1 
                        c_down = cc{ic}(i,1);
                        c_up = cc{ic}(i,j+1);
                    case nyplus
                        c_down = cc{ic}(i,j-1);
                        c_up = cc{ic}(i,nyplus);
                    otherwise
                        c_down = cc{ic}(i,j-1);
                        c_up = cc{ic}(i,j+1);
                end
                % update cdiff
                cdiff = cdiff + Dc(ic)*(c_down-2*cc{ic}(i,j)+c_up)/h^2;
            end
            
            csink=0;
            for ib=1:nb
                gc=cc{ic}(i,j)/(cc{ic}(i,j)+cchar(ic,ib));
                csink = csink + bb{ib}(i,j)*kappa(ic,ib)*gc;
            end
            dcc{ic}(i,j) = cdiff - csink;

        end % j=1:nyplus
    end % i=1:nxplus
        
end

%-- multiple bacteria
%---- correction factor
% btot=zeros(size(bb{1}); % total cell concentration
btot=zeros(nxplus,nxplus);
for ib=1:nb
    btot = btot + bb{ib};
end

correction=zeros(nxplus,nyplus); % initialize correction factor

for i=1:nxplus 
    for j=1:nyplus 
        
        l_cell=(3*epsilon/(4*pi*btot(i,j)*1e-12))^(1/3)-d(ib);
        if l_cell<=0
            correction(i,j)=0;
        else
            if l_cell<=l_c
                correction(i,j)=(l_cell/l_c)^2;
            else
                correction(i,j)=1;
            end
        end
        
    end
end

%---- mass balances for bb{ib}
for ib=1:nb

    for i=1:nxplus 
        for j=1:nyplus 
            switch bc_m
                case 'no_flux'
                    switch i
                        case 1 
                            c_down = cc{ic}(i+1,j);
                            c_up = c_down;
                            b_down = bb{ib}(i+1,j);
                            b_up = b_down;
                        case nxplus
                            c_up = cc{ic}(i-1,j);
                            c_down = c_up;
                            b_up = bb{ib}(i-1,j);
                            b_down = b_up;
                        otherwise
                            c_up = cc{ic}(i-1,j);
                            c_down = cc{ic}(i+1,j);
                            b_up = bb{ib}(i-1,j);
                            b_down = bb{ib}(i+1,j);            
                    end

                case 'periodic'
                    switch i
                        case 1 
                            c_up = cc{ic}(nxplus,j);
                            c_down = cc{ic}(i+1,j);
                            b_up = bb{ib}(nxplus,j);
                            b_down = bb{ib}(i+1,j);
                        case nxplus
                            c_up = cc{ic}(i-1,j);
                            c_down = cc{ic}(1,j);
                            b_up = bb{ib}(i-1,j);
                            b_down = bb{ib}(1,j);
                        otherwise
                            c_up = cc{ic}(i-1,j);
                            c_down = cc{ic}(i+1,j);
                            b_up = bb{ib}(i-1,j);
                            b_down = bb{ib}(i+1,j);            
                    end

                case 'const_flux'

                    switch i
                        case 1 
                            c_down = cc{ic}(i+1,j);
                            c_up = c_down;
                            b_down = bb{ib}(i+1,j);
                            chemo_sens = log10((1+(cc{ic}(i,j)/cminus(ic,ib)))/(1+(cc{ic}(i,j)/cplus(ic,ib))));
                            b_up = b_down-((bflux_in(ib)*para.h)/(Db(ib)+chi0(ib)*chemo_sens));
                        case nxplus
                            c_up = cc{ic}(i-1,j);
                            c_down = c_up;
                            b_up = bb{ib}(i-1,j);
                            chemo_sens = log10((1+(cc{ic}(i,j)/cminus(ic,ib)))/(1+(cc{ic}(i,j)/cplus(ic,ib))));
                            b_down = b_up-((bflux_in(ib)*para.h)/(Db(ib)+chi0(ib)*chemo_sens));

                        otherwise
                            c_up = cc{ic}(i-1,j);
                            c_down = cc{ic}(i+1,j);
                            b_up = bb{ib}(i-1,j);
                            b_down = bb{ib}(i+1,j);            
                    end
                    

            end % bc_m.row
              
            bdiff = Db(ib)*correction(i,j)*(b_down-2*bb{ib}(i,j)+b_up)/h^2;
            
            bchemo=0;
            for ic=1:nc
                logfc_up = log10((1+c_up/cminus(ic,ib))/(1+c_up/cplus(ic,ib)));
                logfc = log10((1+cc{ic}(i,j)/cminus(ic,ib))/(1+cc{ic}(i,j)/cplus(ic,ib)));
                logfc_down = log10((1+c_down/cminus(ic,ib))/(1+c_down/cplus(ic,ib)));
                bchemo_a = (b_down-b_up)/(2*h)*(logfc_down-logfc_up)/(2*h);
                bchemo_b = bb{ib}(i,j)*(logfc_down-2*logfc+logfc_up)/h^2;
                bchemo = bchemo + chi0(ib)*correction(i,j)*(bchemo_a+bchemo_b);
            end
            
            % add additional terms along y-axis if nd > 1
            if nyplus>1
                switch bc_m
                case 'no_flux'
                    switch j
                    case 1 
                        c_right = cc{ic}(i,j+1);
                        c_left = c_right;
                        b_right = bb{ib}(i,j+1);
                        b_left = b_right;
                    case nyplus
                        c_left = cc{ic}(i,j-1);
                        c_right = c_left;
                        b_left = bb{ib}(i,j-1);
                        b_right = b_left;
                    otherwise
                        c_left = cc{ic}(i,j-1);
                        c_right = cc{ic}(i,j+1);
                        b_left = bb{ib}(i,j-1);
                        b_right = bb{ib}(i,j+1);            
                end
                
               case 'periodic'
                   switch j
                    case 1 
                        c_right = cc{ic}(i,j+1);
                        c_left = cc{ic}(i,nyplus);
                        b_right = bb{ib}(i,j+1);
                        b_left = bb{ib}(i,nyplus);
                    case nyplus
                        c_left = cc{ic}(i,j-1);
                        c_right = cc{ic}(i,1);
                        b_left = bb{ib}(i,j-1);
                        b_right = bb{ib}(i,1);
                    otherwise
                        c_left = cc{ic}(i,j-1);
                        c_right = cc{ic}(i,j+1);
                        b_left = bb{ib}(i,j-1);
                        b_right = bb{ib}(i,j+1);            
                    end

                case 'const_flux'
                    
                            switch j
                                case 1 
                                    c_right = cc{ic}(i,j+1);
                                    c_left = c_right;
                                    b_right = bb{ib}(i,j+1);
                                    chemo_sens = log10((1+(cc{ic}(i,j)/cminus(ic,ib)))/(1+(cc{ic}(i,j)/cplus(ic,ib))));
                                    b_left = b_right+((bflux_in(ib)*para.h)/(Db(ib)+chi0(ib)*chemo_sens));
                                case nyplus
                                    c_left = cc{ic}(i,j-1);
                                    c_right = c_left;
                                    b_left = bb{ib}(i,j-1);
                                    chemo_sens = log10((1+(cc{ic}(i,j)/cminus(ic,ib)))/(1+(cc{ic}(i,j)/cplus(ic,ib))));
                                    b_right = b_left-((bflux_out(ib)*para.h)/(Db(ib)+chi0(ib)*chemo_sens));
                                otherwise
                                    c_left = cc{ic}(i,j-1);
                                    c_right = cc{ic}(i,j+1);
                                    b_left = bb{ib}(i,j-1);
                                    b_right = bb{ib}(i,j+1);            
                            end
                        
                 
                end
                
                % update cdiff and bdiff
                
                bdiff = bdiff+Db(ib)*correction(i,j)*(b_left-2*bb{ib}(i,j)+b_right)/h^2;
                
                % update bchemo
                for ic=1:nc
                    logfc_left = log10((1+c_left/cminus(ic,ib))/(1+c_left/cplus(ic,ib)));
                    logfc = log10((1+cc{ic}(i,j)/cminus(ic,ib))/(1+cc{ic}(i,j)/cplus(ic,ib)));
                    logfc_right = log10((1+c_right/cminus(ic,ib))/(1+c_right/cplus(ic,ib)));
                    bchemo_a = (b_right-b_left)/(2*h)*(logfc_right-logfc_left)/(2*h);
                    bchemo_b = bb{ib}(i,j)*(logfc_right-2*logfc+logfc_left)/h^2;
                    bchemo = bchemo + chi0(ib)*correction(i,j)*(bchemo_a+bchemo_b);
                end
  
            end
            
            bsource=0;
            for ic=1:nc
                gc=cc{ic}(i,j)/(cc{ic}(i,j)+cchar(ic,ib));
                bsource = bsource + bb{ib}(i,j)*gamma(ic,ib)*gc;
            end
           
            bdeath = 0;
            bdeath = bb{ib}(i,j)*kd(ib);

            dbb{ib}(i,j) = bdiff - bchemo + bsource-bdeath;

            
        end % j=1:nyplus
    end % i=1:nxplus
    
end % ib=1:nb

correction=correction(:);
dy=[];
for ic=1:nc
    dy=[dy;dcc{ic}(:)];
end
for ib=1:nb
    dy=[dy;dbb{ib}(:)];
end


