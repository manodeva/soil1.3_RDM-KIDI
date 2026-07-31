% Computing relative abundance and biomass conc. of each species from microscopic image data using GlcNAc as C source 
clc
clear
close all

%% Loading the image data

% data_nag and data_chito5 is a struct variable with three fields DAPI,
% BFP, GFP containing the normalized fluorescence intensities of Rhodo,
% Sphingo, and Vario as matrices for different time points 
% time - vector of time points when image data was measured

load('imgdata_feb2023.mat','data_nag','time')

%% Info related to data

[nrow,ncol] = size(data_nag.DAPI{1}); % size of the data matrix
nt = length(time);
b_fields = fieldnames(data_nag);

%% find overall max intensity across all time points and microbes

overall_max_intensity = 0;

for j = 1:numel(b_fields)
    b_field = b_fields{j};

    for k = 1:nt
        overall_max_intensity = max(overall_max_intensity, max(data_nag.(b_field){k}, [], 'all'));
    end
end

%% Concentration of biomass computation (in g/L)
% Assumes same bacterial weight for all 3 microbes, and no. of layers 

ind_bac_wt = 3e-13;  % Individual bacteria weight in g/cell
nl = 3; % no. of layers
for k = 1:nt
    for j = 1:numel(b_fields)

        b_field = b_fields{j};

        intensity = data_nag.(b_field){k};

        % Convert fluorescence to cells/µm³
        cell_density = (intensity./overall_max_intensity) .* (1/nl);

        % Convert cells/µm³ to g/µm³
        biomass_density = cell_density .* ind_bac_wt;

        % Convert g/µm³ to g/L
        data_nag.cell_conc{k}(:,:,j) = biomass_density ./ 1e-15; 

    end
end

%% Relative abundance computation

for k = 1:nt

    intensity = zeros(nrow*ncol, numel(b_fields));

    for j = 1:numel(b_fields)
        b_field = b_fields{j};
        intensity(:,j) = data_nag.(b_field){k}(:);
    end

    total_intensity = sum(intensity, 2);
    rela_abd_nag{k} = zeros(size(intensity));
    nonzero = total_intensity > 0; % ignoring 0/0

    rela_abd_nag{k}(nonzero,:) = intensity(nonzero,:) ./ total_intensity(nonzero,:);

end