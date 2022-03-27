
GRAPH THEORY MEASURE -RANDOM 

clc
clear
path='C:\Research\ODP\';
cd(path)
SUBJlist=dir('Subj*');
%%
for i=1:length(SUBJlist)
    tic
    SUBJname=SUBJlist(i).name;
    path1=([path SUBJname])
    cd(path1);
    %% 
    filelist= ([ SUBJname ]);  
    for i1=1%%%%:length(filelist)
        filename=filelist; 
 %       prefix_name=filename(1:end);
                
        final_data=load(['dataanalysis_ODP' filename(5:end)]);
        final_data=final_data.y_roi_regressed; %y_roi_regressed_filtered;

        GT_corr_data=corr(final_data);
        GT_corr_data_abs = abs(GT_corr_data);
        chanlocs = size(final_data,2); % No  of ROIs                      
        %% Thresholding
         sparsity_val=0.1:0.05:0.52; %0.1:0.025:0.52; %sparsity_val=0.01:0.05:0.52; %sparsity_val=0.01:0.025:1;
        for i2=1:length (sparsity_val)                   %i2=1:40
            %% %%%Network properties/ network measurement 
            GT_sparsity(i2)=sparsity_val(i2); % define sparsity
            corr_data_thr1=threshold_proportional(GT_corr_data_abs,GT_sparsity(i2)); % calcualte binary matrix
            corr_data_thr_bin=weight_conversion(corr_data_thr1,'binarize'); % binary weight matrix
            corr_data_thr = weight_conversion(corr_data_thr_bin, 'autofix'); % removing NaN & Inf
                       
            for random_number=1:40   % random_number=1:50
                             
                random_network=randmio_und(corr_data_thr,5);
                %%
                GT_corr_data_rand_thr(i2,random_number,:,:)=random_network;% asign to different array
                
                GT_degree_rand(i2,random_number,:)=degrees_und(random_network); % calculate degree
                
                GT_clust_coeff_rand(i2,random_number,:)=clustering_coef_bu(random_network);%% calculate clustering coeff
                                
                GT_local_eff_rand(i2,random_number,:)=efficiency_bin(random_network,1); % local efficiency
                
                GT_global_eff(i2,random_number,:)=efficiency_bin(random_network,0); % global efficiency
                
                GT_distance_matrix_rand(i2,random_number,:,:)=distance_bin(random_network); % distance matrix
                
                GT_path_length_rand(i2,random_number) =charpath(squeeze(GT_distance_matrix_rand(i2,random_number,:,:)),1,0); % path length

                GT_betweenness_rand(i2,random_number,:) = betweenness_bin(corr_data_thr); %betweenness centrality
                GT_strengths_rand(i2,random_number,:) = strengths_und(corr_data_thr);
                [~,~,a]=density_und(corr_data_thr);
                GT_density_rand(i2,random_number,:) = a; clear a
                GT_assortativity_rand(i2,random_number,:) = assortativity_bin(corr_data_thr,0);
            %%%%% Participation coefficient and modspan
                param.heuristic=50;
                for i = 1:param.heuristic
                    [Ci, allQ(i2,random_number,i)] = community_louvain(random_network);
                     allCi(i2,random_number,i,:) = Ci;
                     allpc(i2,random_number,i,:) = participation_coef(random_network,Ci); 
                end
                GT_modularity_rand(i2,random_number)= mean(allQ(i2,random_number,:));  % modularity
                GT_community_structure_rand(i2,random_number,1:chanlocs) = squeeze(allCi(i2,random_number,1,:)); % community structure
                GT_participation_coeff_rand(i2,random_number,1:chanlocs) = mean(squeeze(allpc(i2,random_number,:,:)));  %participation coefficient
                %%%%%    
            end
        end
        varname=([SUBJname '_RAND'])
        save(varname);
    end
    cd ..
end
