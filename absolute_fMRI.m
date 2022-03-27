


GRAPH THEORY MEASURE - ABSOLUTE 

clc
clear
path='C:\Research\ODN\';
cd(path)
SUBJlist=dir('Subj*');
%%
for i=1:length(SUBJlist)
    %tic
    SUBJname=SUBJlist(i).name;
    path1=([path SUBJname])
    cd(path1);
    %% 
    filelist= ([ SUBJname ]);  
    for i1=1%%%%:length(filelist)
        filename=filelist; 
 %       prefix_name=filename(1:end);
                
        final_data=load(['dataanalysis_ODN' filename(5:end)]);
        final_data=final_data.y_roi_regressed_filtered;

        GT_corr_data=corr(final_data);
        GT_corr_data_abs = abs(GT_corr_data);
        chanlocs = size(final_data,2); % No  of ROIs       
        %% Thresholding
         sparsity_val=0.1:0.05:0.52; %0.1:0.025:0.52; %sparsity_val=0.01:0.025:10.1:0.025:0.52; %sparsity_val=0.01:0.05:0.52; %sparsity_val=0.01:0.025:1;
        for i2=1:length (sparsity_val)                   %i2=1:40
            %% %%%Network properties/ network measurement 
            GT_sparsity(i2)=sparsity_val(i2); % define sparsity
            corr_data_thr1=threshold_proportional(GT_corr_data_abs,GT_sparsity(i2)); % calcualte binary matrix
            corr_data_thr_bin=weight_conversion(corr_data_thr1,'binarize'); % binary weight matrix
            corr_data_thr = weight_conversion(corr_data_thr_bin, 'autofix'); % removing NaN & Inf
            
            %%
            GT_corr_data_thr(i2,:,:)=corr_data_thr;% asign to different array
            
            GT_degree(i2,:)=degrees_und(corr_data_thr); % calculate degree
            
            GT_clust_coeff(i2,:)=clustering_coef_bu(corr_data_thr);%% calculate clustering coeff
            
            GT_local_eff(i2,:)=efficiency_bin(corr_data_thr,1); % local efficiency
            
            GT_global_eff(i2,:)=efficiency_bin(corr_data_thr,0); % global efficiency
            
            GT_distance_matrix(i2,:,:)=distance_bin(corr_data_thr); % distance matrix
            
            GT_path_length(i2)=charpath(squeeze(GT_distance_matrix(i2,:,:)),1,0); % path length

            GT_betweenness(i2,:) = betweenness_bin(corr_data_thr); %betweenness centrality
            GT_strengths(i2,:) = strengths_und(corr_data_thr);
            [~,~,a]=density_und(corr_data_thr);
            GT_density(i2,:) = a; clear a
            GT_assortativity(i2,:) = assortativity_bin(corr_data_thr,0);
            %GT_assortativity(i2,:,:) = search_information(corr_data_thr);
            %%%%% Participation coefficient and modspan
            param.heuristic=50;
            
            for i = 1:param.heuristic
                [Ci, allQ(i2,i)] = community_louvain(corr_data_thr);
            
                allCi(i2,i,:) = Ci;
   
                allpc(i2,i,:) = participation_coef(corr_data_thr,Ci); 
            end
        
            GT_modularity(i2)= mean(allQ(i2,:));  % modularity
            GT_community_structure(i2,1:chanlocs) = squeeze(allCi(i2,1,:)); % community structure
            GT_participation_coeff(i2,1:chanlocs) = mean(squeeze(allpc(i2,:,:)));  %participation coefficient
            
            %%%%%
            
        end
        
        varname=([SUBJname '_ABS'])
        save(varname);
    end
    cd ..
    %toc
end