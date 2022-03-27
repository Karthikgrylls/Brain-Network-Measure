


NORMALISED GRAPH THEORY -  PLOT RESULTS 


%% clear 
clc; clear; close all

%% Group/Condition one

path = 'C:\Research\CTL\'
cd(path);
SUBJlist_Group1 = dir('Subj*');

%% Absolute and Random and Normalised CC, PL, SW for Group1_ET data extraction
for i = 1:length(SUBJlist_Group1)
    %%
    SUBJname = SUBJlist_Group1(i).name;
    path1=([path SUBJname]);
    cd(path1);
    Group1fix_sub_name=SUBJname(1:end);
    data=load([Group1fix_sub_name '_ABS.mat'])
    Group1_CC(i,:,:)= data.GT_clust_coeff;         
    Group1_PL(i,:,:)= data.GT_path_length;
    Group1_LE(i,:,:)= data.GT_local_eff;         
    Group1_GE(i,:,:)= data.GT_global_eff;
    Group1_Degree(i,:,:)= data.GT_degree;
    Group1_PC(i,:,:)= data.GT_participation_coeff;  
    Group1_corr(i,:,:)=data.GT_corr_data;
     
    data_rand=load([Group1fix_sub_name '_RAND.mat'])
    Group1_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group1_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group1_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group1_GE_rand(i,:,:,:)= data_rand.GT_global_eff ;  %GT_global_eff_rand;
    Group1_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group1_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
    
    Group1_BC(i,:,:)= data.GT_betweenness; 
    Group1_ST(i,:,:)= data.GT_strengths;     
    Group1_Ass(i,:,:)= data.GT_assortativity;
    Group1_Dens(i,:,:)= data.GT_density;
    Group1_Modu(i,:,:)= data.GT_modularity;
    Group1_BC_rand(i,:,:,:)= data_rand.GT_betweenness_rand;     
    Group1_ST_rand(i,:,:,:)= data_rand.GT_strengths_rand; 
    Group1_Ass_rand(i,:,:,:)= data_rand.GT_assortativity_rand;     
    Group1_Dens_rand(i,:,:,:)= data_rand.GT_density_rand; 
    Group1_Modu_rand(i,:,:,:)= data_rand.GT_modularity_rand; 
end

%%
Spartcity_rng=length(data_rand.sparsity_val);
Group1_CC_50=Group1_CC(:,1:Spartcity_rng,:);                          
Group1_CC_rand_squ = squeeze(mean(Group1_CC_rand,3));        
AvgGroup1_CC=mean(mean(Group1_CC_50,3));

sparsity_CC_Group1_50 = (mean(Group1_CC_50,3));
sparsity_CC_rand_Group1 = mean(Group1_CC_rand_squ,3);
sparsity_CC_rand_Group1_50 = sparsity_CC_rand_Group1 (:,1:Spartcity_rng);

Group1_PL_50=Group1_PL(:,:,1:Spartcity_rng);
AvgGroup1_PL=squeeze(mean(Group1_PL_50));

Sparsity_PL_Group1_50=squeeze(Group1_PL_50);
Sparsity_PL_Group1_rand_squ = squeeze(mean(Group1_PL_rand,3)); 
Sparsity_PL_Group1_rand_50 = Sparsity_PL_Group1_rand_squ (:,1:Spartcity_rng);
Group1_Degree_50=Group1_Degree(:,1:Spartcity_rng,:);                             
AvgGroup1_Degree=mean(mean(Group1_Degree_50,3));
sparsity_Degree_Group1_50 = (mean(Group1_Degree_50,3));

Group1_PC_50=Group1_PC(:,1:Spartcity_rng,:);                          
Group1_PC_rand_squ = squeeze(mean(Group1_PC_rand,3));        
AvgGroup1_PC=mean(mean(Group1_PC_50,3));

sparsity_PC_Group1_50 = (mean(Group1_PC_50,3));
sparsity_PC_rand_Group1 = mean(Group1_PC_rand_squ,3);
sparsity_PC_rand_Group1_50 = sparsity_PC_rand_Group1 (:,1:Spartcity_rng);
%%%
Group1_BC_50=Group1_BC(:,1:Spartcity_rng,:);                          
Group1_BC_rand_squ = squeeze(mean(Group1_BC_rand,3));        
AvgGroup1_BC=mean(mean(Group1_BC_50,3));

sparsity_BC_Group1_50 = (mean(Group1_BC_50,3));
sparsity_BC_rand_Group1 = mean(Group1_BC_rand_squ,3);
sparsity_BC_rand_Group1_50 = sparsity_BC_rand_Group1 (:,1:Spartcity_rng);

Group1_ST_50=Group1_ST(:,1:Spartcity_rng,:);                          
Group1_ST_rand_squ = squeeze(mean(Group1_ST_rand,3));        
AvgGroup1_ST=mean(mean(Group1_ST_50,3));

sparsity_ST_Group1_50 = (mean(Group1_ST_50,3));
sparsity_ST_rand_Group1 = mean(Group1_ST_rand_squ,3);
sparsity_ST_rand_Group1_50 = sparsity_ST_rand_Group1 (:,1:Spartcity_rng);

Group1_Ass_50=Group1_Ass(:,1:Spartcity_rng);
AvgGroup1_Ass=squeeze(mean(Group1_Ass_50));
Sparsity_Ass_Group1_50=squeeze(Group1_Ass_50);
Sparsity_Ass_Group1_rand_squ = squeeze(mean(Group1_Ass_rand,3)); 
Sparsity_Ass_Group1_rand_50 = Sparsity_Ass_Group1_rand_squ (:,1:Spartcity_rng);

Group1_Dens_50=Group1_Dens(:,1:Spartcity_rng);
AvgGroup1_Dens=squeeze(mean(Group1_Dens_50));
Sparsity_Dens_Group1_50=squeeze(Group1_Dens_50);
Sparsity_Dens_Group1_rand_squ = squeeze(mean(Group1_Dens_rand,3)); 
Sparsity_Dens_Group1_rand_50 = Sparsity_Dens_Group1_rand_squ (:,1:Spartcity_rng);

Group1_Modu_50=Group1_Modu(:,:,1:Spartcity_rng);
AvgGroup1_Modu=squeeze(mean(Group1_Modu_50));
Sparsity_Modu_Group1_50=squeeze(Group1_Modu_50);
Sparsity_Modu_Group1_rand_squ = squeeze(mean(Group1_Modu_rand,3)); 
Sparsity_Modu_Group1_rand_50 = Sparsity_Modu_Group1_rand_squ (:,1:Spartcity_rng);
%%
for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_CC_normalised_Group1(i,j) = sparsity_CC_Group1_50(i,j)/sparsity_CC_rand_Group1_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_PL_normalised_Group1(i,j) = Sparsity_PL_Group1_50(i,j)/Sparsity_PL_Group1_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            SmallWorldNess_Group1(i,j) = sparsity_CC_normalised_Group1(i,j)/sparsity_PL_normalised_Group1(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_PC_normalised_Group1(i,j) = sparsity_PC_Group1_50(i,j)/sparsity_PC_rand_Group1_50(i,j); 
    end
end
%%%
for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_BC_normalised_Group1(i,j) = sparsity_BC_Group1_50(i,j)/sparsity_BC_rand_Group1_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_ST_normalised_Group1(i,j) = sparsity_ST_Group1_50(i,j)/sparsity_ST_rand_Group1_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_Ass_normalised_Group1(i,j) = Sparsity_Ass_Group1_50(i,j)/Sparsity_Ass_Group1_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_Dens_normalised_Group1(i,j) = Sparsity_Dens_Group1_50(i,j)/Sparsity_Dens_Group1_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_Modu_normalised_Group1(i,j) = Sparsity_Modu_Group1_50(i,j)/Sparsity_Modu_Group1_rand_50(i,j); 
    end
end
%%
Group1_LE_50=Group1_LE(:,1:Spartcity_rng,:);                          
Group1_LE_rand_squ = squeeze(mean(Group1_LE_rand,3));        
AvgGroup1_LE=mean(mean(Group1_LE_50,3));

sparsity_LE_Group1_50 = (mean(Group1_LE_50,3));
sparsity_LE_rand_Group1 = mean(Group1_LE_rand_squ,3);
sparsity_LE_rand_Group1_50 = sparsity_LE_rand_Group1 (:,1:Spartcity_rng);

Group1_GE_50=Group1_GE(:,1:Spartcity_rng);
AvgGroup1_GE=squeeze(mean(Group1_GE_50));

Sparsity_GE_Group1_50=squeeze(Group1_GE_50);
Sparsity_GE_Group1_rand_squ = squeeze(mean(Group1_GE_rand,3)); 
Sparsity_GE_Group1_rand_50 = Sparsity_GE_Group1_rand_squ (:,1:Spartcity_rng); 

%%
for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_LE_normalised_Group1(i,j) = sparsity_LE_Group1_50(i,j)/sparsity_LE_rand_Group1_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group1); 
    for j = 1:Spartcity_rng; 
            sparsity_GE_normalised_Group1(i,j) = Sparsity_GE_Group1_50(i,j)/Sparsity_GE_Group1_rand_50(i,j); 
    end
end




%% Group/Condition two study

path = 'C:\Research\ODN\'
cd(path);
SUBJlist_Group2 = dir('Subj*');
%%
for i = 1:length(SUBJlist_Group2)
    SUBJname = SUBJlist_Group2(i).name;
    path1=([path SUBJname]);
    cd(path1);
    Group1fix_sub_name=SUBJname(1:end);
    data=load([Group1fix_sub_name '_ABS.mat'])
    Group2_CC(i,:,:)= data.GT_clust_coeff;         
    Group2_PL(i,:,:)= data.GT_path_length;
    Group2_LE(i,:,:)= data.GT_local_eff;         
    Group2_GE(i,:,:)= data.GT_global_eff;
    Group2_Degree(i,:,:)= data.GT_degree;
    Group2_PC(i,:,:)= data.GT_participation_coeff;
    Group2_corr(i,:,:)=data.GT_corr_data;
    
    data_rand=load([Group1fix_sub_name '_RAND.mat'])
    Group2_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group2_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group2_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group2_GE_rand(i,:,:,:)= data_rand.GT_global_eff;
    Group2_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group2_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
    
    Group2_BC(i,:,:)= data.GT_betweenness; 
    Group2_ST(i,:,:)= data.GT_strengths;     
    Group2_Ass(i,:,:)= data.GT_assortativity;
    Group2_Dens(i,:,:)= data.GT_density;
    Group2_Modu(i,:,:)= data.GT_modularity;
    Group2_BC_rand(i,:,:,:)= data_rand.GT_betweenness_rand;     
    Group2_ST_rand(i,:,:,:)= data_rand.GT_strengths_rand; 
    Group2_Ass_rand(i,:,:,:)= data_rand.GT_assortativity_rand;     
    Group2_Dens_rand(i,:,:,:)= data_rand.GT_density_rand; 
    Group2_Modu_rand(i,:,:,:)= data_rand.GT_modularity_rand; 
end
%%
Group2_CC_50=Group2_CC(:,1:Spartcity_rng,:);                          
Group2_CC_rand_squ = squeeze(mean(Group2_CC_rand,3));        
AvgGroup2_CC=mean(mean(Group2_CC_50,3));

sparsity_CC_Group2_50 = (mean(Group2_CC_50,3));
sparsity_CC_rand_Group2 = mean(Group2_CC_rand_squ,3);
sparsity_CC_rand_Group2_50 = sparsity_CC_rand_Group2 (:,1:Spartcity_rng);

Group2_PL_50=Group2_PL(:,:,1:Spartcity_rng);
AvgGroup2_PL=squeeze(mean(Group2_PL_50));

Sparsity_PL_Group2_50=squeeze(Group2_PL_50);
Sparsity_PL_Group2_rand_squ = squeeze(mean(Group2_PL_rand,3)); 
Sparsity_PL_Group2_rand_50 = Sparsity_PL_Group2_rand_squ (:,1:Spartcity_rng);
Group2_Degree_50=Group2_Degree(:,1:Spartcity_rng,:);                             
AvgGroup2_Degree=mean(mean(Group2_Degree_50,3));
sparsity_Degree_Group2_50 = (mean(Group2_Degree_50,3));

Group2_PC_50=Group2_PC(:,1:Spartcity_rng,:);                          
Group2_PC_rand_squ = squeeze(mean(Group2_PC_rand,3));        
AvgGroup2_PC=mean(mean(Group2_PC_50,3));

sparsity_PC_Group2_50 = (mean(Group2_PC_50,3));
sparsity_PC_rand_Group2 = mean(Group2_PC_rand_squ,3);
sparsity_PC_rand_Group2_50 = sparsity_PC_rand_Group2 (:,1:Spartcity_rng);
%%%
Group2_BC_50=Group2_BC(:,1:Spartcity_rng,:);                          
Group2_BC_rand_squ = squeeze(mean(Group2_BC_rand,3));        
AvgGroup2_BC=mean(mean(Group2_BC_50,3));

sparsity_BC_Group2_50 = (mean(Group2_BC_50,3));
sparsity_BC_rand_Group2 = mean(Group2_BC_rand_squ,3);
sparsity_BC_rand_Group2_50 = sparsity_BC_rand_Group2 (:,1:Spartcity_rng);

Group2_ST_50=Group2_ST(:,1:Spartcity_rng,:);                          
Group2_ST_rand_squ = squeeze(mean(Group2_ST_rand,3));        
AvgGroup2_ST=mean(mean(Group2_ST_50,3));

sparsity_ST_Group2_50 = (mean(Group2_ST_50,3));
sparsity_ST_rand_Group2 = mean(Group2_ST_rand_squ,3);
sparsity_ST_rand_Group2_50 = sparsity_ST_rand_Group2 (:,1:Spartcity_rng);

Group2_Ass_50=Group2_Ass(:,1:Spartcity_rng);
AvgGroup2_Ass=squeeze(mean(Group2_Ass_50));
Sparsity_Ass_Group2_50=squeeze(Group2_Ass_50);
Sparsity_Ass_Group2_rand_squ = squeeze(mean(Group2_Ass_rand,3)); 
Sparsity_Ass_Group2_rand_50 = Sparsity_Ass_Group2_rand_squ (:,1:Spartcity_rng);

Group2_Dens_50=Group2_Dens(:,1:Spartcity_rng);
AvgGroup2_Dens=squeeze(mean(Group2_Dens_50));
Sparsity_Dens_Group2_50=squeeze(Group2_Dens_50);
Sparsity_Dens_Group2_rand_squ = squeeze(mean(Group2_Dens_rand,3)); 
Sparsity_Dens_Group2_rand_50 = Sparsity_Dens_Group2_rand_squ (:,1:Spartcity_rng);

Group2_Modu_50=Group2_Modu(:,:,1:Spartcity_rng);
AvgGroup2_Modu=squeeze(mean(Group2_Modu_50));
Sparsity_Modu_Group2_50=squeeze(Group2_Modu_50);
Sparsity_Modu_Group2_rand_squ = squeeze(mean(Group2_Modu_rand,3)); 
Sparsity_Modu_Group2_rand_50 = Sparsity_Modu_Group2_rand_squ (:,1:Spartcity_rng);
%%
for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_CC_normalised_Group2(i,j) = sparsity_CC_Group2_50(i,j)/sparsity_CC_rand_Group2_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_PL_normalised_Group2(i,j) = Sparsity_PL_Group2_50(i,j)/Sparsity_PL_Group2_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            SmallWorldNess_Group2(i,j) = sparsity_CC_normalised_Group2(i,j)/sparsity_PL_normalised_Group2(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_PC_normalised_Group2(i,j) = sparsity_PC_Group2_50(i,j)/sparsity_PC_rand_Group2_50(i,j); 
    end
end
%%%
for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_BC_normalised_Group2(i,j) = sparsity_BC_Group2_50(i,j)/sparsity_BC_rand_Group2_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_ST_normalised_Group2(i,j) = sparsity_ST_Group2_50(i,j)/sparsity_ST_rand_Group2_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_Ass_normalised_Group2(i,j) = Sparsity_Ass_Group2_50(i,j)/Sparsity_Ass_Group2_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_Dens_normalised_Group2(i,j) = Sparsity_Dens_Group2_50(i,j)/Sparsity_Dens_Group2_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_Modu_normalised_Group2(i,j) = Sparsity_Modu_Group2_50(i,j)/Sparsity_Modu_Group2_rand_50(i,j); 
    end
end
%%
Group2_LE_50=Group2_LE(:,1:Spartcity_rng,:);                          
Group2_LE_rand_squ = squeeze(mean(Group2_LE_rand,3));        
AvgGroup2_LE=mean(mean(Group2_LE_50,3));

sparsity_LE_Group2_50 = (mean(Group2_LE_50,3));
sparsity_LE_rand_Group2 = mean(Group2_LE_rand_squ,3);
sparsity_LE_rand_Group2_50 = sparsity_LE_rand_Group2 (:,1:Spartcity_rng);

Group2_GE_50=Group2_GE(:,1:Spartcity_rng);
AvgGroup2_GE=squeeze(mean(Group2_GE_50));

Sparsity_GE_Group2_50=squeeze(Group2_GE_50);
Sparsity_GE_Group2_rand_squ = squeeze(mean(Group2_GE_rand,3)); 
Sparsity_GE_Group2_rand_50 = Sparsity_GE_Group2_rand_squ (:,1:Spartcity_rng); 

%%
for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_LE_normalised_Group2(i,j) = sparsity_LE_Group2_50(i,j)/sparsity_LE_rand_Group2_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group2); 
    for j = 1:Spartcity_rng; 
            sparsity_GE_normalised_Group2(i,j) = Sparsity_GE_Group2_50(i,j)/Sparsity_GE_Group2_rand_50(i,j); 
    end
end
%%  %%%%group 3 study %%%%%%

path = 'C:\Research\ODP\'
cd(path);
SUBJlist_Group3 = dir('Subj*');
%%
for i = 1:length(SUBJlist_Group3)
    SUBJname = SUBJlist_Group3(i).name;
    path1=([path SUBJname]);
    cd(path1);
    Group1fix_sub_name=SUBJname(1:end);
    data=load([Group1fix_sub_name '_ABS.mat'])
    Group3_CC(i,:,:)= data.GT_clust_coeff;         
    Group3_PL(i,:,:)= data.GT_path_length;
    Group3_LE(i,:,:)= data.GT_local_eff;         
    Group3_GE(i,:,:)= data.GT_global_eff;
    Group3_Degree(i,:,:)= data.GT_degree;
    Group3_PC(i,:,:)= data.GT_participation_coeff;
    Group3_corr(i,:,:)=data.GT_corr_data;
    
    data_rand=load([Group1fix_sub_name '_RAND.mat'])
    Group3_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group3_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group3_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group3_GE_rand(i,:,:,:)= data_rand.GT_global_eff;
    Group3_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group3_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
    
    Group3_BC(i,:,:)= data.GT_betweenness; 
    Group3_ST(i,:,:)= data.GT_strengths;     
    Group3_Ass(i,:,:)= data.GT_assortativity;
    Group3_Dens(i,:,:)= data.GT_density;
    Group3_Modu(i,:,:)= data.GT_modularity;
    Group3_BC_rand(i,:,:,:)= data_rand.GT_betweenness_rand;     
    Group3_ST_rand(i,:,:,:)= data_rand.GT_strengths_rand; 
    Group3_Ass_rand(i,:,:,:)= data_rand.GT_assortativity_rand;     
    Group3_Dens_rand(i,:,:,:)= data_rand.GT_density_rand; 
    Group3_Modu_rand(i,:,:,:)= data_rand.GT_modularity_rand; 
end
%%
Group3_CC_50=Group3_CC(:,1:Spartcity_rng,:);                          
Group3_CC_rand_squ = squeeze(mean(Group3_CC_rand,3));        
AvgGroup3_CC=mean(mean(Group3_CC_50,3));

sparsity_CC_Group3_50 = (mean(Group3_CC_50,3));
sparsity_CC_rand_Group3 = mean(Group3_CC_rand_squ,3);
sparsity_CC_rand_Group3_50 = sparsity_CC_rand_Group3 (:,1:Spartcity_rng);

Group3_PL_50=Group3_PL(:,:,1:Spartcity_rng);
AvgGroup3_PL=squeeze(mean(Group3_PL_50));

Sparsity_PL_Group3_50=squeeze(Group3_PL_50);
Sparsity_PL_Group3_rand_squ = squeeze(mean(Group3_PL_rand,3)); 
Sparsity_PL_Group3_rand_50 = Sparsity_PL_Group3_rand_squ (:,1:Spartcity_rng);
Group3_Degree_50=Group3_Degree(:,1:Spartcity_rng,:);                             
AvgGroup3_Degree=mean(mean(Group3_Degree_50,3));
sparsity_Degree_Group3_50 = (mean(Group3_Degree_50,3));

Group3_PC_50=Group3_PC(:,1:Spartcity_rng,:);                          
Group3_PC_rand_squ = squeeze(mean(Group3_PC_rand,3));        
AvgGroup3_PC=mean(mean(Group3_PC_50,3));

sparsity_PC_Group3_50 = (mean(Group3_PC_50,3));
sparsity_PC_rand_Group3 = mean(Group3_PC_rand_squ,3);
sparsity_PC_rand_Group3_50 = sparsity_PC_rand_Group3 (:,1:Spartcity_rng);
%%%
Group3_BC_50=Group3_BC(:,1:Spartcity_rng,:);                          
Group3_BC_rand_squ = squeeze(mean(Group3_BC_rand,3));        
AvgGroup3_BC=mean(mean(Group3_BC_50,3));

sparsity_BC_Group3_50 = (mean(Group3_BC_50,3));
sparsity_BC_rand_Group3 = mean(Group3_BC_rand_squ,3);
sparsity_BC_rand_Group3_50 = sparsity_BC_rand_Group3 (:,1:Spartcity_rng);

Group3_ST_50=Group3_ST(:,1:Spartcity_rng,:);                          
Group3_ST_rand_squ = squeeze(mean(Group3_ST_rand,3));        
AvgGroup3_ST=mean(mean(Group3_ST_50,3));

sparsity_ST_Group3_50 = (mean(Group3_ST_50,3));
sparsity_ST_rand_Group3 = mean(Group3_ST_rand_squ,3);
sparsity_ST_rand_Group3_50 = sparsity_ST_rand_Group3 (:,1:Spartcity_rng);

Group3_Ass_50=Group3_Ass(:,1:Spartcity_rng);
AvgGroup3_Ass=squeeze(mean(Group3_Ass_50));
Sparsity_Ass_Group3_50=squeeze(Group3_Ass_50);
Sparsity_Ass_Group3_rand_squ = squeeze(mean(Group3_Ass_rand,3)); 
Sparsity_Ass_Group3_rand_50 = Sparsity_Ass_Group3_rand_squ (:,1:Spartcity_rng);

Group3_Dens_50=Group3_Dens(:,1:Spartcity_rng);
AvgGroup3_Dens=squeeze(mean(Group3_Dens_50));
Sparsity_Dens_Group3_50=squeeze(Group3_Dens_50);
Sparsity_Dens_Group3_rand_squ = squeeze(mean(Group3_Dens_rand,3)); 
Sparsity_Dens_Group3_rand_50 = Sparsity_Dens_Group3_rand_squ (:,1:Spartcity_rng);

Group3_Modu_50=Group3_Modu(:,:,1:Spartcity_rng);
AvgGroup3_Modu=squeeze(mean(Group3_Modu_50));
Sparsity_Modu_Group3_50=squeeze(Group3_Modu_50);
Sparsity_Modu_Group3_rand_squ = squeeze(mean(Group3_Modu_rand,3)); 
Sparsity_Modu_Group3_rand_50 = Sparsity_Modu_Group3_rand_squ (:,1:Spartcity_rng);
%%
for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_CC_normalised_Group3(i,j) = sparsity_CC_Group3_50(i,j)/sparsity_CC_rand_Group3_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_PL_normalised_Group3(i,j) = Sparsity_PL_Group3_50(i,j)/Sparsity_PL_Group3_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            SmallWorldNess_Group3(i,j) = sparsity_CC_normalised_Group3(i,j)/sparsity_PL_normalised_Group3(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_PC_normalised_Group3(i,j) = sparsity_PC_Group3_50(i,j)/sparsity_PC_rand_Group3_50(i,j); 
    end
end
%%%
for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_BC_normalised_Group3(i,j) = sparsity_BC_Group3_50(i,j)/sparsity_BC_rand_Group3_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_ST_normalised_Group3(i,j) = sparsity_ST_Group3_50(i,j)/sparsity_ST_rand_Group3_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_Ass_normalised_Group3(i,j) = Sparsity_Ass_Group3_50(i,j)/Sparsity_Ass_Group3_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_Dens_normalised_Group3(i,j) = Sparsity_Dens_Group3_50(i,j)/Sparsity_Dens_Group3_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_Modu_normalised_Group3(i,j) = Sparsity_Modu_Group3_50(i,j)/Sparsity_Modu_Group3_rand_50(i,j); 
    end
end
%%
Group3_LE_50=Group3_LE(:,1:Spartcity_rng,:);                          
Group3_LE_rand_squ = squeeze(mean(Group3_LE_rand,3));        
AvgGroup3_LE=mean(mean(Group3_LE_50,3));

sparsity_LE_Group3_50 = (mean(Group3_LE_50,3));
sparsity_LE_rand_Group3 = mean(Group3_LE_rand_squ,3);
sparsity_LE_rand_Group3_50 = sparsity_LE_rand_Group3 (:,1:Spartcity_rng);

Group3_GE_50=Group3_GE(:,1:Spartcity_rng);
AvgGroup3_GE=squeeze(mean(Group3_GE_50));

Sparsity_GE_Group3_50=squeeze(Group3_GE_50);
Sparsity_GE_Group3_rand_squ = squeeze(mean(Group3_GE_rand,3)); 
Sparsity_GE_Group3_rand_50 = Sparsity_GE_Group3_rand_squ (:,1:Spartcity_rng); 

%%
for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_LE_normalised_Group3(i,j) = sparsity_LE_Group3_50(i,j)/sparsity_LE_rand_Group3_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group3); 
    for j = 1:Spartcity_rng; 
            sparsity_GE_normalised_Group3(i,j) = Sparsity_GE_Group3_50(i,j)/Sparsity_GE_Group3_rand_50(i,j); 
    end
end

%% %Ploting Normalized CC, PC, PL, SW, GE, LE and other graph measures Images
CI = 1.96; % Z values of 95% Confidence Interval
sparsity_val = [0.1000    0.1500    0.2000    0.2500    0.3000    0.3500    0.4000    0.4500    0.5000]
close all
figure;
y1_CC = mean (sparsity_CC_normalised_Group1);
z1_CC = std (sparsity_CC_normalised_Group1)/sqrt (length (sparsity_CC_normalised_Group1)); 
errorbar (y1_CC,CI*z1_CC, 'k','LineWidth',3); grid on; 
hold on
y2_CC = mean (sparsity_CC_normalised_Group2);
z2_CC = std (sparsity_CC_normalised_Group2)/sqrt (length (sparsity_CC_normalised_Group2));
errorbar (y2_CC,CI*z2_CC, 'b','LineWidth',3); grid on; 
y3_CC = mean (sparsity_CC_normalised_Group3);
z3_CC = std (sparsity_CC_normalised_Group3)/sqrt (length (sparsity_CC_normalised_Group3));
errorbar (y3_CC,CI*z3_CC, 'r','LineWidth',3); grid on; 
title('Whole Brain Segregation (Clustering Coefficient)');
xlabel('Sparsity', 'FontSize', 24);
ylabel('Average Clustering Coefficient', 'FontSize', 24);
xticklabels({sparsity_val})
legend('HC', 'ODN', 'ODP', 'FontWeight', 'bold')
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
%
figure;
y1_PC = mean (sparsity_PC_normalised_Group1);
z1_PC = std (sparsity_PC_normalised_Group1)/sqrt (length (sparsity_PC_normalised_Group1)); 
errorbar (y1_PC,CI*z1_PC, 'k','LineWidth',3); grid on; 
hold on
y2_PC = mean (sparsity_PC_normalised_Group2);
z2_PC = std (sparsity_PC_normalised_Group2)/sqrt (length (sparsity_PC_normalised_Group2));
errorbar (y2_PC,CI*z2_PC, 'b','LineWidth',3); grid on; 
y3_PC = mean (sparsity_PC_normalised_Group3);
z3_PC = std (sparsity_PC_normalised_Group3)/sqrt (length (sparsity_PC_normalised_Group3));
errorbar (y3_PC,CI*z3_PC, 'r','LineWidth',3); grid on; 
title('Whole Brain Integration (Participation Coefficient)');
xlabel('Sparsity', 'FontSize', 24);
ylabel('Average Participation Coefficient', 'FontSize', 24);
xticklabels({sparsity_val})
legend('HC', 'ODN', 'ODP', 'FontWeight', 'bold')
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
%
figure;
hold on
y1_PL = mean (sparsity_PL_normalised_Group1);
z1_PL = std (sparsity_PL_normalised_Group1)/sqrt (length (sparsity_PL_normalised_Group1)); 
errorbar (y1_PL,CI*z1_PL, 'k','LineWidth',3); grid on; 
hold on
y2_PL = mean (sparsity_PL_normalised_Group2);
z2_PL = std (sparsity_PL_normalised_Group2)/sqrt (length (sparsity_PL_normalised_Group2));
errorbar (y2_PL,CI*z2_PL, 'b','LineWidth',3); grid on; 
y3_PL = mean (sparsity_PL_normalised_Group3);
z3_PL = std (sparsity_PL_normalised_Group3)/sqrt (length (sparsity_PL_normalised_Group3));
errorbar (y3_PL,CI*z3_PL, 'r','LineWidth',3); grid on; 
title('Whole Brain Average Path Length');
xlabel('Sparsity', 'FontSize', 24);
ylabel('Average Path Length', 'FontSize', 24);
xticklabels({sparsity_val})
legend('HC', 'ODN', 'ODP', 'FontWeight', 'bold')
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')

figure;
y1_SW = mean (SmallWorldNess_Group1);
z1_SW = std (SmallWorldNess_Group1)/sqrt (length (SmallWorldNess_Group1));
errorbar (y1_SW,CI*z1_SW, 'k','LineWidth',3); grid on; 
hold on
y2_SW = mean (SmallWorldNess_Group2);
z2_SW = std (SmallWorldNess_Group2)/sqrt (length (SmallWorldNess_Group2));
errorbar (y2_SW,CI*z2_SW, 'b','LineWidth',3); grid on; 
y3_SW = mean (SmallWorldNess_Group3);
z3_SW = std (SmallWorldNess_Group3)/sqrt (length (SmallWorldNess_Group3));
errorbar (y3_SW,CI*z3_SW, 'r','LineWidth',3); grid on; 
title('Brain Small Worldness');
xlabel('Sparsity', 'FontSize', 24);
ylabel('Small Worldness', 'FontSize', 24);
xticklabels({sparsity_val})
legend('HC', 'ODN', 'ODP', 'FontWeight', 'bold')
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')

figure;
y1_CC = mean (sparsity_LE_normalised_Group1);
z1_CC = std (sparsity_LE_normalised_Group1)/sqrt (length (sparsity_LE_normalised_Group1)); 
errorbar (y1_CC,CI*z1_CC, 'k','LineWidth',3); grid on; 
hold on
y2_CC = mean (sparsity_LE_normalised_Group2);
z2_CC = std (sparsity_LE_normalised_Group2)/sqrt (length (sparsity_LE_normalised_Group2));
errorbar (y2_CC,CI*z2_CC, 'b','LineWidth',3); grid on; 
y3_CC = mean (sparsity_LE_normalised_Group3);
z3_CC = std (sparsity_LE_normalised_Group3)/sqrt (length (sparsity_LE_normalised_Group3));
errorbar (y3_CC,CI*z3_CC, 'r','LineWidth',3); grid on; 
title('Whole Brain Local Eficency');
xlabel('Sparsity', 'FontSize', 24);
ylabel('Average Local Eficency', 'FontSize', 24);
xticklabels({sparsity_val})
legend('HC', 'ODN', 'ODP', 'FontWeight', 'bold')
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')

figure;
hold on
y1_PL = mean (sparsity_GE_normalised_Group1);
z1_PL = std (sparsity_GE_normalised_Group1)/sqrt (length (sparsity_GE_normalised_Group1)); 
errorbar (y1_PL,CI*z1_PL, 'k','LineWidth',3); grid on; 
hold on
y2_PL = mean (sparsity_GE_normalised_Group2);
z2_PL = std (sparsity_GE_normalised_Group2)/sqrt (length (sparsity_GE_normalised_Group2));
errorbar (y2_PL,CI*z2_PL, 'b','LineWidth',3); grid on; 
y3_PL = mean (sparsity_GE_normalised_Group3);
z3_PL = std (sparsity_GE_normalised_Group3)/sqrt (length (sparsity_GE_normalised_Group3));
errorbar (y3_PL,CI*z3_PL, 'r','LineWidth',3); grid on; 
title('Whole Brain Global Eficency');
xlabel('Sparsity', 'FontSize', 24);
ylabel('Average Global Eficency', 'FontSize', 24);
xticklabels({sparsity_val})
legend('HC', 'ODN', 'ODP', 'FontWeight', 'bold')
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
%
figure;
y1_Modu = mean (sparsity_Modu_normalised_Group1);
z1_Modu = std (sparsity_Modu_normalised_Group1)/sqrt (length (sparsity_Modu_normalised_Group1)); 
errorbar (y1_Modu,CI*z1_Modu, 'k','LineWidth',3); grid on; 
hold on
y2_Modu = mean (sparsity_Modu_normalised_Group2);
z2_Modu = std (sparsity_Modu_normalised_Group2)/sqrt (length (sparsity_Modu_normalised_Group2));
errorbar (y2_Modu,CI*z2_Modu, 'b','LineWidth',3); grid on; 
y3_Modu = mean (sparsity_Modu_normalised_Group3);
z3_Modu = std (sparsity_Modu_normalised_Group3)/sqrt (length (sparsity_Modu_normalised_Group3));
errorbar (y3_Modu,CI*z3_Modu, 'r','LineWidth',3); grid on; 
title('Whole Brain Modularity');
xlabel('Sparsity', 'FontSize', 24);
ylabel('Average Modularity', 'FontSize', 24);
xticklabels({sparsity_val})
legend('HC', 'ODN', 'ODP', 'FontWeight', 'bold')
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')

figure;
y1_BC = mean (sparsity_BC_Group1_50);
z1_BC = std (sparsity_BC_Group1_50)/sqrt (length (sparsity_BC_Group1_50)); 
errorbar (y1_BC,CI*z1_BC, 'k','LineWidth',3); grid on; 
hold on
y2_BC = mean (sparsity_BC_Group2_50);
z2_BC = std (sparsity_BC_Group2_50)/sqrt (length (sparsity_BC_Group2_50));
errorbar (y2_BC,CI*z2_BC, 'b','LineWidth',3); grid on;  
y3_BC = mean (sparsity_BC_Group3_50);
z3_BC = std (sparsity_BC_Group3_50)/sqrt (length (sparsity_BC_Group3_50));
errorbar (y3_BC,CI*z3_BC, 'r','LineWidth',3); grid on; 
title('Whole Brain Network Centrality');
xlabel('Sparsity', 'FontSize', 24);
ylabel('Average Betweenness Centrality', 'FontSize', 24);
xticklabels({sparsity_val})
legend('HC', 'ODN', 'ODP', 'FontWeight', 'bold')
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
%
figure;
y1_Ass = mean (Sparsity_Ass_Group1_50);
z1_Ass = std (Sparsity_Ass_Group1_50)/sqrt (length (Sparsity_Ass_Group1_50)); 
errorbar (y1_Ass,CI*z1_Ass, 'k','LineWidth',3); grid on;  
hold on
y2_Ass = mean (Sparsity_Ass_Group2_50);
z2_Ass = std (Sparsity_Ass_Group2_50)/sqrt (length (Sparsity_Ass_Group2_50));
errorbar (y2_Ass,CI*z2_Ass, 'b','LineWidth',3); grid on; 
y3_Ass = mean (Sparsity_Ass_Group3_50);
z3_Ass = std (Sparsity_Ass_Group3_50)/sqrt (length (Sparsity_Ass_Group3_50));
errorbar (y3_Ass,CI*z3_Ass, 'r','LineWidth',3); grid on; 
title('Brain Assortativity');
xlabel('Sparsity', 'FontSize', 24);
ylabel('Average Assortativity', 'FontSize', 24);
xticklabels({sparsity_val})
legend('HC', 'ODN', 'ODP', 'FontWeight', 'bold')
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
%%
close all
SR=1:9;
Grp_con = ones(length(SUBJlist_Group1),1); Grp_mcs = 2*(ones(length(SUBJlist_Group2),1)); Grp_uws = 3*(ones(length(SUBJlist_Group3),1));
Grp = [Grp_con; Grp_mcs; Grp_uws]; 
Grp_CC = [mean(sparsity_CC_normalised_Group1(:,SR),2); mean(sparsity_CC_normalised_Group2(:,SR),2); mean(sparsity_CC_normalised_Group3(:,SR),2)];
%%
figure();
subplot(1,3,1); notBoxPlot(Grp_CC,Grp,0.5,'patch',ones(length(Grp_CC),1));
title('Whole Brain Segregation'); ylabel('Average Clustering Coefficient'); xticklabels({'HC','ODN','ODP'})
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
Grp_PC = [mean(sparsity_PC_normalised_Group1(:,:),2); mean(sparsity_PC_normalised_Group2(:,:),2); mean(sparsity_PC_normalised_Group3(:,:),2)];
subplot(1,3,2); notBoxPlot(Grp_PC,Grp,0.5,'patch',ones(length(Grp_PC),1));
title('Whole Brain Integration'); ylabel('Average Participation Coefficient'); xticklabels({'HC','ODN','ODP'})
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
Grp_SW = [mean(SmallWorldNess_Group1(:,SR),2); mean(SmallWorldNess_Group2(:,SR),2); mean(SmallWorldNess_Group3(:,SR),2)];
subplot(1,3,3); notBoxPlot(Grp_SW,Grp,0.5,'patch',ones(length(Grp_CC),1));
title('Brain Small-Worldness'); ylabel('Small-Worldness'); xticklabels({'HC','ODN','ODP'})
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
%% corelation matrix
Group1_corr_m=squeeze(mean(Group1_corr,1));
Group2_corr_m=squeeze(mean(Group2_corr,1));
Group3_corr_m=squeeze(mean(Group3_corr,1));
figure; imagesc(Group1_corr_m)
title('HC'); xlabel('ROIs'); ylabel('ROIs'); set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
figure; imagesc(Group2_corr_m)
title('ODN'); xlabel('ROIs'); ylabel('ROIs'); set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
figure; imagesc(Group3_corr_m)
title('ODP'); xlabel('ROIs'); ylabel('ROIs'); set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
%% --------------t-stats between Group2 and Group1 sparsity level------------%

for i = 1:Spartcity_rng
    [h_CC1(i),p_CC1(i)] = ttest2(sparsity_CC_normalised_Group2(:,i),sparsity_CC_normalised_Group1(:,i),0.05,'right');
end

h_CC1
p_CC1
%%
for i = 1:Spartcity_rng
    [h_GE1(i),p_GE1(i)] = ttest2(sparsity_GE_normalised_Group3(:,i),sparsity_GE_normalised_Group2(:,i),0.05,'left');
end
h_GE1
p_GE1
%%
for i = 1:Spartcity_rng
    [h_SW1(i),p_SW1(i)] = ttest2(SmallWorldNess_Group2(:,i),SmallWorldNess_Group1(:,i),0.05,'right');
end

h_SW1
p_SW1
%%
for i = 1:Spartcity_rng
    [h_PC1(i),p_PC1(i)] = ttest2(sparsity_PC_normalised_Group3(:,i),sparsity_PC_normalised_Group2(:,i),0.05,'left');
end

h_PC1
p_PC1
%%
%----------------------%% Brain resion significant Computations for CC %-------------%


sparsity_CC_Group1_ROI = squeeze(mean(Group1_CC_50,2)); 
sparsity_CC_Group2_ROI = squeeze(mean(Group2_CC_50,2));
sparsity_CC_Group3_ROI = squeeze(mean(Group3_CC_50,2));
sparsity_CC_rand_Group1_ROI = mean(Group1_CC_rand_squ,2);
sparsity_CC_rand_Group2_ROI = mean(Group2_CC_rand_squ,2);
sparsity_CC_rand_Group3_ROI = mean(Group3_CC_rand_squ,2);


for i = 1:length(SUBJlist_Group1); 
    for j = 1:160; 
            sparsity_CC_normalised_Group1_ROI(i,j) = sparsity_CC_Group1_ROI(i,j)/sparsity_CC_rand_Group1_ROI(i,j);
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:160; 
            sparsity_CC_normalised_Group2_ROI(i,j) = sparsity_CC_Group2_ROI(i,j)/sparsity_CC_rand_Group2_ROI(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:160; 
            sparsity_CC_normalised_Group3_ROI(i,j) = sparsity_CC_Group3_ROI(i,j)/sparsity_CC_rand_Group3_ROI(i,j); 
    end
end
%%
%----------------------%% Brain resion significant Computations for PC %-------------%

sparsity_PC_Group1_ROI = squeeze(mean(Group1_PC_50,2)); 
sparsity_PC_Group2_ROI = squeeze(mean(Group2_PC_50,2)); 
sparsity_PC_Group3_ROI = squeeze(mean(Group3_PC_50,2));
sparsity_PC_rand_Group1_ROI = mean(Group1_PC_rand_squ,2);
sparsity_PC_rand_Group2_ROI = mean(Group2_PC_rand_squ,2);
sparsity_PC_rand_Group3_ROI = mean(Group3_PC_rand_squ,2);



for i = 1:length(SUBJlist_Group1); 
    for j = 1:160; 
            sparsity_PC_normalised_Group1_ROI(i,j) = sparsity_PC_Group1_ROI(i,j)/sparsity_PC_rand_Group1_ROI(i,j);
    end
end



for i = 1:length(SUBJlist_Group2); 
    for j = 1:160; 
            sparsity_PC_normalised_Group2_ROI(i,j) = sparsity_PC_Group2_ROI(i,j)/sparsity_PC_rand_Group2_ROI(i,j); 
    end
end

for i = 1:length(SUBJlist_Group3); 
    for j = 1:160; 
            sparsity_PC_normalised_Group3_ROI(i,j) = sparsity_PC_Group3_ROI(i,j)/sparsity_PC_rand_Group3_ROI(i,j); 
    end
end
%% Brain resion significant Computations for CC
clc
for i = 1:160
    [h_CC_ROI_Normalised(i),p_CC_ROI_Normalised(i)] = ttest2(sparsity_CC_normalised_Group1_ROI(:,i),sparsity_CC_normalised_Group2_ROI(:,i),0.005,'right');
end

h_CC_ROI_Normalised
p_CC_ROI_Normalised
%%
clc
for i = 1:160
    [h_PC_ROI_Normalised(i),p_PC_ROI_Normalised(i)] = ttest2(sparsity_PC_normalised_Group2_ROI(:,i),sparsity_PC_normalised_Group3_ROI(:,i),0.001,'left');
end

h_PC_ROI_Normalised   
p_PC_ROI_Normalised
%%[p,h]=fdr(p_CC_ROI_Normalised,0.05);
%p

%%  Graphical Plot od GT measures %%%
% % % tpz_cc_control = trapz(sparsity_CC_normalised_Control(:,1:30),2)/30;
% % % tpz_cc_Patient = trapz(sparsity_CC_normalised_Patient(:,1:30),2)/30;
% % % tpz_pl_control = trapz(sparsity_PL_normalised_Control(:,1:30),2)/30;
% % % tpz_pl_Patient = trapz(sparsity_PL_normalised_Patient(:,1:30),2)/30;
% % % tpz_sw_control = trapz(SmallWorldNess_Control(:,1:30),2)/30;
% % % tpz_sw_Patient = trapz(SmallWorldNess_Patient(:,1:30),2)/30;
% % % tpz_le_control = trapz(sparsity_LE_normalised_Control(:,1:30),2)/30;
% % % tpz_le_Patient = trapz(sparsity_LE_normalised_Patient(:,1:30),2)/30;
% % % tpz_ge_control = trapz(sparsity_GE_normalised_Control(:,1:30),2)/30;
% % % tpz_ge_Patient = trapz(sparsity_GE_normalised_Patient(:,1:30),2)/30;
% % % 
% % % tpz_GT = [tpz_cc_control tpz_cc_Patient tpz_pl_control tpz_pl_Patient tpz_sw_control tpz_sw_Patient tpz_le_control tpz_le_Patient tpz_ge_control tpz_ge_Patient]
% % % 
% % % [p,h]=ttest2(tpz_GT(:,1),tpz_GT(:,2), 0.05,'left')
% % % 
% % % 
