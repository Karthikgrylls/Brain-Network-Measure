


FEATURE RANKING OF GRAPH THEORY



cc_g1 = trapz(sparsity_CC_normalised_Group1(:,:),2)/size(sparsity_PC_normalised_Group1,2);
cc_g2 = trapz(sparsity_CC_normalised_Group3(:,:),2)/size(sparsity_PC_normalised_Group1,2);
cc_g3 = trapz(sparsity_CC_normalised_Group2(:,:),2)/size(sparsity_PC_normalised_Group1,2);

pc_g1 = trapz(sparsity_PC_normalised_Group1(:,:),2)/size(sparsity_PC_normalised_Group1,2);
pc_g2 = trapz(sparsity_PC_normalised_Group3(:,:),2)/size(sparsity_PC_normalised_Group1,2);
pc_g3 = trapz(sparsity_PC_normalised_Group2(:,:),2)/size(sparsity_PC_normalised_Group1,2);

sw_g1 = trapz(SmallWorldNess_Group1(:,:),2)/size(sparsity_PC_normalised_Group1,2);
sw_g2 = trapz(SmallWorldNess_Group2(:,:),2)/size(sparsity_PC_normalised_Group1,2);
sw_g3 = trapz(SmallWorldNess_Group3(:,:),2)/size(sparsity_PC_normalised_Group1,2);

pl_g1 = trapz(sparsity_PL_normalised_Group1(:,:),2)/size(sparsity_PC_normalised_Group1,2);
pl_g2 = trapz(sparsity_PL_normalised_Group3(:,:),2)/size(sparsity_PC_normalised_Group1,2);
pl_g3 = trapz(sparsity_PL_normalised_Group2(:,:),2)/size(sparsity_PC_normalised_Group1,2);

le_g1 = trapz(sparsity_LE_normalised_Group1(:,:),2)/size(sparsity_PC_normalised_Group1,2);
le_g2 = trapz(sparsity_LE_normalised_Group3(:,:),2)/size(sparsity_PC_normalised_Group1,2);
le_g3 = trapz(sparsity_LE_normalised_Group2(:,:),2)/size(sparsity_PC_normalised_Group1,2);

ge_g1 = trapz(sparsity_GE_normalised_Group1(:,:),2)/size(sparsity_PC_normalised_Group1,2);
ge_g2 = trapz(sparsity_GE_normalised_Group3(:,:),2)/size(sparsity_PC_normalised_Group1,2);
ge_g3 = trapz(sparsity_GE_normalised_Group2(:,:),2)/size(sparsity_PC_normalised_Group1,2);

md_g1 = trapz(sparsity_Modu_normalised_Group1(:,:),2)/size(sparsity_PC_normalised_Group1,2);
md_g2 = trapz(sparsity_Modu_normalised_Group3(:,:),2)/size(sparsity_PC_normalised_Group1,2);
md_g3 = trapz(sparsity_Modu_normalised_Group2(:,:),2)/size(sparsity_PC_normalised_Group1,2);

bc_g1 = trapz(sparsity_BC_Group1_50(:,:),2)/size(sparsity_PC_normalised_Group1,2);
bc_g2 = trapz(sparsity_BC_Group3_50(:,:),2)/size(sparsity_PC_normalised_Group1,2);
bc_g3 = trapz(sparsity_BC_Group2_50(:,:),2)/size(sparsity_PC_normalised_Group1,2);

as_g1 = trapz(Sparsity_Ass_Group1_50(:,:),2)/size(sparsity_PC_normalised_Group1,2);
as_g2 = trapz(Sparsity_Ass_Group3_50(:,:),2)/size(sparsity_PC_normalised_Group1,2);
as_g3 = trapz(Sparsity_Ass_Group2_50(:,:),2)/size(sparsity_PC_normalised_Group1,2);

corr_g1 = mean(mean(Group1_corr,3),2);
corr_g2 = mean(mean(Group2_corr,3),2);
corr_g3 = mean(mean(Group3_corr,3),2);

g1 = [cc_g1 pc_g1 sw_g1 pl_g1 le_g1 ge_g1 md_g1 bc_g1 as_g1 corr_g1];
g2 = [cc_g2 pc_g2 sw_g2 pl_g2 le_g2 ge_g2 md_g2 bc_g2 as_g2 corr_g1];
g3 = [cc_g3 pc_g3 sw_g3 pl_g3 le_g3 ge_g3 md_g3 bc_g3 as_g3 corr_g1];


% % gPar = [g2; g3];
% % gX = gPar;
% % species={'CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT' ...
% %     'PNT1', 'PNT1', 'PNT1', 'PNT1', 'PNT1','PNT1', 'PNT1', 'PNT1', 'PNT1', 'PNT1','PNT1', 'PNT1', 'PNT1', 'PNT1', 'PNT1'}';
% % gy = categorical(species);
%%
gPar = [g1; g2; g3];
gX = gPar;
species={'CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT','CNT' ...
    'PNT1', 'PNT1', 'PNT1', 'PNT1', 'PNT1','PNT1', 'PNT1', 'PNT1', 'PNT1', 'PNT1','PNT1', 'PNT1', 'PNT1', 'PNT1', 'PNT1' ...
    'PNT2', 'PNT2', 'PNT2', 'PNT2', 'PNT2','PNT2', 'PNT2', 'PNT2', 'PNT2', 'PNT2','PNT2', 'PNT2', 'PNT2', 'PNT2', 'PNT2'}';
gy = categorical(species);
%%  Feature weights by Diagonal adaptation of neighborhood component analysis (NCA)
mdl_fw=zeros(100,size(gX,2));
for i =1:length(mdl_fw)
mdl = fscnca(gX,gy,'Solver','sgd','Verbose',1);
mdl_fw(i,:)=mdl.FeatureWeights;
end
%% %
clc
close all
figure()
plot(mean(mdl_fw,1),'ro') %plot(mdl.FeatureWeights,'ro')
grid on; xlabel('Feature index'); ylabel('Feature weight')
%%%
figure ()
fr = sort(mean(mdl_fw,1),'descend');
stem(fr,'bo')
grid on; xlabel('Feature index'); ylabel('Feature weight')
title('Feature Ranking');
xlim([0 11])
g3 = [cc_g3 pc_g3 sw_g3 pl_g3 le_g3 ge_g3 md_g3 bc_g3 as_g3 corr_g1];
xticklabels({'', 'CC', 'SW', 'LE', 'GE', 'Ass', 'PL','PC', 'FC', 'Mod', 'BC'})
set(gca,'FontSize', 24, 'FontName','Calibri', 'FontWeight', 'bold')
