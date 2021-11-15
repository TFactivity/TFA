
%%0
clear;
clc;

result_dir='..\results';
mkdir(result_dir);

% import P53 target genes and their half lives
P53_target=readtable(['.\p53 target half-lives\P53_target_half_life.csv']);
target_name=table2array(P53_target(:,2));
target_half_life=table2array(P53_target(:,3));
target_num=length(target_name);

target_chosen_for_plot={target_name{1},target_name{round((target_num+1)/2)},target_name{target_num}};

alpha_min=100; %unit: uM/h
alpha_max=200; %unit: uM/h

splicing_time_min=5; %unit: min
splicing_time_max=10; %unit: min

kd_A_min=1/10;
kd_A_max=1/2;

alpha_vector=power(10,0:0.5:3);
alpha_vector=100;
n_len = numel(alpha_vector);
splicing_time_vector=zeros(numel(alpha_vector),1)+splicing_time_max;
kd_A_vector=zeros(numel(alpha_vector),1)+kd_A_max;

simulated_cell_num=5;
result_all_cell={};

% ~500 min
for i_cell=1:simulated_cell_num
    time_matrix=[];
    TF_dynamics_matrix=[];
    unspliced_matrix=[];
    spliced_matrix=[];
    cor_matrix={'target name','half_life','corr_TF_unspliced','corr_TF_spliced','corr_unspliced_spliced','corr_TF_unspliced_spar','corr_TF_spliced_spar'};
    
    for i_target=1:numel(alpha_vector)
        
        target_i=target_name{i_target};
        half_life_i=100;
        
        alpha_0=alpha_vector(i_target); %unit: uM/h
        splicing_time=splicing_time_vector(i_target); %unit: min
        beta_0=log(10)/(splicing_time/60); %unit: /h
        
        half_life=half_life_i; %unit: min
        
        gamma_0=log(2)/(half_life/60); %unit: /h
        n_0=2;
        
        A_basal_0=0.06; %unit: uM
        A_max_0=0.5; %unit: uM
        T_0=5.5; %unit: h
        sigma_0=A_basal_0/5; %unit: uM
        
        kd_0=A_max_0*kd_A_vector(i_target); %unit: uM
        
        tau_0=0.001; %unit: h
        timelimit_0=20; %unit: h
        timelag_record_0=0; %unit: h
        
        cellnum_0=1;
        sample_name=target_i;
        result_dir=result_dir;
        
        result=P53_decoding_model(alpha_0,beta_0,gamma_0,n_0,kd_0,...
            A_basal_0,A_max_0,T_0,sigma_0,...
            tau_0,timelimit_0,timelag_record_0,cellnum_0,...
            sample_name,result_dir);
        
        time_matrix=[time_matrix result{1}.time];
        TF_dynamics_matrix=[TF_dynamics_matrix result{1}.TF];
        unspliced_matrix=[unspliced_matrix result{1}.unspliced];
        spliced_matrix=[spliced_matrix result{1}.spliced];
        
        % simulate dropout
        p = 0.1;
        un_mat_spar = zeros(size(result{1}.spliced,1),size(result{1}.spliced,2));
        sp_mat_spar = zeros(size(result{1}.spliced,1),size(result{1}.spliced,2));
        for i = 1:size(result{1}.spliced,1)
            for j = 1:size(result{1}.spliced,2)
                un_mat_spar(i,j) = sim_sparsity(result{1}.unspliced(i,j),p);
                sp_mat_spar(i,j) = sim_sparsity(result{1}.spliced(i,j),p);
            end
        end
        
        
        cor_matrix=[cor_matrix;...
            {target_i,half_life_i,corr(result{1}.TF,result{1}.unspliced),...
            corr(result{1}.TF,result{1}.spliced),corr(result{1}.unspliced,result{1}.spliced),...
            corr(result{1}.TF,un_mat_spar),corr(result{1}.TF,sp_mat_spar)}];
        
        result_all_cell{i_cell}.time_matrix=time_matrix;
        result_all_cell{i_cell}.TF_dynamics_matrix=TF_dynamics_matrix;
        result_all_cell{i_cell}.unspliced_matrix=unspliced_matrix;
        result_all_cell{i_cell}.spliced_matrix=spliced_matrix;
        result_all_cell{i_cell}.cor_matrix=cor_matrix;
        
        
    end
    str=['Simulating...',num2str(i_cell/simulated_cell_num*100),'%'];
    delete('*.mat');
    
    fprintf('%f\n',i_cell/simulated_cell_num);
end

save([result_dir,'\','result_sparsity.mat'],'result_all_cell');


load([result_dir,'\','result_sparsity.mat']);
%% calculate sparsity matrix
% simulate dropout
ps = 0.05:0.05:0.3;
cor_un_spar = zeros(numel(ps)*numel(ps),simulated_cell_num);
cor_sp_spar = zeros(numel(ps)*numel(ps),simulated_cell_num);
ps_mat = zeros(numel(ps)*numel(ps),2);
for cell_id = 1:simulated_cell_num
    un_mat= result_all_cell{1,cell_id}.unspliced_matrix;
    sp_mat= result_all_cell{1,cell_id}.spliced_matrix;
    TF_mat = result_all_cell{1,cell_id}.TF_dynamics_matrix;
    for p1_id = 1:numel(ps)
        for p2_id = 1:numel(ps)
            p_id= (p1_id-1)*numel(ps)+p2_id;
            p1 = ps(p1_id);
            p2 = ps(p2_id);
            ps_mat(p_id,:) = [p1,p2];
            
            un_mat_spar = un_mat;
            sp_mat_spar = sp_mat;
            for i = 1:size(un_mat,1)
                for j = 1:size(un_mat,2)
                    un_mat_spar(i,j) = sim_sparsity(un_mat(i,j),p1);
                    sp_mat_spar(i,j) = sim_sparsity(sp_mat(i,j),p2);
                end
            end   
            
            cor_un_spar(p_id,cell_id) = corr(un_mat_spar,TF_mat);
            cor_sp_spar(p_id,cell_id) = corr(sp_mat_spar,TF_mat);

        end
    end
end
% calculate mean correlation for all cells
cor_un_spar_mean = mean(cor_un_spar,2);
cor_sp_spar_mean = mean(cor_sp_spar,2);
cor_ratio = cor_un_spar_mean/cor_sp_spar_mean;
plot3(ps_mat(:,1),ps_mat(:,2),cor_un_spar_mean,'.r');
hold on
plot3(ps_mat(:,1),ps_mat(:,2),cor_sp_spar_mean,'.b');
xlabel('p_{unspliced}'),ylabel('p_{spliced}');

% save data as csv file
csvwrite('cor_un.csv',cor_un_spar_mean);
csvwrite('cor_sp.csv',cor_sp_spar_mean);
csvwrite('pinfo.csv',ps_mat);
        
%% 2
% combine correlation data of all cells
corr_TF_unspliced_matrix=[];
corr_TF_spliced_matrix=[];
corr_unspliced_spliced_matrix=[];
corr_TF_unspliced_matrix_spar= [];
corr_TF_spliced_matrix_spar = [];
for i_cell=1:simulated_cell_num
    cor_matrix_i_cell=result_all_cell{i_cell}.cor_matrix;
    cor_matrix_i_cell(1,:)=[];
    if i_cell==1
        target_name_vector=cor_matrix_i_cell(:,1)';
        target_half_life_vector=cor_matrix_i_cell(:,2)';
    end
    
    corr_TF_unspliced_matrix=[corr_TF_unspliced_matrix;(cell2mat(cor_matrix_i_cell(:,3)))'];
    corr_TF_spliced_matrix=[corr_TF_spliced_matrix;(cell2mat(cor_matrix_i_cell(:,4)))'];
    corr_unspliced_spliced_matrix=[corr_unspliced_spliced_matrix;(cell2mat(cor_matrix_i_cell(:,5)))'];
    corr_TF_unspliced_matrix_spar=[corr_TF_unspliced_matrix_spar;(cell2mat(cor_matrix_i_cell(:,6)))'];
    corr_TF_spliced_matrix_spar=[corr_TF_spliced_matrix_spar;(cell2mat(cor_matrix_i_cell(:,7)))'];
end

corr_TF_unspliced_mean=mean(corr_TF_unspliced_matrix);
corr_TF_unspliced_se=std(corr_TF_unspliced_matrix)/sqrt(simulated_cell_num);

corr_TF_spliced_mean=mean(corr_TF_spliced_matrix);
corr_TF_spliced_se=std(corr_TF_spliced_matrix)/sqrt(simulated_cell_num);

corr_unspliced_spliced_mean=mean(corr_unspliced_spliced_matrix);
corr_unspliced_spliced_se=std(corr_unspliced_spliced_matrix)/sqrt(simulated_cell_num);

corr_TF_unspliced_spar_mean=mean(corr_TF_unspliced_matrix_spar);
corr_TF_spliced_spar_mean=mean(corr_TF_spliced_matrix_spar);

% semilogx(alpha_vector, corr_TF_unspliced_mean,'r');
% hold on
% semilogx(alpha_vector, corr_TF_spliced_mean,'b');
% xlabel('alpha'), ylabel('correlation with TF');
% 
% unspliced_m = mean(unspliced_matrix,1);
% spliced_m = mean(spliced_matrix,1);
% 
% semilogx(spliced_m, corr_TF_unspliced_mean,'r');
% hold on
% semilogx(spliced_m, corr_TF_spliced_mean,'b');
% xlabel('mean count'), ylabel('correlation with TF');


Y = [corr_TF_unspliced_mean,corr_TF_unspliced_spar_mean;corr_TF_spliced_mean,corr_TF_spliced_spar_mean];
h = bar(Y,0.75);
ylabel('correlation with TF');
set(gca,'XTickLabel',{'unspliced','spliced','dropout: unspliced','dropout: spliced'});
legend('No subsampling','Subsampling with p=0.1');
%% simulate dropout
% function for dropout, input value -> output value
p = 0.1;
input = x;
y = 0;
for i = 1:x
    tmp = rand;
    if tmp < p
        y = y+1;
    end
end
sim_sparsity(10,0.1)

