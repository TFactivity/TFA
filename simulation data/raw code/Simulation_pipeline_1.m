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

rand('seed',2021012502);
alpha_vector=alpha_min+rand(target_num)*(alpha_max-alpha_min);
splicing_time_vector=splicing_time_min+rand(target_num)*(splicing_time_max-splicing_time_min);
kd_A_vector=kd_A_min+rand(target_num)*(kd_A_max-kd_A_min);

simulated_cell_num=100;
result_all_cell={};

% ~500 min
h=waitbar(0,'Simulating');
for i_cell=1:simulated_cell_num
    time_matrix=[];
    TF_dynamics_matrix=[];
    unspliced_matrix=[];
    spliced_matrix=[];
    cor_matrix={'target name','half_life','corr_TF_unspliced','corr_TF_spliced','corr_unspliced_spliced'};
    
    for i_target=1:target_num
        
        target_i=target_name{i_target};
        half_life_i=target_half_life(i_target);
        
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
        
        cor_matrix=[cor_matrix;...
            {target_i,half_life_i,corr(result{1}.TF,result{1}.unspliced),...
            corr(result{1}.TF,result{1}.spliced),corr(result{1}.unspliced,result{1}.spliced)}];
        
        result_all_cell{i_cell}.time_matrix=time_matrix;
        result_all_cell{i_cell}.TF_dynamics_matrix=TF_dynamics_matrix;
        result_all_cell{i_cell}.unspliced_matrix=unspliced_matrix;
        result_all_cell{i_cell}.spliced_matrix=spliced_matrix;
        result_all_cell{i_cell}.cor_matrix=cor_matrix;
        
        if i_cell==1 && ismember(target_i,target_chosen_for_plot)
            figure('Name','Results');
            subplot(3,1,1)
            plot(result{1}.time,result{1}.TF);
            xlim([min(result{1}.time),max(result{1}.time)]);
            xlabel('time (h)','Fontname', 'Arial','FontSize',15);
            ylabel('TF dynamics (uM)','Fontname', 'Arial','FontSize',15);
            title([ 'target:' target_i ', half life (min):' num2str(half_life_i)],'Fontname', 'Arial','FontSize',15);
            subplot(3,1,2)
            plot(result{1}.time,result{1}.unspliced);
            xlim([min(result{1}.time),max(result{1}.time)]);
            xlabel('time (h)','Fontname', 'Arial','FontSize',15);
            ylabel('unspliced target (uM)','Fontname', 'Arial','FontSize',15);
            subplot(3,1,3)
            plot(result{1}.time,result{1}.spliced);
            xlim([min(result{1}.time),max(result{1}.time)]);
            xlabel('time (h)','Fontname', 'Arial','FontSize',15);
            ylabel('spliced target (uM)','Fontname', 'Arial','FontSize',15);
            
            set(gcf,'Position',[500,300,500,500])
            saveas(gcf,[result_dir,'\',sample_name,'_finalResults','.pdf']);
            close(gcf);
        end
    end
    str=['Simulating...',num2str(i_cell/simulated_cell_num*100),'%'];
    waitbar(i_cell/simulated_cell_num,h,str);
    delete('*.mat');
end
close(h);

save([result_dir,'\','result_all_cell.mat'],'result_all_cell');

load([result_dir,'\','result_all_cell.mat']);
% combine correlation data of all cells
corr_TF_unspliced_matrix=[];
corr_TF_spliced_matrix=[];
corr_unspliced_spliced_matrix=[];
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
end

corr_TF_unspliced_mean=mean(corr_TF_unspliced_matrix);
corr_TF_unspliced_se=std(corr_TF_unspliced_matrix)/sqrt(simulated_cell_num);

corr_TF_spliced_mean=mean(corr_TF_spliced_matrix);
corr_TF_spliced_se=std(corr_TF_spliced_matrix)/sqrt(simulated_cell_num);

corr_unspliced_spliced_mean=mean(corr_unspliced_spliced_matrix);
corr_unspliced_spliced_se=std(corr_unspliced_spliced_matrix)/sqrt(simulated_cell_num);

corr_combined=[target_name_vector;target_half_life_vector;
    num2cell(corr_TF_unspliced_mean);num2cell(corr_TF_unspliced_se);
    num2cell(corr_TF_spliced_mean);num2cell(corr_TF_spliced_se);
    num2cell(corr_unspliced_spliced_mean);num2cell(corr_unspliced_spliced_se)];

rownames={'name';'half life';'corr_TF_unspliced_mean';'corr_TF_unspliced_se';
    'corr_TF_spliced_mean';'corr_TF_spliced_se';'corr_unspliced_spliced_mean';'corr_unspliced_spliced_se'};
corr_combined=[rownames,corr_combined];
xlswrite([result_dir,'\','corr_combined','.xls'],corr_combined);

% save parameters
paras = ([alpha_vector(1:target_num);splicing_time_vector(1:target_num);kd_A_vector(1:target_num)])';
rownames = ['alpha','splice_time','kd'];
csvwrite([result_dir,'\','alpha_splicetime_kd','.csv'],paras);
