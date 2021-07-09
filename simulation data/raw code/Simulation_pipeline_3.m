clear;
clc;


result_dir='..\results';

% import data
load([result_dir,'\','result_all_cell.mat']);
cellnum=size(result_all_cell);
cellnum=cellnum(2);

time_matrix_average=[];
time_matrix_se=[];
TF_dynamics_matrix_average=[];
TF_dynamics_matrix_se=[];
unspliced_matrix_average=[];
unspliced_matrix_se=[];
spliced_matrix_average=[];
spliced_matrix_se=[];

cor_mat=[];

for i_cell=1:cellnum
    result_i_cell=result_all_cell{i_cell};
    target_num=size(result_i_cell.time_matrix);
    target_num=target_num(2);
    
    time_matrix_average=[time_matrix_average,mean(result_i_cell.time_matrix,2)];
    time_matrix_se=[time_matrix_se,std(result_i_cell.time_matrix,0,2)/sqrt(target_num)];
    
    TF_dynamics_matrix_average=[TF_dynamics_matrix_average,mean(result_i_cell.TF_dynamics_matrix,2)];
    TF_dynamics_matrix_se=[TF_dynamics_matrix_se,std(result_i_cell.TF_dynamics_matrix,0,2)/sqrt(target_num)];
    
    unspliced_matrix_average=[unspliced_matrix_average,mean(result_i_cell.unspliced_matrix,2)];
    unspliced_matrix_se=[unspliced_matrix_se,std(result_i_cell.unspliced_matrix,0,2)/sqrt(target_num)];
    
    spliced_matrix_average=[spliced_matrix_average,mean(result_i_cell.spliced_matrix,2)];
    spliced_matrix_se=[spliced_matrix_se,std(result_i_cell.spliced_matrix,0,2)/sqrt(target_num)];
    
    cor_mat=[cor_mat,[corr(mean(result_i_cell.TF_dynamics_matrix,2),mean(result_i_cell.unspliced_matrix,2));corr(mean(result_i_cell.TF_dynamics_matrix,2),mean(result_i_cell.spliced_matrix,2))]];
    
    fprintf('%d\n',i_cell);
end

csvwrite([result_dir,'\','time_matrix_average','.csv'],time_matrix_average);
csvwrite([result_dir,'\','time_matrix_se','.csv'],time_matrix_se);
csvwrite([result_dir,'\','TF_dynamics_matrix_average','.csv'],TF_dynamics_matrix_average);
csvwrite([result_dir,'\','TF_dynamics_matrix_se','.csv'],TF_dynamics_matrix_se);
csvwrite([result_dir,'\','unspliced_matrix_average','.csv'],unspliced_matrix_average);
csvwrite([result_dir,'\','unspliced_matrix_se','.csv'],unspliced_matrix_se);
csvwrite([result_dir,'\','spliced_matrix_average','.csv'],spliced_matrix_average);
csvwrite([result_dir,'\','spliced_matrix_se','.csv'],spliced_matrix_se);
xlswrite([result_dir,'\','cor_mat','.xls'],[{'unspliced';'spliced'},num2cell(cor_mat)]);
