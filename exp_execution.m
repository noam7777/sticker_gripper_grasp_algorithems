% %% complete poteto experiment
% SO = [0,20,40,60,80,96];
% l_init = 40.0;
% s_tot = 37;
% gamma = 1;
% n_film = 50;
% peal_vel_model = 1;
% full_peal_animation = 1;
% plot_all_streamlines = 0;
% plot_optimal_sticker_state = 1;
% path_to_exp_csv_folder = "C:\Users\noam_post\OneDrive - post.bgu.ac.il\matlab\sticker analysis6_numerical_analysis-comp10277\exp_results_csv_files";
% n_seeds = 20;
% sigma = 4;
% csv_files_list = ["p_so0.csv",'p_so20.csv','p_so40.csv','p_so60.csv','p_so80.csv','p_so94.csv'];
% csv_path_list = strcat(path_to_exp_csv_folder,'\',csv_files_list);
% save_figures_flag = 1;
% save_figure_path = 'C:\Users\noam_post\OneDrive - post.bgu.ac.il\matlab\sticker analysis6_numerical_analysis-comp10277\figures\poteto';
% for ii = 1:length(SO)
%     clf
%     sO = SO(ii);
%     path_to_exp_csv = csv_path_list(ii);
%     max_load = numerical_analisis_poteto(sO,l_init,s_tot,gamma,n_film,peal_vel_model,full_peal_animation,path_to_exp_csv,plot_all_streamlines,plot_optimal_sticker_state,sigma,n_seeds);
%     if save_figures_flag == 1 
%         saveas(gcf,strcat(save_figure_path,'\','poteto_so',num2str(sO),'.png'))
%     end
%     disp("maximum load of potato at sO = " + sO + " is: " + max_load + "[N]");
% end

% %% complete experiment peeked_cylinder
% SO = [2,40,80,120,160,180];
% l_init = 35.36;
% s_tot = 85;
% gamma = 2;
% n_film = 30;  
% peal_vel_model = 1;
% full_peal_animation =0;
% plot_all_streamlines = 0;
% plot_optimal_sticker_state = 1;
% n_seeds = 20;
% sigma = 4;
% path_to_exp_csv_folder = "C:\Users\noam_post\OneDrive - post.bgu.ac.il\matlab\sticker analysis6_numerical_analysis-comp10277\exp_results_csv_files";
% csv_files_list = ["pc_so0.csv",'pc_so40.csv','pc_so80.csv','pc_so120.csv','pc_so160.csv','pc_so180.csv'];
% csv_path_list = strcat(path_to_exp_csv_folder,'\',csv_files_list);
% 
% save_figures_flag = 1;
% save_figure_path = 'C:\Users\noam_post\OneDrive - post.bgu.ac.il\matlab\sticker analysis6_numerical_analysis-comp10277\figures\peaked_cylinder';
% 
% for ii = 1:length(SO)
%     clf
%     sO = SO(ii);
%     path_to_exp_csv = csv_path_list(ii);
%     max_load = numerical_analisis_peeked_cylinder(sO,l_init,s_tot,gamma,n_film,peal_vel_model,full_peal_animation,path_to_exp_csv,plot_all_streamlines,plot_optimal_sticker_state,,sigma,n_seeds);
%     if save_figures_flag == 1 
%         saveas(gcf,strcat(save_figure_path,'\','peaked_cylinder_so',num2str(sO),'.png'))
%     end
%     disp("maximum load of peeked_cylinder at sO = " + sO + " is: " + max_load + "[N]");
% end

% %% complete experiment crazy_square
% SO = [0,40,77,120,150,170]; % 150 on sim is 160 inreality. 170 on sim is 180 in reality
% l_init = 80.36;
% s_tot = 60;
% gamma = 1;
% n_film =40;
% peal_vel_model = 1; %1 for the v_dot proprtional to mc_cr. 2 for simple A or B pealing. 
% full_peal_animation = 1;
% n_seeds = 20;
% sigma = 2; 
% path_to_exp_csv_folder = "C:\Users\noam_post\OneDrive - post.bgu.ac.il\matlab\sticker analysis6_numerical_analysis-comp10277\exp_results_csv_files";
% csv_files_list = ["cs_so0.csv",'cs_so40.csv','cs_so80.csv','cs_so120.csv','cs_so160.csv','cs_so180.csv'];
% csv_path_list = strcat(path_to_exp_csv_folder,'\',csv_files_list);
% plot_all_streamlines = 0;
% plot_optimal_sticker_state = 1;
% 
% save_figures_flag = 1;
% save_figure_path = 'C:\Users\noam_post\OneDrive - post.bgu.ac.il\matlab\sticker analysis6_numerical_analysis-comp10277\figures\crazy_square';
% 
% for ii = 1:length(SO)
% 
%     clf
%     sO = SO(ii);
%     path_to_exp_csv = csv_path_list(ii);
%     max_load = numerical_analisis_crazy_square(sO,l_init,s_tot,gamma,n_film,peal_vel_model,full_peal_animation,path_to_exp_csv,plot_all_streamlines,plot_optimal_sticker_state,,sigma,n_seeds);
%     if save_figures_flag == 1 
%         saveas(gcf,strcat(save_figure_path,'\','crazy_square_so',num2str(sO),'.png'))
%     end
%     disp("maximum load of crazy_square at sO = " + sO + " is: " + max_load + "[N]");
% end

%% compare the peal models

sO = 20;
l_init = 40.0;
s_tot = 37;
gamma = 1;
full_peal_animation = 0;
plot_all_streamlines = 0;
plot_optimal_sticker_state = 0;
path_to_exp_csv_folder = "C:\Users\noam_post\OneDrive - post.bgu.ac.il\matlab\sticker analysis6_numerical_analysis-comp10277\exp_results_csv_files";
csv_files_list = 'p_so20.csv';
save_figures_flag = 1;
sigma = 0;
n_seeds = 1;

save_figure_path = 'C:\Users\noam_post\OneDrive - post.bgu.ac.il\matlab\sticker analysis6_numerical_analysis-comp10277\figures\compare_peal_models';
path_to_exp_csv = strcat(path_to_exp_csv_folder,'\',csv_files_list); % <-------
n_film_array = 50:5:800;


% model 1 proportional:
peal_vel_model = 1;
max_load_array = zeros(length(n_film_array),1);
for ii = 1:length(n_film_array)
    n_film = n_film_array(ii);
    clf
    max_load = numerical_analisis_poteto(sO,l_init,s_tot,gamma,n_film,peal_vel_model,full_peal_animation,path_to_exp_csv,plot_all_streamlines,plot_optimal_sticker_state,sigma,n_seeds);
    max_load_array(ii) = max_load;
    disp("maximum load of potato at sO = " + num2str(sO) + " is: " + num2str(max_load) + "[N], n_film = " + num2str(n_film) + 'peal model = '+ num2str(peal_vel_model));
end
max_load_array_model1 = max_load_array;
% model 2 descreate:

peal_vel_model = 2;
max_load_array = zeros(length(n_film_array),1);
for ii = 1:length(n_film_array)
    n_film = n_film_array(ii);
    clf
    max_load = numerical_analisis_poteto(sO,l_init,s_tot,gamma,n_film,peal_vel_model,full_peal_animation,path_to_exp_csv,plot_all_streamlines,plot_optimal_sticker_state,sigma,n_seeds);
    max_load_array(ii) = max_load;
    disp("maximum load of potato at sO = " + num2str(sO) + " is: " + num2str(max_load) + "[N], n_film = " + num2str(n_film) + 'peal model = '+ num2str(peal_vel_model));
end
max_load_array_model2 = max_load_array;

figure(2)
plot(n_film_array,max_load_array_model1,n_film_array,max_load_array_model2)
xlabel("Discretization Resolution",'interpreter','latex')
ylabel('Grasp Value ${\left[\frac{N}{J/m^2}\right]}$','interpreter','latex')
% xlabel('${\mu}$','interpreter','latex'
legend('Peal Model 1', 'Peal Model 2')
if save_figures_flag == 1 
    saveas(gcf,strcat(save_figure_path,'\','poteto_so',num2str(sO),'grasp_value_vs_resolution.png'))
end

figure(3)
error_model1 = abs((max_load_array_model1-max_load_array_model1(end))./max_load_array_model1(end)).*100;
error_model2 = abs((max_load_array_model2-max_load_array_model2(end))./max_load_array_model2(end)).*100;
plot(n_film_array(1:end-1),error_model1(1:end-1),n_film_array(1:end-1),error_model2(1:end-1))
xlabel("Discretization Resolution",'interpreter','latex')
ylabel('Error [%]','interpreter','latex')
legend('Peal Model 1', 'Peal Model 2')
if save_figures_flag == 1 
    saveas(gcf,strcat(save_figure_path,'\','poteto_so',num2str(sO),'error_vs_resolution.png'))
end
