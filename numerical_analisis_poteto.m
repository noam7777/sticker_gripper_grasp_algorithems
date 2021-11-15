function max_load = numerical_analisis_poteto(sO,l_init,s_tot,gamma,n_film,peal_vel_model,full_peal_animation,path_to_exp_csv,plot_all_streamlines,plot_optimal_sticker_state,sigma,n_seeds)

hold off
%% bodies definition
%bodies are defined counter clockwise!!!!
%peeked cilynder
     n=1000;
    t = linspace(0,2*pi,n);
    t=t(1:n-1);
    b = pi;
    c = 0.3;
    g = 0.4*(exp(-(t).^2/c^2)+exp(-(t-2*b).^2/c^2)+exp(-(t-b).^2/c^2));
    r = 100*(1+g);
    x = cos(t+pi).*r;
    y = sin(t+pi).*r;
    body1 = body(x,y,0,0);    
%poteto:
    n=1000;
    t = linspace(0,2*pi,n);
    t=t(1:n-1);
    r = 50.*(cos(t).^2+0.4.*sin(t).^2);
    % r = 1.5.*cos(t).^2+0.4.*sin(t).^2;
    x = cos(-t+pi).*50.*(1.5.*cos(t).^2+0.4.*sin(t).^2);
    y = sin(t+pi).*50.*(1.5.*cos(t).^2+0.4.*sin(t).^2);
    body2 = body(x,y,0,0);
    



    
    
% Dodecagon % !!!!!!!!!!!!!!!!!!!!!!!!TODO
    %t = linspace(0,2*pi,4);
    %t(end) = [];
    %x = cos(t)+0.5*cos(2*t);
    %y = sin(t);
    x = [2,1.9, -1.9,-2,-2, -2,-2,-1.9, 1.9,2,2, 2];
    y = [2,2  2,2,1.9, -1.9,-2,-2, -2,-2,-1.9, 1.9];
    body3 = body(x,y,0,0);
bodies = [body1,body2,body3];

body_i = bodies(2);
%% grasp definition

% sO = 170;
%     l_init = 1;
%     s_tot = 50;
%     gamma = 1;
%     n_film = 16;
    filmA = film(l_init,s_tot,gamma,n_film);
    filmB = filmA;
    film_direction = 1;
gr11 = grasp_single_film(body_i,sO,filmA,film_direction,s_tot); %for film A, the counter clockwise filme assuming that the body is defined ccw
gr12 = grasp_single_film(body_i,sO,filmB,film_direction*(-1),s_tot); %for film B
gr2 = grasp2(body_i,sO,filmA,filmB,s_tot);


%% numerical analisis parameters
% n_seeds = 20;
% sigma = 4; % related to the distribution of the sticker's initial state at the grasp point
seeds = abs(normrnd(0,sigma,n_seeds,2));

%% calculate the peeling velocity

if peal_vel_model == 1 %% streamlines , version1: velocity direction proportional to the
%% devision between mg_crA/mg_crB.    
    s_dot = (gr2.mg_cr_A.^2+gr2.mg_cr_B.^2).^0.5;
    s_dot_A = gr2.mg_cr_B./s_dot;
    s_dot_B = gr2.mg_cr_A./s_dot;

   
%% streamlines , version2: velocity direction is where mg_cr is smaller.
elseif peal_vel_model == 2
    s_dot = min(gr2.mg_cr_A,gr2.mg_cr_B);
    s_dot_A = gr2.mg_cr_A<gr2.mg_cr_B;
    s_dot_B = gr2.mg_cr_B<gr2.mg_cr_A;
else 
    disp("choose a valid peal_vel_model (1 for the simple one direction or 2 for v_dot proportional to mg_cr")
end
    
% this part is true for both methods of s_dot estimation:    
s_dot(gr2.LA_loose_cond | gr2.LB_loose_cond) = 1;
s_dot_A(gr2.LA_loose_cond) = 0;
s_dot_A(gr2.LB_loose_cond) = 1;
s_dot_B(gr2.LB_loose_cond) = 0;
s_dot_B(gr2.LA_loose_cond) = 1;
%% estimating the quality of the grasp
streamlines_samples = stream2(gr2.SB,gr2.SA,s_dot_B,s_dot_A,seeds(:,1),seeds(:,2));
mg_cr_max_samples = zeros(n_seeds,1);
end_points = zeros(n_seeds,2);
pealing_index_max_gr2 = zeros(n_seeds,1);
mg_cr_gr2 = zeros(n_seeds,1);
mg_cr_gr1 = zeros(n_seeds,1);
pealing_index_max_gr1 = zeros(n_seeds,1);
optimal_state_ufter_1_film_failure_flag = zeros(n_seeds,1);
sA_optimal = zeros(n_seeds,1); 
sB_optimal = zeros(n_seeds,1); 
filmA_rip_1st = boolean(zeros(n_seeds,1)); %condition array. did film A failed 1st?
entrance_index = zeros(n_seeds,1); %represents the index on the single film
                                   %plot where the streamline enters it. 
                                   %it holds for both filmA and filmB cases
for ii = 1:n_seeds %run on each stereamline sample
    strm = cell2mat(streamlines_samples(ii));
    [mg_cr_gr2(ii),pealing_index_max_gr2(ii)] = max(interp2(gr2.SB,gr2.SA,gr2.mg_cr,strm(:,1),strm(:,2)));
    % chaking the maximum value of the streamline extention on the 1-film
    
    % pealing_index_max_gr2 is the index along the streamline where mgcr is maximal
    % we keep sa and sb of the state where mg_cr is maximal:
    sB_optimal(ii) = strm(pealing_index_max_gr2(ii),1);
    sA_optimal(ii) = strm(pealing_index_max_gr2(ii),2);
    
    % grasp curve:
    end_point = strm(end,:);% end pOint is [sB_end,sA_end] of the streamline
                            % it is uniqe to the streamline.
    end_points(ii,:) = end_point; %for later use..
    % in this case, s_tot is equal for both filmA and filmB.
    s_totA = s_tot;
    s_totB = s_tot;
    filmA_rip_1st(ii) = (s_totA - end_point(2))<(s_totB - end_point(1));
    if filmA_rip_1st(ii)
        entrance_index(ii) = round(n_film*(end_point(1)/s_totB));
        [mg_cr_gr1(ii),pealing_index_max_gr1(ii)] = max(gr12.mg_cr(entrance_index(ii):end));
        if mg_cr_gr1(ii) > mg_cr_gr2 %if the best state of the single grasp ufter the failiour of the other is grater then the double film optimal grasp, the optimal state is updated
            optimal_state_ufter_1_film_failure_flag(ii) = 1; % true if the optimal grasp state is on a curve of a single grasp ufter the other one has failed
            sB_optimal(ii) = gr12.sE(pealing_index_max_gr1(ii)+entrance_index(ii)-1);
            sA_optimal(ii) = nan;
        end
    else
        entrance_index(ii) = round(n_film*(end_point(2)/s_totA));
        [mg_cr_gr1(ii),pealing_index_max_gr1(ii)] = max(gr11.mg_cr(entrance_index(ii):end));
        if mg_cr_gr1(ii) > mg_cr_gr2 %if the best state of the single grasp ufter the failiour of the other is grater then the double film optimal grasp, the optimal state is updated
            optimal_state_ufter_1_film_failure_flag(ii) = 1; % true if the optimal grasp state is on a curve of a single grasp ufter the other one has failed
            sA_optimal(ii) = gr11.sE(pealing_index_max_gr1(ii)+entrance_index(ii)-1);
            sB_optimal(ii) = nan;
        end
    end
    if mg_cr_gr1(ii) > mg_cr_gr2
        optimal_state_ufter_1_film_failure_flag(ii) = 1; % true if the optimal grasp state is on a curve of a single grasp ufter the other one has failed
    end
    mg_cr_max_samples(ii)  = max(mg_cr_gr1(ii),mg_cr_gr2(ii)); %the maximum over all the course of the streamline
end

mg_max_grasp = mean(mg_cr_max_samples); %finally, we estimate the reward value for this grasp
max_load = mg_max_grasp;



%% plots of the body, the griper positions, and the mg_cr plots
figure(1)
subplot(3,6,[1,2,3,7,8,9,13,14,15])
print(body_i); 
hold on
% scatter(gr11.E(:,1),gr11.E(:,2),'k');
% scatter(gr11.Dx(:),gr11.Dy(:),1,'m');
% scatter(gr12.E(:,1),gr12.E(:,2),'k');
% scatter(gr12.Dx(:),gr12.Dy(:),1,'m');
scatter(gr2.A(:,1),gr2.A(:,2),'g');
scatter(gr2.B(:,1),gr2.B(:,2),'b');
% scatter(gr2.Dx(:),gr2.Dy(:),1,'k');
% xlim([-7 7])
% ylim([-7 7])
body_i.plot_body_point(sO);


if full_peal_animation == 1
    for ii = 1:40:size(strm,1) % we animate the last streamline created
        sA_k = strm(ii,2);
        sB_k = strm(ii,1);
        [~,sA_k_index] = min(abs(sA_k-gr2.sA));
        [~,sB_k_index] = min(abs(sB_k-gr2.sB));
        if gr2.LA_loose_cond(sA_k_index,sB_k_index) == 1
            gr12.plot_sticker_at_state(sB_k)
        elseif gr2.LB_loose_cond(sA_k_index,sB_k_index) == 1
            gr11.plot_sticker_at_state(sA_k)
        else
            gr2.plot_sticker_at_state(sA_k,sB_k)
        end
        pause(0.1);
    end
end

% plot the optimal state 
if plot_optimal_sticker_state == 1
    if optimal_state_ufter_1_film_failure_flag(1) == 0% we enter this if the optimal state acuure ufter the first filme hase failed. we do this only for the first seed
        [~,sA_optimal_index] = min(abs(sA_optimal(1)-gr2.sA));
        [~,sB_optimal_index] = min(abs(sB_optimal(1)-gr2.sB));
        if gr2.LA_loose_cond(sA_optimal_index,sB_optimal_index) == 1
            gr12.plot_sticker_at_state_thick(sB_optimal)
        elseif gr2.LB_loose_cond(sA_optimal_index,sB_optimal_index) == 1
            gr11.plot_sticker_at_state_thick(sA_optimal(1))
        else
            gr2.plot_sticker_at_state_thick(sA_optimal(1),sB_optimal(1))
        end
    else % plot the optimal state of a single film ufter the other one has failed
        if filmA_rip_1st(1) 
            gr12.plot_sticker_at_state_thick(gr12.sE(pealing_index_max_gr1(1)+entrance_index(1)-1))
        else
            gr11.plot_sticker_at_state_thick(gr11.sE(pealing_index_max_gr1(1)+entrance_index(1)-1))
            % we need the sB optimal for the single film grasp
        end
    end
end
axis equal

%% plotting mg_cr in the sAsB plane in simulation and experiment
subplot(3,6,[10,11,16,17])
contourf(gr2.sB,gr2.sA,gr2.mg_cr,80, 'LineColor', 'none')
xlabel('$s_B$','interpreter','latex')%x axiss is the colomn index - meaning B direction
ylabel('$s_A$','interpreter','latex')%y axiss corespond to row - A direction
zlabel('$mg_{critical}$','interpreter','latex')%y axiss corespond to row - A direction
hp4 = get(subplot(3,6,[10,11,16,17]),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.15  hp4(2)  0.05  hp4(2)+hp4(3)*3])
hold on

if plot_all_streamlines == 1
    streamline(gr2.SB,gr2.SA,s_dot_B,s_dot_A,gr2.SA,gr2.SB)
end
str = streamline(streamlines_samples);
set(str,'color','red','LineWidth',3)
scatter(seeds(:,1),seeds(:,2),'k') % marking the seeds
%quiver(gr2.SB,gr2.SA,s_dot_B,s_dot_A)

if plot_optimal_sticker_state == 1
    for ii = 1:n_seeds %run on each stereamline sample
        if optimal_state_ufter_1_film_failure_flag(ii) == 0
            scatter(sB_optimal,sA_optimal,200,'m*') 
        end
    end
end

EXP = readtable(path_to_exp_csv);
sA1 = EXP.triel1_sA_left;
sB1 = EXP.triel1_sB_right;
sA2 = EXP.triel2_sA_left;
sB2 = EXP.triel2_sB_right;
sA3 = EXP.triel3_sA_left;
sB3 = EXP.triel3_sB_right;
plot(sB1,sA1,'y--o',sB2,sA2,'y--o',sB3,sA3,'y--o')
axis equal

    
%% display the single film grasp curve
%in case filmB was disattached first, we plot the sA vs mg_cr curve
subplot(3,6,[12,18])
hold on
plot(gr11.mg_cr,gr11.sE) 
ylabel('$s_A$','interpreter','latex')%x axiss is the colomn index - meaning B direction
xlabel('$mg_{critical}$','interpreter','latex')%y axiss corespond to mg_cr - A direction
B_curve_entrance = entrance_index(~filmA_rip_1st);
scatter(gr11.mg_cr(B_curve_entrance),gr11.sE(B_curve_entrance))
plot(gr11.mg_cr(min(B_curve_entrance):end),gr11.sE(min(B_curve_entrance):end),'r','LineWidth',1)
if plot_optimal_sticker_state == 1
    for ii = 1:n_seeds %run on each stereamline sample
        if optimal_state_ufter_1_film_failure_flag(ii) == 1 && filmA_rip_1st(ii) == 0
            scatter(mg_cr_gr1(ii),sA_optimal(ii),200,'m*') 
        end
    end
end
hold off

%in case filmA was disattached first, we plot the sB vs mg_cr curve
subplot(3,6,[4,5])
hold on
plot(gr12.sE,gr12.mg_cr)
xlabel('$s_B$','interpreter','latex')%x axiss is the colomn index - meaning B direction
ylabel('$mg_{critical}$','interpreter','latex')%y axiss corespond to mg_cr - A direction
A_curve_entrance = entrance_index(filmA_rip_1st);
scatter(gr12.sE(A_curve_entrance),gr12.mg_cr(A_curve_entrance))
plot(gr12.sE(min(A_curve_entrance):end),gr12.mg_cr(min(A_curve_entrance):end),'r','LineWidth',1)
if plot_optimal_sticker_state == 1
    for ii = 1:n_seeds %run on each stereamline sample
        if optimal_state_ufter_1_film_failure_flag(ii) == 1 && filmA_rip_1st(ii) == 1
            scatter(sB_optimal(ii),mg_cr_gr1(ii),200,'m*') 
        end
    end
end
hold off
axis equal

end