function [y_pred, optimize_iter_SF_act, y, best_view, d_P_all, iter_P] = Cluster(dataname, X, y, cluster_num, views_num, isGraph, select_best_view_style, anchor_num, k ,optimize_iter_P, optimize_iter_SF, alpha, delta, disrupt_index_all)


%% Find the most reliable view
% best_view_all = zeros(20,1);
% for i = 1:20
%     [best_view] = select_best_view(X, y, views_num, cluster_num, select_best_view_style);
%     best_view_all(i) = best_view;
% end

if strcmp(dataname,'NH_csmsc')
    best_view = 1;
elseif strcmp(dataname,'Caltech101-7')
    best_view = 6;
elseif strcmp(dataname,'Pascal')
    best_view = 2;
else
    [best_view] = select_best_view(X, y, views_num, cluster_num, select_best_view_style);
end



%% Normalization
for v = 1:views_num
    X{v} = zscore(X{v});
end



%% Randomly disrupted data
[X, y] = disrupt_data(X, y, views_num, best_view, disrupt_index_all);



%% Generate bipartite graphs
if isGraph == 1
    B = X;
else
    n = size(X{1},1);
    B = cell(views_num,1);

    [anchor] = gen_anchor_score(X,anchor_num,views_num);



    %% Initialization of graphs
    for v = 1:views_num
        D_v = L2_distance_1(X{v}', anchor{v}');
        [~, index] = sort(D_v, 2);
        B{v} = zeros(n,anchor_num);
        for j = 1:n
            index_11 = index(j,1:k+1);
            d_11 = D_v(j, index_11);
            B{v}(j,index_11) = (d_11(k+1)-d_11)/(k*d_11(k+1)-sum(d_11(1:k))+eps);
        end
    end
end



%% Optimization of bipartite graphs
B_cv0 = B{best_view};
[n,m] = size(B{1});
B_all = zeros(n,m*views_num);
for i = 1:views_num
    B_all(:,m*(i-1)+1:m*i) = B{i};
end

[n,m] = size(B_all);
B_all_new = B_all;
optimize_iter_SF_act = zeros(optimize_iter_P,1);
d_P_all = zeros(optimize_iter_P,1);
P_all_0 = zeros(n,n*views_num);

for iter_P = 1:optimize_iter_P

    %% Optimize S,F
    [y_pred,~,S,optimize_iter_SF_act] = cluster_optimize(B_all_new, cluster_num, optimize_iter_SF, views_num, optimize_iter_SF_act, iter_P);

    %% Optimize P
    I_n = eye(n);
    I_m = eye(m/views_num);
    for i = 1:views_num
        if i == best_view
            P_all(:,(i-1)*n+1:i*n) = eye(n);
        else
            B_all = sparse(B_all);
            B_v = B_all(:,(i-1)*(m/views_num)+1:i*m/views_num);
            S_v = S(:,(i-1)*(m/views_num)+1:i*m/views_num);
            P_v = (alpha*B_cv0-S_v) * B_v' * (I_n/delta-(alpha+1)/delta^2*B_v/((alpha+1)/delta*(B_v'*B_v)+I_m)*B_v'); 
            P_all(:,(i-1)*n+1:i*n) = P_v;
            B_all_new(:,(i-1)*(m/views_num)+1:i*m/views_num) = P_v * B_all(:,(i-1)*(m/views_num)+1:i*m/views_num);
        end
    end



    %% Normalization
    B_all_new = max(B_all_new,0);
    a = sum(B_all_new,2);
    B_all_new = views_num*B_all_new./a;



    % Calculate the change in P and jump out of the loop
    d_P = norm(P_all-P_all_0,2);
    d_P_all(iter_P) = d_P;
    d_P_change = [];
    for i = 1:iter_P-1
        d_P_change = [d_P_change,abs((d_P-d_P_all(i))/d_P_all(i))];
    end
%     fprintf('%.4f\t',d_P_change)
    if d_P < 10^-3
        break;
    elseif min(d_P_change) < 0.03
        break;
    end
    P_all_0 = P_all;
end