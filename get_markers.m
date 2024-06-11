function [T] = get_markers(Xnm_use,gene_names,cluster_idxs_use,cluster_names_use,global_enrichment_thresh,local_enrichment_thresh,fraction_expressing_threshold,write_to_file)

% Identify variable genes
[S CV_beta CV_input h1 h2] = get_single_cell_Vscores(Xnm_use, [],'fit_CVeff',true,'show_plot',true); close; close;
[sortedS,sortingIndices] = sort(S,'descend');
i0=find(~isnan(sortedS),1,'first');
topVarGenes = sortingIndices(i0:end);
[f x] = histcounts(log10(S),100);
mode_S = x(find(f==max(f),1,'first'));
ix_low = log10(S)<=mode_S;
sigma_low_S = sqrt(sum((log10(S(ix_low))-mode_S).^2)/(sum(ix_low)-1));
thresh_S = 10.^(mode_S + 2*sigma_low_S);
num_varGenes = sum(S>thresh_S*0.1);
topVarGenes = topVarGenes(1:num_varGenes); % indices for top num_varGenes
Xnm_var = Xnm_use(topVarGenes,:);
gene_names_top = gene_names(topVarGenes);
% Calculate the average and fraction expression per node for var genes
avg_node_expr = []; % genes vs clusters, avg expression level
fraction_node_expr = []; % genes vs clusters, fraction of cells with detectable expression
for i = 1:max(cluster_idxs_use)
    cell_idx = cluster_idxs_use==i;
    N = sum(cell_idx);
    tmp = Xnm_var(:,cell_idx);
    avg_expr = sum(tmp,2)./N;
    avg_node_expr = [avg_node_expr avg_expr];
    fraction_expr = sum(tmp>0.02,2)./N;
    fraction_node_expr = [fraction_node_expr fraction_expr];
end

% Filter by 'good marker' thresholds
candidate_markers = zeros(size(avg_node_expr,1),max(cluster_idxs_use));
global_enrichment = zeros(size(avg_node_expr,1),max(cluster_idxs_use));
local_enrichment = zeros(size(avg_node_expr,1),max(cluster_idxs_use));
for i = 1:size(avg_node_expr,1)
    avg_temp = avg_node_expr(i,:);
    fraction_temp = fraction_node_expr(i,:);  
    for j = 1:max(cluster_idxs_use)
        a = avg_temp(j);
        b = avg_temp;
        b(j) = [];
        b = sum(b)/length(b);
        global_enrichment(i,j) = a/b;
        global_enriched = (a/b)>global_enrichment_thresh;
        c = sort(avg_temp,'descend');
        local_enrichment(i,j) = a/c(2);
        local_enriched = (a/c(2))>local_enrichment_thresh;
        fraction_enriched = fraction_temp(j)>fraction_expressing_threshold;
        good_marker = global_enriched & local_enriched & fraction_enriched;
        if (match(gene_names_top(i),'ndrg1'))
          disp([global_enriched local_enriched fraction_enriched]);  
        end    
        candidate_markers(i,j) = good_marker;
    end
end

% Patch gaps by relaxing parameters 2-fold
for j = 1:max(cluster_idxs_use)
    if sum(candidate_markers(:,j))==0
        disp(['Relaxing parameters to try and find marker genes for ' cluster_names_use{j}])
    end
end

for i = 1:size(avg_node_expr,1)
    avg_temp = avg_node_expr(i,:);
    fraction_temp = fraction_node_expr(i,:);
    for j = 1:max(cluster_idxs_use)
        if sum(candidate_markers(:,j))==0
            global_enriched = global_enrichment(i,j)>(0.5*global_enrichment_thresh);
            local_enriched = local_enrichment(i,j)>(0.5*local_enrichment_thresh);
            fraction_enriched = fraction_temp(j)>fraction_expressing_threshold;
            good_marker = global_enriched & local_enriched & fraction_enriched;
            candidate_markers(i,j) = good_marker;
        end
    end
end

% Split results summary by cluster, sort by local enrichment, and write outputs
clear node_summary_stats marker_gene_idx
for i = 1:max(cluster_idxs_use)
    idx = logical(candidate_markers(:,i));
    if sum(idx) == 0
        disp(['Warning: no marker genes found for cluster "' num2str(i) ' = ' cluster_names_use{i} '". Consider relaxing parameters.'])
    end
    marker_gene_idx = topVarGenes(idx);
    tmp = [marker_gene_idx local_enrichment(idx,i) global_enrichment(idx,i) avg_node_expr(idx,i) fraction_node_expr(idx,i)];
    tmp = sortrows(tmp,-2);
    
    marker_gene_outputs{i}{1} = cluster_names_use{i};
    marker_gene_outputs{i}{2} = gene_names(tmp(:,1));
    marker_gene_outputs{i}{3} = tmp;
end

states = {};
marker_genes = {};
local_enrichment = [];
global_enrichment = [];
avg_expr = [];
fraction_expr = [];

for i = 1:length(marker_gene_outputs)
    count = length(marker_gene_outputs{i}{2});
    for j = 1:count
        states = [states; cluster_names_use{i}]; 
    end
    marker_genes = [marker_genes; marker_gene_outputs{i}{2}];
    local_enrichment = [local_enrichment; marker_gene_outputs{i}{3}(:,2)];
    global_enrichment = [global_enrichment; marker_gene_outputs{i}{3}(:,3)];
    avg_expr = [avg_expr; marker_gene_outputs{i}{3}(:,4)];
    fraction_expr = [fraction_expr; marker_gene_outputs{i}{3}(:,5)];
end

local_enrichment = num2cell(local_enrichment);
global_enrichment = num2cell(global_enrichment);
avg_expr = num2cell(avg_expr);
fraction_expr = num2cell(fraction_expr);

% Store in structure
s = struct('State',states,'Marker_genes',marker_genes,'Local_enrichment',local_enrichment,'Global_enrichment',global_enrichment,'Avg_expression',avg_expr,'Fraction_expression',fraction_expr);

% Convert to table
T = struct2table(s);

% Write table to csv file
if write_to_file
    writetable(T,'marker_genes.txt')
end

end
