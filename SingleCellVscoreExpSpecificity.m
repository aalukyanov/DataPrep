function ClutchSpec = ...
               SingleCellVscoreExpSpecificity()

%
% Load Mapping to Human Symbols 
%

isUseBlast                     = 2;    % Use or not BLAST gene symbols assigment ?
isGlobalBACIQ                  = false; 
MS_detrend_method              = 8;    % Normalization factor
isUniqueOnly                   = false;
%
isFiltering_69_Mapping         = -1;
SBH_69_Eps_Mapping             = 1.e-50;
RBH_69_Eps_Mapping             = 1.e-50;
isFROG2Xelaevs_Mapping         = -1;
SBH_FR2Xel_Eps_Mapping         = 1.e-50;
RBH_FR2Xel_Eps_Mapping         = 1.e-50;
%//
isGeneSymbolAssigmentFiltering = -1; 
SBH_GeneSAF_Eps_Mapping        = 1.e-50;
RBH_GeneSAF_Eps_Mapping        = 1.e-50;      
%
Map                            = Mapping (isUseBlast,                      ...
                                          isFiltering_69_Mapping,          ...
                                          SBH_69_Eps_Mapping,              ...
                                          RBH_69_Eps_Mapping,              ...
                                          isFROG2Xelaevs_Mapping,          ...
                                          SBH_FR2Xel_Eps_Mapping,          ...
                                          RBH_FR2Xel_Eps_Mapping,          ...
                                          isGeneSymbolAssigmentFiltering,  ... 
                                          SBH_GeneSAF_Eps_Mapping,         ...
                                          RBH_GeneSAF_Eps_Mapping);
%
%
%
BName1                         = Map.UnFilteredBName;
BProtein                       = Map.UnFilteredBProtein;
BProteinID                     = Map.UnFilteredBProteinID;
BDescription                   = Map.UnFilteredBDescription;
BSBH                           = Map.UnFilteredBSBH;
BRBH                           = Map.UnFilteredBRBH;
%
%
% Load Mapping
TableNameLevisTropicalis   = char('input/LevisTropicalisMapping.csv');
Table_User_LevisTropicalis = readtable(TableNameLevisTropicalis, 'ReadVariableNames', false);

GeneLevNames               = table2array(Table_User_LevisTropicalis(1:end, 'Var1'));
GeneTropNames              = table2array(Table_User_LevisTropicalis(1:end, 'Var2'));
SBH_LevTrop_Map            = table2array(Table_User_LevisTropicalis(1:end, 'Var5'));
RBH_LevTrop_Map            = table2array(Table_User_LevisTropicalis(1:end, 'Var6'));


EpsValue                   = 1e-10; 

FilteredIndMappingLevTrop  = (SBH_LevTrop_Map < EpsValue) & (RBH_LevTrop_Map  < EpsValue); 

FGeneLevNames              = GeneLevNames(FilteredIndMappingLevTrop);
FGeneTropNames             = GeneTropNames(FilteredIndMappingLevTrop);
FSBH_LevTrop_Map           = SBH_LevTrop_Map(FilteredIndMappingLevTrop);
FRBH_LevTrop_Map           = RBH_LevTrop_Map(FilteredIndMappingLevTrop);

%
%

PROThours    = [-2.00;   ... %// oo56 
                 0.00;   ... %// egg   Yes (included)
                 7.00;   ... %// st09  Yes (included)
                13.25;   ... %// st12  Yes (included)
                18.75;   ... %// st17
                24.00;   ... %// st22  Yes (included)
                26.25;   ... %// st24
                29.50;   ... %// st26
                35.00;   ... %// st30  Yes (included)
                80.00;   ... %// st42
                360      ... %// st50
             ];

%
%
%cd /Users/peshaG4/Research/SingleCell/Xenopus/tree_browser
%addpath(genpath('~/Research/SingleCell/code'))   % add path where I have needed routines like tot_normalize_sc_gene_counts
%

load   gene_names.mat;
%//load A_Xnm_all.mat

load   Xnm_node_avgs.mat;
load   tp_map.mat;

%
load   X22.mat; load X20.mat; load X18.mat; load X16.mat; load X14.mat; load X12.mat; %load  X10.mat
%//Xall = [X22 X20 X18 X16 X14 X12];
Xall    = X22;

clear X22 X20 X18 X16 X14 X12; 
%//X22 = Xnm_node_avgs; % select stage 22 data 

prcLoFrac    = 0.15;  % percentile of proteins in low variability 
SigFracCutOf = 0.05;

%//

[XNrm, umicounts]            = tot_normalize_sc_gene_counts(Xall); %  normalize to total expression per cell; 
[V, CV_beta, CV_input, ~, ~] = get_single_cell_Vscores(XNrm, [],'fit_CVeff',true,'show_plot', true); close; close;

%//save  saveVscore V CV_beta CV_input

%//
%// This is the original code
%// [X22Nrm, umicounts]          = tot_normalize_sc_gene_counts(X22); %  normalize to total expression per cell; 
%// [V, CV_beta, CV_input, ~, ~] = get_single_cell_Vscores(X22Nrm, [],'fit_CVeff',true,'show_plot', true); close; close;
%//

load saveVscore.mat

Nans        = isnan(V);
V(isnan(V)) = median(V(~isnan(V)));          % set NaN to 1 not to get in the way of sorting 
[i,j] = sort(V,'descend'); i(1:10), j(1:10); % get ten most specific genes and check 

%

SpecGeneNames  = gene_names(j(1:10));        % specific names
LSpecGeneNames = gene_names(j(end-10:end));  % less specific names

%
%
%

%figure; 
%  hist(log10(V(IndV)),71);
%  xlabel('V score (log10 units)'); ylabel('counts')  % plot histogram of V-score 

%//  
%//figure;
%//  plot(log10(sum(Xall,2)), sum(Xall<=0.5,2),'ro');
%//  xlabel('Comulative expression in all cells  (log10 units)'); 
%//  ylabel('#cells expressing ');  
%//  IndexG  = find(match(gene_names,'mpo'));
%//  X       = log10(sum(Xall,2));
%//  Y       = sum(Xall<=0.5,2);
%//  hold on;
%//  plot(X(IndexG),      Y(IndexG), 'b*');
%//  text(X(IndexG)*1.01, Y(IndexG)*1.01, 'mpo', 'FontSize',16);
%//
 
%gname(gene_names)  %plot the over-dispersion scatterplot 


%
% (1) Find all good values without NaN
%

IndV                       =  ~isnan(V);
VV                         =  V(IndV);
VGene_names                =  gene_names(IndV);

[IntGenes, IA, IJ]         =  intersect(FGeneTropNames, VGene_names);
SFVGeneLevNames            =  FGeneLevNames(IA); 
SFVGene_names              =  VGene_names(IJ);
SFVGeneTropNames           =  FGeneTropNames(IA);
SFV                        =  V(IJ);

%
% Find intersect with the experimentally measured
% Load Mapping
%
TableNameLevisTropicalis   = char('input/LevisTropicalisMapping.csv');
Table_User_LevisTropicalis = readtable(TableNameLevisTropicalis, 'ReadVariableNames', false);
%
%
TableNameProteinsData      = char('input/user_proteins_out_case_final.csv');
Table_User_Proteins        = readtable(TableNameProteinsData, 'ReadVariableNames', true);
% Data
ProtSym                    = table2array(Table_User_Proteins(1:end, 'ProteinID'));
NumPeps                    = table2array(Table_User_Proteins(1:end, 'NumPeps'));
CProteinID                 = table2array(Table_User_Proteins(1:end, 'ProteinID'));
UData                      = table2array(Table_User_Proteins(:, {'x126_Me',  ...
                                                                 'x127n_Me', ...
                                                                 'x127c_Me', ...
                                                                 'x128n_Me', ...
                                                                 'x128c_Me', ...
                                                                 'x129n_Me', ...
                                                                 'x129c_Me', ...
                                                                 'x130n_Me', ...
                                                                 'x130c_Me', ...
                                                                 'x131_Me'}));                                                               
%
% Overlap with protein experimental data
%

[IALL_LevNames, IA_A, JB_A]  = intersect(CProteinID, SFVGeneLevNames, 'stable');

%
%

ClutchSpec.FVGeneLevNames      =  SFVGeneLevNames  (JB_A); 
ClutchSpec.FVGene_names        =  SFVGene_names    (JB_A);
ClutchSpec.FVGeneTropNames     =  SFVGeneTropNames (JB_A);
ClutchSpec.FV                  =  SFV              (JB_A);
ClutchSpec.FUData              =  UData(IA_A,:);
ClutchSpec.FCProteinID         =  CProteinID(IA_A);

%                                                                                                                      
%













