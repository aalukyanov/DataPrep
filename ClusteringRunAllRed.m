function ClusteringRunAllRed()  
%

clear all;

ClusterRun = 0; 


Input_PROThours       = [-6.00;   ... %// oo56 
                          0.00;   ... %// egg   Yes (included)
                          9;      ... %// st09  Yes (included)
                          12;     ... %// st12  Yes (included)
                          17;     ... %// st17
                          22;     ... %// st22  Yes (included)
                          24;     ... %// st24
                          26;     ... %// st26
                          30;     ... %// st30  Yes (included)
                          42;     ... %// st42
                          50      ... %// st50
                        ];


%//Input_PROThours    = [-2.00;   ... %// oo56 
%//                       0.00;   ... %// egg   Yes (included)
%//                       7.00;   ... %// st09  Yes (included)
%//                      13.25;   ... %// st12  Yes (included)
%//                      18.75;   ... %// st17
%//                      24.00;   ... %// st22  Yes (included)
%//                      26.25;   ... %// st24
%//                      29.50;   ... %// st26
%//                      35.00;   ... %// st30  Yes (included)
%//                      80.00;   ... %// st42
%//                      360      ... %// st50
%//                        ];
                                                                   
NoCluster          = 2;             

%//Number_Cluster =  [24; ... %      1
%//                   16; ... %Done  2
%//                   4;  ... %Done  3
%//                   8;  ... %Done  4
%//                   8;  ... %Done  5
%//                   8;  ... %Done  6
%//                  12;  ... %Done  7
%//                  32;  ... %Done  8
%//                  12;  ... %Done  9
%//                  12;  ... %Done 10
%//                   3;  ... %Done 11
%//                   6;  ... %Done 12
%//                   4;  ... %Done 13
%//                  24;  ... %Done 14
%//                  24;  ... %Done 15
%//                  32;  ... %Done 16
%//                   4;  ... %Done 17
%//                  24;  ... %Done 18
%//                  12;  ... %Done 19
%//                  24];     %Done 20      
%//                        
                                    
%
%
% Load Mapping for tissue specificity computev

TableNameLevisTropicalis   = char('data/LevisTropicalisMapping.csv');
Table_User_LevisTropicalis = readtable(TableNameLevisTropicalis, 'ReadVariableNames', false);

GeneLevNames               = table2array(Table_User_LevisTropicalis(1:end, 'Var1'));
GeneTropNames              = table2array(Table_User_LevisTropicalis(1:end, 'Var2'));
SBH_LevTrop_Map            = table2array(Table_User_LevisTropicalis(1:end, 'Var5'));
TRBH_LevTrop_Map           = table2array(Table_User_LevisTropicalis(1:end, 'Var6'));

RBH_LevTrop_Map            = str2double(TRBH_LevTrop_Map);

EpsValue                   = 1e-10; 

FilteredIndMappingLevTrop  = (SBH_LevTrop_Map < EpsValue) & (RBH_LevTrop_Map  < EpsValue); 

FGeneLevNames              = GeneLevNames(FilteredIndMappingLevTrop);
FGeneTropNames             = GeneTropNames(FilteredIndMappingLevTrop);
FSBH_LevTrop_Map           = SBH_LevTrop_Map(FilteredIndMappingLevTrop);
FRBH_LevTrop_Map           = RBH_LevTrop_Map(FilteredIndMappingLevTrop);

%
%
%
%                        
                                             
ClustInfoInputTable   = readtable('geneset/genesetList.txt', 'ReadVariableNames', false);
AllSets               = table2array(ClustInfoInputTable(1:end, 'Var1')); % Set names 
Number_Cluster        = table2array(ClustInfoInputTable(1:end, 'Var2')); % Number of clusters                    
   
Cluster_Discription   = table2array(ClustInfoInputTable(1:end, 'Var3'));
Cluster_Discription   = upper(Cluster_Discription);

Cluster_NamesID       = table2array(ClustInfoInputTable(1:end, 'Var4'));

% Input_Title           = AllSets{NoCluster};   
% Input_Number_Cluster  = Number_Cluster(NoCluster); 
% 
% Input_Title_F = fullfile('geneset', strcat(Input_Title,'.txt')); 
% T_Table       = readtable(Input_Title_F, 'ReadVariableNames', false);
% Input_List    = table2array(T_Table(1:end, 'Var1'));
% Input_List    = upper(Input_List);

N = 10; % Select how many points to consider

%
% Load appropriate protein file
%

ProtTableAll          = readtable('data/user_proteins_out_case_final.txt', 'ReadVariableNames', true);
ProteinStructure.ID   = table2array(ProtTableAll(1:end, 'ProteinID')); % ID

ProteinStructure.Data   = table2array(ProtTableAll(:, {'x126_Me',     ... % 
                                                       'x127n_Me',    ... %
                                                       'x127c_Me',    ... %
                                                       'x128n_Me',    ... %
                                                       'x128c_Me',    ... %
                                                       'x129n_Me',    ... %
                                                       'x129c_Me',    ... %
                                                       'x130n_Me',    ... %
                                                       'x130c_Me',    ...
                                                       'x131_Me'})) ;     %
% Let's normalized Data                                                             
ProteinStructure.Data(:,1:N) = ProteinStructure.Data(:,1:N)./sum(ProteinStructure.Data(:,1:N),2);
 
%//%
%//% Let's remove the phospho treated piptides
%//%
%// 
%//PhosphoTableAll            = readtable('data/FinalPhosphoCluster.txt', 'ReadVariableNames', true);
%//ProteinStructure.PhosphoID = table2array(PhosphoTableAll(1:end, 'GeneID')); % ID
%//%
%//IndFiltlProtP              = find(~match(ProteinStructure.ID,ProteinStructure.PhosphoID));
%//ProteinStructure.ID        = ProteinStructure.ID(IndFiltlProtP,:);
%//ProteinStructure.Data      = ProteinStructure.Data(IndFiltlProtP,:);
%//

%
% Get all necessary mRNAs names and IDs from Xenbase 
%

mRNATableAll          =  readtable('data/mRNA_Data_New.txt', 'ReadVariableNames', false);

iafmRNA_ID            =  table2array(mRNATableAll(1:end, 'Var1'));
iaGene_Names_mRNA     =  table2array(mRNATableAll(1:end, 'Var2'));

%
% Check if the names are in the correct format
%

for k=1:size(ProteinStructure.ID,1) %
  name            = ProteinStructure.ID{k};
  uns_s           = find(match(name,'|'));
  uns_f           = find(match(name,'.'));
  if (~isempty(uns_s) && ~isempty(uns_f))
    ProteinStructure.ID(k) = cellstr(string(name(uns_s(1)+1:uns_f-1))); % only protein name
  end  
end
  
%
%
%

isFiltering_69_Mapping           = -1;
SBH_69_Eps_Mapping               = 1.e-50;
RBH_69_Eps_Mapping               = 1.e-50;
isFROG2Xelaevs_Mapping           = -1;
SBH_FR2Xel_Eps_Mapping           = 1.e-50;
RBH_FR2Xel_Eps_Mapping           = 1.e-50;
%//
isGeneSymbolAssigmentFiltering   = -1; 
SBH_GeneSAF_Eps_Mapping          = 1.e-50;
RBH_GeneSAF_Eps_Mapping          = 1.e-50; 

isUseBlast    = 2;    % Use which mapping case:
                      %      (0) This is BLAST file
                      %      (1) This is old HMMER file
                      %      (2) This is new HMMER file
%
Map           = Mapping (isUseBlast,                      ...
                         isFiltering_69_Mapping,          ...
                         SBH_69_Eps_Mapping,              ...
                         RBH_69_Eps_Mapping,              ...
                         isFROG2Xelaevs_Mapping,          ...
                         SBH_FR2Xel_Eps_Mapping,          ...
                         RBH_FR2Xel_Eps_Mapping,          ...
                         isGeneSymbolAssigmentFiltering,  ... 
                         SBH_GeneSAF_Eps_Mapping,         ...
                         RBH_GeneSAF_Eps_Mapping);    

UFBProtein                      = Map.UnFilteredBProtein;
UFBProteinID                    = Map.UnFilteredBProteinID;

NotSelectUnFilteredBName        = Map.NotSelectUnFilteredBName;
NotSelectUnFilteredBProtein     = Map.NotSelectUnFilteredBProtein;
NotSelectUnFilteredBProteinID   = Map.NotSelectUnFilteredBProteinID;
NotSelectUnFilteredBDescription = Map.NotSelectUnFilteredBDescription;
NotSelectUnFilteredBSBH         = Map.NotSelectUnFilteredBSBH;
NotSelectUnFilteredBRBH         = Map.NotSelectUnFilteredBRBH;  

%
% Stage 50 analysis
%

TableNameProteinsData   = char('Input/Xen_MS3_17_18.csv');
Table_User_ProteinsMS2  = readtable(TableNameProteinsData, 'ReadVariableNames', true);

SProteinID              = table2array(Table_User_ProteinsMS2(1:end, 'ProteinID'));
l                       = 1;

for k=1:size(SProteinID,1) %
  name            = SProteinID{k};
  uns_s           = find(match(name,'#'));
  if (isempty(uns_s)) 
    %
    ProteinID{l,1} = name;
    Stages(l,:)    = table2array(Table_User_ProteinsMS2(k, {'x126_Me',  ... % 
                                                            'x127n_Me', ... %
                                                            'x127c_Me', ... %
                                                            'x128n_Me', ... %
                                                            'x128c_Me', ... %
                                                            'x129n_Me', ... %
                                                            'x129c_Me', ... %
                                                            'x130n_Me', ... %
                                                            'x130c_Me', ... %
                                                            'x131_Me',  ... %
                                                            'x131C_Me'}));  %
    %
    l = l + 1;
    %
  end
end
                                                             
% Let's normalized Data                                                             
Stages(:,1:11)      = Stages(:,1:11)./sum(Stages(:,1:11),2);

Stage42             = Stages(:,10);
Stage50             = Stages(:,11);

%
% Set up the color pattern
%

% define some colors and plotting options 
plt.MS_clusters  = 1; 
plt.Geometry     = [200,400,590,270]; 
plt.dump         = 0;
plt.Gray         = [.4 .4 .4]; 
plt.GrayLight    = [.6 .6 .6]; 
plt.Cyan         = [.0 1. 1.]; 
plt.Blue         = [.04 .14 .98];   % MATHEMATICA  blue 
plt.Green        = [.16 1. .18];    % MATHEMATICA green
plt.GreenDark    = [.26 .58 0.17];  % Dark Green
plt.Magenta      = [1. 0. 1.];      % Magenta 
plt.Khaki        = [.52 .38 .12];   % Khaki 

%
% Let's chek up data between all cluctes and 17&18 cases  
%

[VComm, VIAs, VIBs] = intersect(ProteinStructure.ID, ProteinID); 

VIDName42  = ProteinStructure.ID(VIAs);
VData42    = ProteinStructure.Data(VIAs,10);
VStage42   = Stages(VIBs,10);
VStage50   = Stages(VIBs,11);

%
% Let's find out human gene symbols (for intersect of ALL & 17 & 18)
%

[FCommIDs, FIAs, FIBs]  = intersect(UFBProteinID,VIDName42);
FProteinSym             = UFBProtein(FIAs);
FVData42                = VData42(FIBs);
FVStage42               = VStage42(FIBs);
FStage50                = VStage50(FIBs);

%//%
%//% Comparison between Xenbase and 42 Stage (17 & 18 data sets)
%//%
%//
%//figure
%//  % 
%//  scatter(FVData42, FVStage42, [], plt.Gray);
%//  hold on;
%//  plot([0 1], [0 1], '--r');
%//  gname(FProteinSym');
%//  hold on;
%//  %
%//  title(char('Comparision of Stage42 (ALL) and Stage42 (17&18)')); 
%//  %
%//  ylabel({'','Stage42 relative abundance'});
%//  % Create ylabel
%//  xlabel({'Stage42 relative abundance';''});
%//  % 
%//  set(gca,'XLim',[0, 1]);
%//  set(gca,'YLim',[0, 1]);
%//  %  
%//  set(gca,'FontSize',17);
%//


%
% Let's have here a comparison with Xenbase data
%  
  
FiltIndex        = FStage50 ./ FVData42 > 5;
FNameEnrStg50    = FProteinSym(FiltIndex);
  
figure
  % 
  scatter(FVData42, FStage50, [], plt.Gray);
  hold on;
  plot([0 1], [0 1], '--r');
  hold on;
  %
  title(char('Comparision of Stage42 (ALL) and Stage50 (17&18)')); 
  %
  ylabel({'','Stage50 relative abundance'});
  % Create ylabel
  xlabel({'Stage42 relative abundance (ALL)';''});
  % 
  set(gca,'XLim',[0, 1]);
  set(gca,'YLim',[0, 1]);
  %  
  set(gca,'FontSize',17);
  hold on;
  %  
  %//gname(FProteinSym');

%
% Let's find out human gene symbols (17 & 18 only)
%

[CommIDs, IAs, IBs]     = intersect(UFBProteinID,ProteinID);
ProteinSym              = UFBProtein(IAs);
NStage42                = Stage42(IBs);
NStage50                = Stage50(IBs);

%
%
%

TEST             = {'FETUB'; 'H1FOO'; 'HBG1'; 'TUBB6'};
                                  
Index            = find(match(UFBProtein,TEST));
FiltIndex        = NStage50 ./ NStage42 > 5;
NameEnrStg50     = ProteinSym(FiltIndex);

%
%
%

figure
  % 
  scatter(NStage42, NStage50, [], plt.Gray);
  hold on;
  plot([0 1], [0 1], '--r');
  hold on;
  %
  title(char('Comparision of Stage42 and Stage50')); 
  %
  ylabel({'','Stage50 relative abundance'});
  % Create ylabel
  xlabel({'Stage42 relative abundance';''});
  % 
  set(gca,'XLim',[0, 0.6]);
  set(gca,'YLim',[0, 1]);
  %  
  set(gca,'FontSize',17);
  %
  %//gname(ProteinSym');
  % 

%
% Let's load up trully xenbase data names
%

TableXenbase2JGImap     = readtable('data/Xenbase2JGImap_xbgene_v3.txt', 'ReadVariableNames', true);
TableXenbase2JGImap_ID  = table2array(TableXenbase2JGImap(1:end, 'PROTEINID'));
TableXenbase2JGImap_Sym = table2array(TableXenbase2JGImap(1:end, 'GENE_SYMBOL'));

%
% Read Levis to Tropicalis to Human Mapping 
%

Xla92XtrTableAll  =  readtable('data/Xla92Xtr.28920187.proteins.genesym.defs.uniq.txt', 'ReadVariableNames', false);
Xla92Xtr.ID       =  table2array(Xla92XtrTableAll(1:end, 'Var1')); % ID
Xla92Xtr.GeneSym  =  table2array(Xla92XtrTableAll(1:end, 'Var3')); % Gene Symbols     

for k=1:size(Xla92Xtr.ID,1) %
  name            = Xla92Xtr.ID{k};
  uns_s           = find(match(name,'|'));
  if (~isempty(uns_s))
    name(uns_s(2:end)) ='_';
    Xla92Xtr.ID(k) = cellstr(string(name(uns_s(1)+1:end))); % only protein name
  end  
end

%
% Since not all genes can have human IDs we need to be carefull here
% Let's loop over all output genes and assign a Human ID if you have
% something
%

 for k = 1:size(ProteinStructure.ID,1)
   %
   %
   
   fIXTMT  = find(match(UFBProteinID,ProteinStructure.ID(k)));
   fIXTMT1 = find(match(NotSelectUnFilteredBProteinID,ProteinStructure.ID(k)));
   fIXTMT2 = find(match(TableXenbase2JGImap_ID,ProteinStructure.ID(k)));
   fIXTMT3 = find(match(Xla92Xtr.ID,ProteinStructure.ID(k)));
   
   %
   % Please check that the mapping option is set to  
   %  isUseBlast = 2;    % Use which mapping case
   %
     
   % This is HMMER namming  
   if (~isempty(fIXTMT))
     FinalHumanNames(k)  = UFBProtein(fIXTMT);
     FinalProteinID(k)   = ProteinStructure.ID(k);
   else
     if (~isempty(fIXTMT3))
       FinalHumanNames(k)  = Xla92Xtr.GeneSym(fIXTMT3);
       FinalProteinID(k)   = ProteinStructure.ID(k);
     else 
       FinalHumanNames(k)  = ProteinStructure.ID(k);
       FinalProteinID(k)   = ProteinStructure.ID(k);
     end
   end
   
   % This is BLAST naming  
   if (~isempty(fIXTMT1))
     FinalHumanNamesBLAST(k) = NotSelectUnFilteredBProtein(fIXTMT1);
     FinalProteinIDBLAST(k)  = ProteinStructure.ID(k);
   else
     if (~isempty(fIXTMT3))
       FinalHumanNamesBLAST(k) = Xla92Xtr.GeneSym(fIXTMT3);
       FinalProteinID(k)       = ProteinStructure.ID(k);
     else        
       FinalHumanNamesBLAST(k) = ProteinStructure.ID(k);
       FinalProteinIDBLAST(k)  = ProteinStructure.ID(k);
     end
   end
   
   % This is Xenbase naming  
   if (~isempty(fIXTMT2))
     FinalHumanNames_Xenbase(k) = upper(TableXenbase2JGImap_Sym(fIXTMT2));
     FinalProteinID_Xenbase(k)  = ProteinStructure.ID(k);
   else
     if (~isempty(fIXTMT3))
       FinalHumanNames_Xenbase(k) = Xla92Xtr.GeneSym(fIXTMT3);
       FinalProteinID(k)          = ProteinStructure.ID(k);
     else               
       FinalHumanNames_Xenbase(k) = ProteinStructure.ID(k);
       FinalProteinID_Xenbase(k)  = ProteinStructure.ID(k);
     end
   end
   %
   %
 end

%save GenesSymbols  FinalHumanNames FinalHumanNamesBLAST FinalHumanNames_Xenbase;
%save GeneIDS       FinalProteinID  FinalProteinIDBLAST  FinalProteinID_Xenbase;

%load GenesSymbols;
%load GeneIDS;

%
%

ProteinStructure.Symbol        = transpose(FinalHumanNames);              % Symbol (HMMER symbols, please check "isUseBlast = 2")
ProteinStructure.SymbolBlast   = transpose(FinalHumanNamesBLAST);         % Symbol (BlAST symbols, please check "isUseBlast = 2")
ProteinStructure.SymbolXenbase = transpose(FinalHumanNames_Xenbase);      % Symbol (BlAST symbols, please check "isUseBlast = 2")

%
% Tissue specificity definition 
%

ProteinStructure.VTissueSpec        = strings(size(ProteinStructure.ID,1),1);
ProteinStructure.VTissueSpec(1:end) = "Unknown"; % this is initial guess 
ProteinStructure.VScore             = zeros(size(ProteinStructure.ID,1),1);

%
% Load up all the concentrations (Contians 11 channels with the last channel for the Follical Cells)
%

ConcTable2016   = readtable('data/PROT_Data_New__16_WebPage.txt', 'ReadVariableNames', false);
ProteinStructure.ConcGeneIDs2016 = table2array(ConcTable2016(:,'Var1')) ;
ProteinStructure.ConcData2016    = table2array(ConcTable2016(:, {'Var2',   ... % 
                                                                 'Var3',   ... %
                                                                 'Var4',   ... %
                                                                 'Var5',   ... %
                                                                 'Var6',   ... %
                                                                 'Var7',   ... %
                                                                 'Var8',   ... %
                                                                 'Var9',   ... %
                                                                 'Var10',  ... %
                                                                 'Var11'})) ;  %
%
%
ConcTable2017   = readtable('data/PROT_Data_New__17_WebPage.txt', 'ReadVariableNames', false);
ProteinStructure.ConcGeneIDs2017 = table2array(ConcTable2017(:,'Var1')) ;
ProteinStructure.ConcData2017    = table2array(ConcTable2017(:, {'Var2',   ... % 
                                                                 'Var3',   ... %
                                                                 'Var4',   ... %
                                                                 'Var5',   ... %
                                                                 'Var6',   ... %
                                                                 'Var7',   ... %
                                                                 'Var8',   ... %
                                                                 'Var9',   ... %          
                                                                 'Var10',  ...
                                                                 'Var11'})) ;  %                                                             
%
%
ConcTable2018   = readtable('data/PROT_Data_New__18_WebPage.txt', 'ReadVariableNames', false);
ProteinStructure.ConcGeneIDs2018 = table2array(ConcTable2018(:,'Var1')) ;
ProteinStructure.ConcData2018    = table2array(ConcTable2018(:, {'Var2',   ... % 
                                                                 'Var3',   ... %
                                                                 'Var4',   ... %
                                                                 'Var5',   ... %
                                                                 'Var6',   ... %
                                                                 'Var7',   ... %
                                                                 'Var8',   ... %
                                                                 'Var9',   ... %
                                                                 'Var10',  ... %
                                                                 'Var11'})) ;  %  
                                                                                                                                                                                                                                                

%
% Load up the cell data
%
% 

GenesCells      = [];
Celltype        = [];
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage8/expr.csv', 'ReadVariableNames', false);
ExprTableName   = table2array(ExprTable(:,'Var1'));
clear ExprTable;         
%
%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage8/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage8/cell_groupings.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(2,2:end))]; % first entry has to be skiped 
%disp("---Stage 8");
%disp(table2array(CellGroupTable(2,2:end)));

if (ClusterRun == 1)

%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage10/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage10/cell_groupings.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(2,2:end))]; % first entry has to be skiped 
%disp("---Stage 10");
%disp(table2array(CellGroupTable(2,2:end)));
%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage11/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage11/cell_groupings2.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(1,2:end))]; % first entry has to be skiped 
%disp("---Stage 11");
%disp(table2array(CellGroupTable(1,2:end)));
%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage12/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage12/cell_groupings.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(2,2:end))]; % first entry has to be skiped 
%disp("---Stage 12");
%disp(table2array(CellGroupTable(2,2:end)));
%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage13/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage13/cell_groupings.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(2,2:end))]; % first entry has to be skiped 
%disp("---Stage 13");
%disp(table2array(CellGroupTable(2,2:end)));
%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage14/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage14/cell_groupings.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(2,2:end))]; % first entry has to be skiped 
%disp("---Stage 14");
%disp(table2array(CellGroupTable(2,2:end)));
%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage16/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage16/cell_groupings.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(2,2:end))]; % first entry has to be skiped
%disp("---Stage 16");
%disp(table2array(CellGroupTable(2,2:end)));
%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage18/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage18/cell_groupings.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(2,2:end))]; % first entry has to be skiped 
%disp("---Stage 18");
%disp(table2array(CellGroupTable(2,2:end)));
%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage20/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage20/cell_groupings.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(1,2:end))]; % first entry has to be skiped 
%disp("---Stage 20");
%disp(table2array(CellGroupTable(1,2:end)));
%
ExprTable       = readtable('.././ClustByGeneSet_Final_Phos/Stage22/expr.csv', 'ReadVariableNames', false);
CellGroupTable  = readtable('.././ClustByGeneSet_Final_Phos/Stage22/cell_groupings.csv', 'ReadVariableNames', false);
GenesCells      = [GenesCells table2array(ExprTable(1:end,2:end))];
%//Celltype     = table2array(CellGroupTable(4,2:end)); % first entry has to be skiped
% Let's use the cell labels (it includes stages)
Celltype        = [Celltype table2array(CellGroupTable(2,2:end))]; % first entry has to be skiped 
%disp("---Stage 22");
%disp(table2array(CellGroupTable(2,2:end)));

end

%
% Sum Check
%

disp("-------Sum Check-------");
disp(size(Celltype,2) == size(GenesCells,2));

%
%
for k=1:size(ExprTableName,1) %
  name            = ExprTableName{k};
  uns_s           = find(match(name,'|'));
  if (~isempty(uns_s))
    ExprTableNameFinal(k,1) = cellstr(string(name(1:uns_s(1)-2)));     % only protein name
  end  
end
%
%% Total count normalize
% To make comparisons between cells with different numbers of UMIs, we
% 'total count normalize', i.e. scale the counts of each cell to the
% average of all cells in the sample.

minCt=1000; % minimum number of UMIs/cell to allow

keepCellIdx=sum(GenesCells,1)>=minCt;
[Xnm_use, tot_counts, tot_counts_orig] = tot_normalize_sc_gene_counts(GenesCells(:,keepCellIdx),true,0.05); %true/false option sets whether genes with total representation over 0.05 are excluded from the normalization calculation
[V, CV_beta, CV_input, ~, ~]           = get_single_cell_Vscores(Xnm_use, [],'fit_CVeff',true,'show_plot', true); close; close;
Nans                                   = isnan(V);
V(isnan(V))                            = median(V(~isnan(V)));                     % set NaN to 1 not to get in the way of sorting 
%
fraction                               = 0.05;
LowSpec                                = quantile(V(:), fraction);
HighSpec                               = quantile(V(:), 1 - fraction);
%
UniqueCelltype                         = unique(Celltype);
IndexCelltype                          = zeros(size(Celltype,2),1);
%
% Let's make index of all tissues (the order of tissues controled by Celltype)
%
for k = 1 : size(UniqueCelltype,2)
  IndCluster                 = find(match(Celltype,UniqueCelltype(k)));
  IndexCelltype(IndCluster)  = k;
end    
%
%
%

%% Define and rank marker genes for each cell state
global_enrichment_thresh      = 4;     % expression ratio threshold relative to rest of embryo
local_enrichment_thresh       = 0.5;   % threshold relative to next highest state
fraction_expressing_threshold = 0.7;   % threshold on number of cells with detected expression
write_to_file                 = false; % option to write outputs to .txt file in current directory

[marker_gene_table]           = get_markers(Xnm_use,                       ...
                                            ExprTableNameFinal,            ...
                                            IndexCelltype,                 ...
                                            UniqueCelltype,                ...
                                            global_enrichment_thresh,      ...
                                            local_enrichment_thresh,       ...
                                            fraction_expressing_threshold, ...
                                            write_to_file);                                     
%
%
disp("-------Number of Tissues--------");
disp(size(UniqueCelltype,2));
disp("-------Names--------");
disp(UniqueCelltype');
%
%
states              = table2array(marker_gene_table(1:end, 'State'));
marker_genes        = table2array(marker_gene_table(1:end, 'Marker_genes'));
local_enrichment    = table2array(marker_gene_table(1:end, 'Local_enrichment'));
global_enrichment   = table2array(marker_gene_table(1:end, 'Global_enrichment'));
avg_expr            = table2array(marker_gene_table(1:end, 'Avg_expression'));
fraction_expr       = table2array(marker_gene_table(1:end, 'Fraction_expression'));
%
%

%{
%
%
% Print Tables
%
fid_ind  = fopen("TableOfIndexes1.txt", 'w');
OK       = 1;
%
if (fid_ind < 0)
 OK = 0;
 disp('Problem to open file'); 
 exit;
end
%
%
%
fprintf(fid_ind, '%20s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
                 "GeneName",      ...
                 "Vscore",        ... 
                 "MI",            ...  
                 "MI_NM",         ...
                 "Gini",          ... 
                 "Gini_CLUST",    ... 
                 "Gini_CLUST_NM", ... 
                 "Min",           ...
                 "Max",           ...
                 "Min_NM",        ... 
                 "Max_NM",        ...
                 "Min_Exp",     ... 
                 "Max_Exp",     ...
                 "Num_Count",     ...
                 "Num_Cells",     ...
                 "Case4",         ...
                 "Case5"); 
%
% Compute Mutual information and Gini Index
%
for n = 1 : size(GenesCells,1)
  %
  disp(n);
  %
  % Compute mutual information
  %
  GenesCellsExpression    = GenesCells(n,1:end);
  MIG(n)                  = mi(GenesCellsExpression',IndexCelltype,size(UniqueCelltype,2));   % symmetrical information
  
  GenesCellsExpressionNM  = Xnm_use(n,1:end);
  MIG_NM(n)               = mi(GenesCellsExpressionNM',IndexCelltype, size(UniqueCelltype,2)); % symmetrical information
  
  avg_node_expr      = []; % genes vs clusters, avg expression level
  avg_node_expr_nm   = [];

  for i = 1:max(IndexCelltype)
    %  
    cell_idx         = IndexCelltype == i;
    Nc               = sum(cell_idx);
    tmp              = GenesCellsExpression(:,cell_idx);   
    avg_expr_loc     = sum(tmp,2)./Nc;
    avg_node_expr    = [avg_node_expr avg_expr_loc];
    %    
    tmp_nm           = GenesCellsExpressionNM(:,cell_idx);   
    avg_expr_nm      = sum(tmp_nm,2)./Nc;
    avg_node_expr_nm = [avg_node_expr_nm avg_expr_nm];
    %
    if (find(match(ExprTableNameFinal(n),'Xetrov90028359m'))>0)   
      disp([num2str(sum(tmp>0)) num2str(Nc)  UniqueCelltype(i) ExprTableNameFinal(n)]);
    end 
    %
  end
  %
  avg_expr_clust(n)          = sum(avg_node_expr)/size(avg_node_expr,2);
  giniCoeff(n)               = GiniCalc(GenesCellsExpression);
  giniCoeff_cluster(n)       = GiniCalc(avg_node_expr);
  giniCoeff_cluster_nm(n)    = GiniCalc(avg_node_expr_nm);
  %
  GeneExprGeneMin(n)         = min(avg_node_expr);
  GeneExprGeneMax(n)         = max(avg_node_expr);
  %
  GeneExprGeneMin_NM(n)      = min(avg_node_expr_nm);
  GeneExprGeneMax_NM(n)      = max(avg_node_expr_nm);  
  %
  Case4(n)                   = sum(avg_node_expr > 2) == size(UniqueCelltype,2);
  Case5(n)                   = sum(avg_node_expr < 2) == size(UniqueCelltype,2);
  
  MaxExp(n)                  = max(GenesCellsExpression);
  MinExp(n)                  = min(GenesCellsExpression);
  NumCount(n)                = sum(GenesCellsExpression > 0)/max(IndexCelltype);
  NumCells(n)                = sum(GenesCellsExpression > 0); 
  
  % 
  %for m = 1 : size(GenesCells,1)
  %  GenesCellsExpression_m  = GenesCells(m,1:end);
  %  CIG(n,m)                = mi(GenesCellsExpression',GenesCellsExpression_m'); % symmetrical information
  %end   
  % print
  %
  fprintf(fid_ind, '%20s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',  ...
          char(ExprTableNameFinal(n)),                                        ...
          num2str(round(V(n),4)),                                             ... 
          num2str(round(MIG(n),4)),                                           ...
          num2str(round(MIG_NM(n),4)),                                        ...
          num2str(round(giniCoeff(n),4)),                                     ... 
          num2str(round(giniCoeff_cluster(n),4)),                             ...
          num2str(round(giniCoeff_cluster_nm(n),4)),                          ...
          num2str(round(GeneExprGeneMin(n),4)),                               ...
          num2str(round(GeneExprGeneMax(n),4)),                               ...
          num2str(round(GeneExprGeneMin_NM(n),4)),                            ...
          num2str(round(GeneExprGeneMax_NM(n),4)),                            ...
          num2str(round(MinExp(n),4)),                                        ... 
          num2str(round(MaxExp(n),4)),                                        ...
          num2str(round(NumCount(n),4)),                                      ...  
          num2str(round(NumCells(n),4)),                                      ... 
          num2str(Case4(n)),                                                  ...
          num2str(Case5(n)));
  %
end

%
%save CovMI.mat CIG; 
%

disp('++++Main Structure Loop');


gl = 1;
for n = 1 : size(FGeneLevNames,1)
  %
  disp(n);
  % 
  Ind = find(match(marker_genes,FGeneTropNames(n)));
  
  if (~isempty(Ind))

   %disp(Ind)
   
   for m = 1 : size(Ind,1)
  
     Ind_MI_GI                  = find(match(ExprTableNameFinal,marker_genes(Ind(m))));

     if (GeneExprGeneMax(Ind_MI_GI) > 3)

       states_g(gl,1)             = states(Ind(m));
       marker_genes_lev_g(gl,1)   = FGeneLevNames(n);
      
       HumInd                     = find(match(UFBProteinID,FGeneLevNames(n)));
   
       if(~isempty(HumInd))
         marker_genes_hum_g(gl,1) = UFBProtein(find(match(UFBProteinID,FGeneLevNames(n))));
       else
         marker_genes_hum_g(gl,1) = FGeneLevNames(n);  
       end

       MIG_g(gl,1)                = {round(MIG(Ind_MI_GI),5)};       % symmetrical information
       MIG_NM_g(gl,1)             = {round(MIG_NM(Ind_MI_GI),5)};
     
       GI_g(gl,1)                 = {round(giniCoeff(Ind_MI_GI),5)};
       CGI_g(gl,1)                = {round(giniCoeff_cluster(Ind_MI_GI),5)};
     
       GeneExprGeneMin_g(gl,1)    = {round(GeneExprGeneMin(Ind_MI_GI),5)};
       GeneExprGeneMax_g(gl,1)    = {round(GeneExprGeneMax(Ind_MI_GI),5)};
       NumCount_g(gl,1)           = {round(NumCount(Ind_MI_GI),5)};     

       marker_genes_trop_g(gl,1)  = upper(marker_genes(Ind(m)));
       local_enrichment_g(gl,1)   = {round(local_enrichment(Ind(m)),5)};
       global_enrichment_g(gl,1)  = {round(global_enrichment(Ind(m)),5)};
       avg_expr_g(gl,1)           = {round(avg_expr(Ind(m)),5)};
       fraction_expr_g(gl,1)      = {round(fraction_expr(Ind(m)),5)};
   
       SDH_g(gl,1)                = {FSBH_LevTrop_Map(n)};
       RBH_g(gl,1)                = {FRBH_LevTrop_Map(n)};

       gl                         = gl + 1;

     end

   end

  end    
    
end    

% Store in structure
s = struct('State', states_g, ...
           'Marker_Laevis_IDs', marker_genes_lev_g, ...    
           'Marker_Laevis_Human_genes', marker_genes_hum_g, ...    
           'Marker_Trop_Human_genes', marker_genes_trop_g, ...
           'SDH', SDH_g, ...
           'RBH', RBH_g, ...
           'Local_enrichment', local_enrichment_g,   ...
           'Global_enrichment', global_enrichment_g, ...
           'Avg_expression', avg_expr_g, ...
           'Fraction_expression', fraction_expr_g, ...
           'Mutual_Info', MIG_g, ...
           'Mutual_NM_Info', MIG_NM_g, ...
           'Gini_Index', GI_g, ...
           'Cluster_Gini_Index', CGI_g, ...
           'Exp_Min', GeneExprGeneMin_g, ...
           'Exp_Max', GeneExprGeneMax_g, ...
           'Num_Count', NumCount_g);

%
%
% Convert to table
T = struct2table(s);

% Write table to csv file
writetable(T,'table_genes_specificity1.txt')

% Write table to tab file
writetable(T,'table_genes_specificity_tab1.txt', 'Delimiter', '\t');

clear T s Xnm_use GenesCells;

%}


%
%
%

TableOfspecificity  = readtable('table_genes_specificity_tab1.txt', 'ReadVariableNames', true);

states_g            = table2array(TableOfspecificity(1:end, 'State'));
marker_genes_lev_g  = table2array(TableOfspecificity(1:end, 'Marker_Laevis_IDs'));   
marker_genes_hum_g  = table2array(TableOfspecificity(1:end, 'Marker_Laevis_Human_genes'));   
marker_genes_trop_g = table2array(TableOfspecificity(1:end, 'Marker_Trop_Human_genes'));
SDH_g               = table2array(TableOfspecificity(1:end, 'SDH'));
RBH_g               = table2array(TableOfspecificity(1:end, 'RBH'));
local_enrichment_g  = table2array(TableOfspecificity(1:end, 'Local_enrichment'));
global_enrichment_g = table2array(TableOfspecificity(1:end, 'Global_enrichment'));
avg_expr_g          = table2array(TableOfspecificity(1:end, 'Avg_expression'));
fraction_expr_g     = table2array(TableOfspecificity(1:end, 'Fraction_expression'));
MIG_g               = table2array(TableOfspecificity(1:end, 'Mutual_Info'));
MIG_NM_g            = table2array(TableOfspecificity(1:end, 'Mutual_NM_Info'));
GI_g                = table2array(TableOfspecificity(1:end, 'Gini_Index'));
CGI_g               = table2array(TableOfspecificity(1:end, 'Cluster_Gini_Index'));
GeneExprGeneMin_g   = table2array(TableOfspecificity(1:end, 'Exp_Min'));
GeneExprGeneMax_g   = table2array(TableOfspecificity(1:end, 'Exp_Max'));
NumCount_g          = table2array(TableOfspecificity(1:end, 'Num_Count'));

%
%
%

disp('++++Main Prot Loop');

%
%
for m = 1 : size(ProteinStructure.ID,1)
  %  
  disp(m);
  %
  IndexGene = find(match(FGeneLevNames,ProteinStructure.ID(m)));  % find index of gene in mapping
  %
  TissueName  = "Unknown";
  TropName    = {'Unknown'};
  %
  if(~isempty(IndexGene))
    %
    TropName    = FGeneTropNames(IndexGene);  
    IndexVGene  = find(match(ExprTableNameFinal,FGeneTropNames(IndexGene))); % find index of gene V-score data 
    %
    if (~isempty(IndexVGene))
      SFV = V(IndexVGene);
      if (SFV > LowSpec)
        IndxMarker = find(match(marker_genes_lev_g,ProteinStructure.ID(m)));  
        if (~isempty(IndxMarker))
          ATissueName    = states_g(IndxMarker);
          TissueName     = ATissueName(1);
          ProteinStructure.VScore(m) = ProteinStructure.VScore(m) + global_enrichment_g(IndxMarker(1));
          for k = 2 : size(ATissueName,1)
            TissueName   = strcat(TissueName,';',{' '},ATissueName(k));
            ProteinStructure.VScore(m) = ProteinStructure.VScore(m) + global_enrichment_g(IndxMarker(k));       
          end
          ProteinStructure.VScore(m) = ProteinStructure.VScore(m) / size(ATissueName,1);
        else
          TissueName   = "Medium";   
        end  
      else
       TissueName   = "Low";   
      end
    end  
  end    
  %
  ProteinStructure.VTissueSpec(m) = TissueName;
  ProteinStructure.TropName{m}    = TropName;
  %
end

disp('++++End Main Prot Loop');


%
% Debugging Code (Please provide the gene of interest)
%

Names  = { 'Xelaev18016058m'; 'Xelaev18045905m'; 'Xelaev18043144m'};
Index  = find(match(ProteinStructure.ID, Names(1)));

Points = ProteinStructure.Data(Index,1:end);

% %
% % Debugging Section
% %
% 
% %TID        = 'Xelaev18018758m';
% %TSymbol    = 'BMP1';
% 
% TID         =  'Xelaev18001521m';
% TSymbol     =  'CTF1';
%     
% DebbugIndex = find(match(ProteinStructure.ID,TID));
% DebbugSymbl = ProteinStructure.Symbol(DebbugIndex);
% 
% if (~match(DebbugSymbl,TSymbol))
%   warning('Erorr');   
% end   
% 
% DebuggingData = ProteinStructure.Data(DebbugIndex,:); 
                                                       
% isAllOutputRequired = true; % Flag to output all proteins                                                       
% if (isAllOutputRequired)
%    % Data Filtering
%    Input_Title           = 'ALL_PROTEINS';
%    Input_Number_Cluster  = 117; 
%    l = 1;
%    for k=1:size(ProteinStructure.Symbol,1)
%       
%       name  = ProteinStructure.ID(k);
%                 
%       cname = char(name);
%       uns1  = find(match(cname(1:end),'##'));
%       uns2  = find(match(cname(1:end),'sp|'));
%       if (isempty(uns1) && isempty(uns2))
%         Input_List(l,1) = ProteinStructure.Symbol(k);
%         l = l+1;
%       end
%       
%    end
%    %
%    
%    AllSets = ['ALL_PROTEINS';AllSets];
%    
%    %
% end  


isAllOutputRequired = true; % Flag to output all proteins                                                       
if (isAllOutputRequired)   
   AllSets               = ['ALL_PROTEINS';AllSets];
   Cluster_Discription   = ['ALL PROTEINS';Cluster_Discription];
   Input_Number_Cluster  = 97;
   ALL_NamesID           = 0;   
   Cluster_NamesID       = [ALL_NamesID;  Cluster_NamesID];   
   Number_Cluster        = [Input_Number_Cluster; Number_Cluster];
end  

%
%   

isNormalizedData      = 0;
isClusterCenetersUsed = 0;

% 
% % '0' is the original data
% % '1' is the row normalized
% % '2' is the median normalized 
% 
% Input_List        = unique(Input_List); % safe guard for unique symbols only
% [MeanDynamis, ok] = ClustByGeneSet(Input_List,            ... 
%                                    Input_Title,           ...
%                                    Input_Number_Cluster,  ... 
%                                    Input_PROThours,       ... 
%                                    ProteinStructure,      ...
%                                    N,                     ...
%                                    Map,                   ...
%                                    isNormalizedData);
%                                               
% %
% % Run multiple realizations
% %
% 
% close all      

%
%
%

FileFolder    =  mfilename('fullpath');
CFileFolder   =  strrep(FileFolder,'ClusteringRunAllRed','');
Directory     =  char(strcat(CFileFolder,'ClustByGeneSet'));  
mkdir(Directory);     

ext           = '.html';
file          = char(strcat(Directory,'/','Index_All_Zip')); 
case_file     = strcat(file, ext);

fid_zip       = fopen(case_file, 'w');
OK            = 1;
%
if (fid_zip < 0)
   OK = 0;
   exit;
end
%
%
% FOR DEBUGGING
%

%AllSets         = AllSets(1);

%
Dynamicity      = zeros(size(AllSets,1),1);
%
GHdist          = [];
GinSet          = zeros(size(AllSets,1),1);
Gdetected       = zeros(size(AllSets,1),1); 
GinGenome       = zeros(size(AllSets,1),1); 
GuniqInGenome   = zeros(size(AllSets,1),1); 
GuniqDetected   = zeros(size(AllSets,1),1);
GunDetected     = zeros(size(AllSets,1),1);
GMedianConcentr = zeros(size(AllSets,1),1);
%
%
%//TableNameSchwanhouser   = readtable('data/Schwanhouser_Mol.txt');
%//SchwanhouserName        = table2array(TableNameSchwanhouser(1:end, 'Var1'));
%//SchwanhouserMoles       = table2array(TableNameSchwanhouser(1:end, 'Var2'));
%

disp('++++Main Clustering');
Ycent_ALL             = [];
%
for i = 1:size(AllSets,1)
   %
   if (~isempty(find(match(AllSets(i),'ALL_PROTEINS'))))
      %
      Input_Title           = char(AllSets(i));
      Input_Number_Cluster  = 117; 
      l                     = 1;
      Ycent_ALL             = [];
      %
      for k=1:size(ProteinStructure.Symbol,1)

        name  = ProteinStructure.ID(k);

        cname = char(name);
        uns1  = find(match(cname(1:end),'##'));
        uns2  = find(match(cname(1:end),'sp|'));
        if (isempty(uns1) && isempty(uns2))
          Input_List(l,1) = ProteinStructure.Symbol(k);
          l = l+1;
        end
      
      end
      %
   else    
     
     Input_Title    = char(AllSets(i));    
     Input_Title_F  = fullfile('geneset', strcat(Input_Title,'.txt')); 

     T_Table        = readtable(Input_Title_F, 'ReadVariableNames', false);
     Input_List     = upper(table2array(T_Table(1:end, 'Var1')));
     if (Cluster_NamesID(i) == 1)
       Input_List     = table2array(T_Table(1:end, 'Var1'));
     end    
   end    
   %    
   disp(Input_Title);
   %
   
   Input_List       = unique(Input_List); % safe guard for unique symbols only
  
   %//if (match(Input_Title, {'RNApolymerase'}))
   %// disp('end');
   %//end    


   %
   % Let's print the discretion 
   %
   
   Input_Title_File = Input_Title;
   Input_Title      = char(Cluster_Discription(i));
   Input_NamesID    = Cluster_NamesID(i);
   %
   [MedianDynamis, inSet, detected, inGenome,           ...
    uniqInGenome, uniqDetected, unDetected,             ...
    MedianConcentr, MedianVscore,                       ...
    Hdist,Ycent, DMean, DVarn, GroupData,               ...
    GroupDataEgg, GroupDataLast,                        ... 
    DMeanEgg, DVarnEgg, DMeanLast, DVarnLast,           ...
    GroupMedian, ok]                                    ...
              = ClustByGeneSet(Input_List,              ... 
                               Input_Title_File,        ...
                               Input_Title,             ...
                               Input_NamesID,           ...
                               Number_Cluster(i),       ... 
                               Input_PROThours(1:N),    ... 
                               ProteinStructure,        ...
                               N,                       ...
                               Map,                     ...
                               isNormalizedData,        ...
                               isClusterCenetersUsed, Ycent_ALL);
   %
   if (~isempty(find(match(AllSets(i),'ALL_PROTEINS'))))
     Ycent_ALL = Ycent;
   end                        
   %   
   GHdist              = [GHdist;Hdist];
   Dynamicity(i)       = MedianDynamis;
   GinSet(i)           = inSet; 
   Gdetected(i)        = detected; 
   GinGenome(i)        = inGenome; 
   GuniqInGenome(i)    = uniqInGenome; 
   GuniqDetected(i)    = uniqDetected;
   GunDetected(i)      = size(unDetected,1);
   GMedianConcentr(i)  = MedianConcentr;
   GMedianVscore(i)    = MedianVscore;
   GDMean(i)           = DMean;
   GDVarn(i)           = DVarn;
   TableRecord(i)      = ok;
   GGroupData(i,:)     = GroupData;
   GGroupDataEgg(i,:)  = GroupDataEgg;
   GGroupDataLast(i,:) = GroupDataLast;
   GGroupMedian(i,:)   = GroupMedian;
   GDMeanEgg(i)        = DMeanEgg; 
   GDVarnEgg(i)        = DVarnEgg; 
   GDMeanLast(i)       = DMeanLast; 
   GDVarnLast(i)       = DVarnLast;
   
   %
   % Check 
   %
   
   %TableNameSchwanhouser   = readtable('data/Schwanhouser_Mol.txt');
   %SchwanhouserName        = table2array(TableNameSchwanhouser(1:end, 'Var1'));
   %SchwanhouserMoles       = table2array(TableNameSchwanhouser(1:end, 'Var2'));
   %
   %[SchwanhouserGroup, IND_A, IND_B] = intersect(SchwanhouserName,Input_List);
   %
   %GSchwanMean       = mean(log10(SchwanhouserMoles(IND_A)+1.e-11));
   %GSchwanVarn       = std(log10(SchwanhouserMoles(IND_A)+1.e-11));
   %
   %[m,s]             = normfit(log10(SumCRes+1.e-11));
   %X                 = 0:0.01:max(log(SchwanhouserMoles+1.e-11));
   %GSchwanGroupData  = normpdf(X,m,s);

   %  
   %
   ext            = '.zip';
   file           = char(strcat('http://xenopus.hms.harvard.edu/Clustring/ClustByGeneSet','/',Input_Title,'/',Input_Title)); 
   case_folder    = strcat(file, ext);
   file_zip       = strcat(Input_Title, ext);
   %
   html_string = strcat('<a href="',case_folder,'"','>',file_zip,'</a><br>');
   fprintf(fid_zip, '%s\n',char(html_string)); 
   %
end


%
close all
fclose('all');

[SSGDMean, SSIndx] = sort(GDMean,'descend');
disp([char(AllSets(SSIndx)),strcat(num2str(round(transpose(GDMean(SSIndx)),3)),'....'),num2str(round(transpose(GDVarn(SSIndx)),3))]); 

for i=1:size(AllSets,1)
   [m,s]         = normfit(GGroupData(i,:));
   MeanFitted(i) = m; 
   StdFitted(i)  = s; 
end

[TSSGDMean, TSSIndx] = sort(MeanFitted,'descend');
disp([char(AllSets(SSIndx)),strcat(num2str(round(transpose(MeanFitted(TSSIndx)),3)),'....'),num2str(round(transpose(StdFitted(TSSIndx)),3))]); 

%
%
%

TimeTrend         = [-2; ... %// oo56 
                      0; ... %// egg   Yes (included)
                      9; ... %// st09  Yes (included)
                     12; ... %// st12  Yes (included)
                     17; ... %// st17
                     22; ... %// st22  Yes (included)
                     24; ... %// st24
                     26; ... %// st26
                     30; ... %// st30  Yes (included)
                     42; ... %// st42
                    ];
                    
for i=1:size(GGroupMedian,1)
 %
 %
 LineWidth = 2;
 hold on;
 %plot(Input_PROThours(1:10),GGroupMedian(i,:),'LineWidth',LineWidth); 
 plot(TimeTrend(1:N),GGroupMedian(i,:),'LineWidth',LineWidth); 
 %leg = AllSets(i);
 leg  = strcat(AllSets(i),{' '},'(',num2str(GuniqDetected(i)),'/', ...
                                    num2str(GuniqInGenome(i)),')');         
 %
 legendInfo{i} = char(leg);
 % 
end
%
%
%
legend1 = legend(legendInfo,'Location','northeastoutside'); 
set(gca,'XTick',[-2 0 9 12 17 22 24 26 30 42],'XTickLabel', {'','2','9','12','17','','24','26','30','42'});                    
% Create ylabel
ylabel('Relative Abundance');
% Create xlabel
xlabel('Time, hpf');
% Set the remaining axes properties
set(gca,'FontSize',12);
%
title(char('Protein Relative abundance by genesets (median)'), 'FontWeight', 'normal');  

set(gca,'FontSize',14);
    
%
% Sets of out interests
%
%
% Define some colors and plotting options 
plt.dump         = 0;
plt.Gray         = [.4 .4 .4]; 
plt.GrayLight    = [.6 .6 .6]; 
plt.Cyan         = [.0 1. 1.]; 
plt.Blue         = [.04 .14 .98];   % MATHEMATICA  blue 
plt.Green        = [.16 1. .18];    % MATHEMATICA green
plt.GreenDark    = [.26 .58 0.17];  % Dark Green
plt.Magenta      = [1. 0. 1.];      % Magenta 
plt.Khaki        = [.52 .38 .12];   % Khaki 
%
%
rnb.Gray         = [192 192 192]/255; 
rnb.Blue         = [0 0 255]/255;
rnb.Orange       = [255 165 0]/255;
rnb.BlueLight    = [65 105 225]/255;
rnb.Green        = [0 128 0]/255;
rnb.Magenta      = [139 0 139]/255;
rnb.OrangeRed    = [255 69 0]/255; 
%
%
Temp = [90;  ...
        89;  ...
        88;  ... 
        87;  ...
        86;  ...
        85;  ...
        84;  ...
        83;  ...
        82;  ...
        81;  ...
        80;  ...
        79;  ...
        78;  ...
        77;  ...
        76;  ...
        75;  ...
        74;  ...
        73;  ...
        72;  ...
        71;  ... 
        70;  ...
        69;  ...
        68;  ...
        67;  ...
        66;  ...
        65;  ...
        64;  ... 
        63;  ...
        62;  ...
        61;  ...
        60;  ...
        59;  ...
        58;  ...
        57;  ...
        56;  ...
        55;  ...
        54;  ... 
        53;  ...
        52;  ...
        51;  ...
        50;  ...
        49;  ...
        48;  ...
        47;  ...
        46;  ...
        45;  ...
        44;  ...
        43;  ...
        42;  ...
        41;  ...
        40;  ...
        39;  ...
        38;  ...
        37;  ...
        36;  ...
        35;  ...
        34;  ...
        33;  ...
        32;  ...
        31;  ...
        30;  ...
        29;  ...
        28;  ...
        27;  ...
        26;  ...
        25;  ...
        24;  ...
        23;  ...
        22;  ...
        21;  ...
        20;  ...
        19;  ...
        18;  ... 
        17;  ...
        16;  ...
        15;  ...
        14;  ...
        13;  ...
        12;  ...
        11;  ...
        10;  ...
        9;   ...
        8;   ...
        7;   ...
        6;   ...
        5;   ...
        4;   ...
        3;   ...
        2;   ...
        1]; 
%
%
% %
% OArrayColors = [ 255 14 240;	 ...
%                 255 13 240;	 ...
%                 255 12 240;	 ...
%                 255 11 240;	 ...
%                 255 10 240;	 ...
%                 255 9  240;	 ...
%                 255 8  240;	 ...
%                 255 7  240;	 ...
%                 255 6  240;	 ...
%                 255 5 240;	 ...
%                 255 4 240;	 ...
%                 255 3 240;	 ...
%                 255 2 240;	 ...
%                 255 1 240;	 ...
%                 255 0 240;   ...
%                 255 0 224;	 ...
%                 255 0 208;	 ...
%                 255 0 192;	 ...
%                 255 0 176;	 ...
%                 255 0 160;	 ...
%                 255 0 144;	 ...
%                 255 0 128;	 ...
%                 255 0 112;   ...
%                 255 0 96;	 ...
%                 255 0 80;	 ...
%                 255 0 64;	 ...
%                 255 0 48;	 ...
%                 255 0 32;	 ...
%                 255 0 16;	 ...
%                 255 0 0;     ...
%                 255 10 0;	 ...
%                 255 20 0;	 ...
%                 255 30 0;	 ...
%                 255 40 0;	 ...
%                 255 50 0;	 ...
%                 255 60 0;	 ...
%                 255 70 0;	 ...
%                 255 80 0;	 ...
%                 255 90 0;	 ...
%                 255 100 0;	 ...
%                 255 110 0;	 ...
%                 255 120 0;	 ...
%                 255 130 0;	 ...
%                 255 140 0;	 ...
%                 255 150 0;	 ...
%                 255 160 0;	 ...
%                 255 170 0;	 ...
%                 255 180 0;	 ...
%                 255 190 0;	 ...
%                 255 200 0;	 ...
%                 255 210 0;	 ...
%                 255 220 0;	 ...
%                 255 230 0;	 ...
%                 255 240 0;	 ...
%                 255 250 0;	 ...
%                 253 255 0;	 ...
%                 215 255 0;	 ...
%                 176 255 0;	 ...
%                 138 255 0;	 ...
%                 101 255 0;	 ...
%                 62 255 0;	 ...
%                 23 255 0;	 ...
%                 0 255 16;	 ...
%                 0 255 54;	 ...
%                 0 255 92;	 ...
%                 0 255 13;    ...
%                 0 255 168;	 ...
%                 0 255 208;	 ...
%                 0 255 244;	 ...
%                 0 228 255;	 ...
%                 0 212 255;	 ...
%                 0 196 255;	 ...
%                 0 180 255;	 ...
%                 0 164 255;   ...
%                 0 148 255;	 ...
%                 0 132 255;	 ...
%                 0 116 255;	 ...
%                 0 100 255;	 ...
%                 0 84 255;    ...
%                 0 68 255;    ...
%                 0 50 255;	 ...
%                 0 34 255;    ...
%                 0 18 255;	 ...
%                 0 2 255;	 ...
%                 0 0 255;	 ...
%                 1 0 255;     ...
%                 2 0 255;     ...
%                 3 0 255;     ...
%                 4 0 255;	 ...
%                 5 0 255]/255;
%             

OArrayColors = [85  0   0;   ... 
                175 0   0;   ...      
                255 0   176; ...
                255 0   16;	 ...
                255 150 0;	 ...
                255 220 0;	 ...
                215 255 0;	 ...
                176 255 0;	 ...   
                138 255 0;	 ...
                253 255 0;   ... % netrual point
                176 255 0;   ...
                0   128 0;   ...
                0   255 168; ... 
                138 255 0;   ...
                23  255 0;   ...    
                0   255 131; ... 
                0   255 244; ... 
                0   212 255; ...
                0   132 255; ...
                0   0   255]/255; 
                      
 %

 LOOK_SET = {'ALL_PROTEINS';      ...
             'Glycolysis';        ...
             'TFs';               ...
             'Ribosome';          ...
             'Proteosome';        ...
             'Spliceosome';       ...
             'TCA';               ... 
             'DNArepli';          ...
             'E3s';               ... 
             'Kinome';            ...
             'SecretedFactors';   ...
             'MitoCartaMouseV2';  ...
             'KEGG_ENDOCYTOSIS'};
        
         
%     {'Glycolysis'                    }
%     {'TFs'                           }
%     {'Lysosome'                      }
%     {'PentosePhosphPath'             }
%     {'Proteosome'                    }
%     {'Ribosome'                      }
%     {'Spliceosome'                   }
%     {'TCA'                           }
%     {'CellCycle'                     }
%     {'DNArepli'                      }
%     {'RNAdegrad'                     }
%     {'RNApolymerase'                 }
%     {'E3s'                           }
%     {'Receptors'                     }
%     {'Ligands'                       }
%     {'FattyAcidMet'                  }
%     {'Kinome'                        }
%     {'SecretedFactors'               }
%     {'ECM'                           }
%     {'MitoCartaMouseV2'              }
%     {'Collagens'                     }
%     {'OxidativePhosphorylation'      }
%     {'GO_PHOSPHATASE'                }
%     {'GO_CHROMATIN_REMODELING'       }
%     {'HALLMARK_MITOTIC_SPINDLE'      }
%     {'HALLMARK_XENOBIOTIC_METABOLISM'}
%     {'KEGG_ENDOCYTOSIS'              }
%     {'KEGG_PEROXISOME'               }
%     {'CPC'                           }         
         
        
%LOOK_SET = AllSets;        
 
Ind1             = find(match(AllSets,LOOK_SET));
%
Vnames           = AllSets(Ind1);
Vmean            = GDMean(Ind1);
[SVmean, SInd]   = sort(Vmean,'descend');
SVnames          = Vnames(SInd);

InD_ALL          = find(match(SVnames,'ALL_PROTEINS')); 


DHColor          = round(0.5*size(OArrayColors,1)/size(Vnames(1:InD_ALL-1),1),-1);
DCColor          = round(0.5*size(OArrayColors,1)/size(Vnames(InD_ALL+1:end),1),-1);


if (DHColor > 1) 
  ArrayColorsH   = OArrayColors(1:DHColor:min(9,DHColor*size(Vnames(1:InD_ALL-1))),:);
else
  ArrayColorsH   = OArrayColors(1:1:size(Vnames(1:InD_ALL-1,1)),:);
end

if (DCColor > 1)
  ArrayColorsC   = OArrayColors(11:DCColor:min(size(OArrayColors,1),11+DCColor*size(Vnames(InD_ALL+1:end),1)),:);
else
  ArrayColorsC   = OArrayColors(end-size(Vnames(InD_ALL+1:end),1):1:end,:);
end    

ArrayColors      = [ArrayColorsH; ArrayColorsC];

%
figureConDist    = figure; 
legendInfo       = {};
gl               = 0;
dx               = (max(Ycent_ALL)-min(Ycent_ALL))/1000;
X                = min(Ycent_ALL):dx:max(Ycent_ALL);
%
%
%
for i=1:size(SVnames,1)
    
  %
  %
  %
  
  if (~isempty(find(match(AllSets,SVnames(i)))))
    %
    %
    Ind   = find(match(AllSets,SVnames(i)));  
    X     = -3:0.01:5;  
    Y     = GGroupData(Ind,:); 
    %
    %Y(1:size(X,2)) = interp1(Ycent_ALL,GHdist(Ind,:),X,'cubic','extrap'); % extrapolate
    %Y(Y<0)         = 0.0;
    % 
    gColor          = ArrayColors(i,:);
    %
    %
    %
    LineWidth       = 2;
    leg             = strcat(Cluster_Discription(Ind),{' '},'(',num2str(GuniqDetected(Ind)),'/', ...
                                                                num2str(GuniqInGenome(Ind)),')');         
    %
    if (match(SVnames(i),'ALL_PROTEINS')) 
      gColor         = plt.Gray;
      LineWidth      = 6;
      leg            = strcat(Cluster_Discription(Ind),{' '},num2str(GuniqDetected(Ind)));         
    end
    %
    legendInfo{i} = char(leg);
    %
    plot(X, Y, 'LineWidth',LineWidth, 'Color', gColor);
    hold on;
    % queuees
    %
  end
  %    
end

legend1 = legend(legendInfo,'Location','northeastoutside'); 

% Create ylabel
ylabel('Probability density function');

% Create xlabel
xlabel('Protein concentration (nM, log_{10} units)');
    
% Set the remaining axes properties
set(gca,'FontSize',12);

%
title(char('Protein concentration by genesets'), 'FontWeight', 'normal');  

set(gca,'FontSize',14);


%
close all
fclose('all');

%
% Make Kolmogorov-Smirnov Test to check if the come from the same
% distribution
%

for k  = 1 : size(GGroupDataEgg,1) 
 %
 Y1    = GGroupDataEgg (k,:);
 Y2    = GGroupDataLast(k,:); 
 %
 [h,p] = kstest2(Y1,Y2,'Alpha',0.05);
 StrKolSmirnov =  char(strcat('For group: ', char(AllSets(k)), ' Kolmogorov-Smirnov Test',', ',' h = ',{' '}, num2str(h),', ', ...
                              ' p = ',{' '}, num2str(p),', ',' \alpha = ',{' '}, num2str(0.05),'  ')); 
 %
 disp(StrKolSmirnov);
 %
end 

%
% This is the place to plot distribution in stage 2
%

%//[SSGDMean, SSIndx] = sort(GDMeanEgg,'descend');
%//
%//disp('This is table for egg');
%//disp([char(AllSets(SSIndx)),strcat(num2str(round(transpose(GDMeanEgg(SSIndx)),3)),'....'),num2str(round(transpose(GDVarnEgg(SSIndx)),3))]); 

for i=1:size(AllSets,1)
  [m,s]            = normfit(GGroupDataEgg(i,:)); 
  MeanFittedEgg(i) = m; 
  StdFittedEgg(i)  = s; 
end


[TSSGDMean, TSSIndx] = sort(MeanFittedEgg,'descend');

disp('This is table for egg');
disp([char(AllSets(TSSIndx)),strcat(num2str(round(transpose(MeanFittedEgg(TSSIndx)),3)),'....'),num2str(round(transpose(StdFittedEgg(TSSIndx)),3))]); 

%
%
%

Ind1             = find(match(AllSets,LOOK_SET));
%
Vnames           = AllSets(Ind1);
Vmean            = GDMeanEgg(Ind1);
[SVmean, SInd]   = sort(Vmean,'descend');
SVnames          = Vnames(SInd);

InD_ALL          = find(match(SVnames,'ALL_PROTEINS')); 


DHColor          = round(0.5*size(OArrayColors,1)/size(Vnames(1:InD_ALL-1),1),-1);
DCColor          = round(0.5*size(OArrayColors,1)/size(Vnames(InD_ALL+1:end),1),-1);


if (DHColor > 1) 
  ArrayColorsH   = OArrayColors(1:DHColor:min(9,DHColor*size(Vnames(1:InD_ALL-1))),:);
else
  ArrayColorsH   = OArrayColors(1:1:size(Vnames(1:InD_ALL-1,1)),:);
end

if (DCColor > 1)
  ArrayColorsC   = OArrayColors(11:DCColor:min(size(OArrayColors,1),11+DCColor*size(Vnames(InD_ALL+1:end),1)),:);
else
  ArrayColorsC   = OArrayColors(end-size(Vnames(InD_ALL+1:end),1):1:end,:);
end    

ArrayColors      = [ArrayColorsH; ArrayColorsC];

%
figureConDist    = figure; 
legendInfo       = {};
gl               = 0;
dx               = (max(Ycent_ALL)-min(Ycent_ALL))/1000;
X                = min(Ycent_ALL):dx:max(Ycent_ALL);
%
%
%
for i=1:size(SVnames,1)
    
  %
  %
  %
  
  if (~isempty(find(match(AllSets,SVnames(i)))))
    %
    %
    Ind   = find(match(AllSets,SVnames(i)));  
    X     = -3:0.01:5;  
    Y     = GGroupDataEgg(Ind,:); 
    %
    %Y(1:size(X,2)) = interp1(Ycent_ALL,GHdist(Ind,:),X,'cubic','extrap'); % extrapolate
    %Y(Y<0)         = 0.0;
    % 
    gColor          = ArrayColors(i,:);
    %
    %
    %
    LineWidth       = 2;
    leg             = strcat(Cluster_Discription(Ind),{' '},'(',num2str(GuniqDetected(Ind)),'/', ...
                                                                num2str(GuniqInGenome(Ind)),')');         
    %
    if (match(SVnames(i),'ALL_PROTEINS')) 
      gColor         = plt.Gray;
      LineWidth      = 6;
      leg            = strcat(Cluster_Discription(Ind),{' '},num2str(GuniqDetected(Ind)));         
    end
    %
    legendInfo{i} = char(leg);
    %
    plot(X, Y, 'LineWidth',LineWidth, 'Color', gColor);
    hold on;
    % queuees
    %
  end
  %    
end

legend1 = legend(legendInfo,'Location','northeastoutside'); 

% Create ylabel
ylabel('Probability density function');

% Create xlabel
xlabel('Protein concentration (nM, log_{10} units)');
    
% Set the remaining axes properties
set(gca,'FontSize',12);

%
title(char('Protein concentration in stage 2 by genesets'), 'FontWeight', 'normal');  

set(gca,'FontSize',14);

%
% This is the place to plot distribution in stage 42
%

%//
%//[SSGDMean, SSIndx] = sort(GDMeanLast,'descend');
%//
%//disp('This is table for last');
%//disp([char(AllSets(SSIndx)),strcat(num2str(round(transpose(GDMeanLast(SSIndx)),3)),'....'),num2str(round(transpose(GDVarnLast(SSIndx)),3))]); 
%//

%
%
%

for i=1:size(AllSets,1)
  [m,s]             = normfit(GGroupDataLast(i,:)); 
  MeanFittedLast(i) = m; 
  StdFittedLast(i)  = s; 
end

[TSSGDMean, TSSIndx] = sort(MeanFittedLast,'descend');

disp('This is table for egg');
disp([char(AllSets(TSSIndx)),strcat(num2str(round(transpose(MeanFittedLast(TSSIndx)),3)),'....'),num2str(round(transpose(StdFittedLast(TSSIndx)),3))]); 

%
%
%

Ind1             = find(match(AllSets,LOOK_SET));
%
Vnames           = AllSets(Ind1);
Vmean            = GDMeanLast(Ind1);
[SVmean, SInd]   = sort(Vmean,'descend');
SVnames          = Vnames(SInd);

InD_ALL          = find(match(SVnames,'ALL_PROTEINS')); 


DHColor          = round(0.5*size(OArrayColors,1)/size(Vnames(1:InD_ALL-1),1),-1);
DCColor          = round(0.5*size(OArrayColors,1)/size(Vnames(InD_ALL+1:end),1),-1);

if (DHColor > 1) 
  ArrayColorsH   = OArrayColors(1:DHColor:min(9,DHColor*size(Vnames(1:InD_ALL-1))),:);
else
  ArrayColorsH   = OArrayColors(1:1:size(Vnames(1:InD_ALL-1,1)),:);
end

if (DCColor > 1)
  ArrayColorsC   = OArrayColors(11:DCColor:min(size(OArrayColors,1),11+DCColor*size(Vnames(InD_ALL+1:end),1)),:);
else
  ArrayColorsC   = OArrayColors(end-size(Vnames(InD_ALL+1:end),1):1:end,:);
end    

ArrayColors      = [ArrayColorsH; ArrayColorsC];

%
figureConDist    = figure; 
legendInfo       = {};
gl               = 0;
dx               = (max(Ycent_ALL)-min(Ycent_ALL))/1000;
X                = min(Ycent_ALL):dx:max(Ycent_ALL);
%
%
%
for i=1:size(SVnames,1)
    
  %
  %
  %
  
  if (~isempty(find(match(AllSets,SVnames(i)))))
    %
    %
    Ind   = find(match(AllSets,SVnames(i)));  
    X     = -3:0.01:5;  
    Y     = GGroupDataLast(Ind,:); 
    %
    %Y(1:size(X,2)) = interp1(Ycent_ALL,GHdist(Ind,:),X,'cubic','extrap'); % extrapolate
    %Y(Y<0)         = 0.0;
    % 
    gColor          = ArrayColors(i,:);
    %
    %
    %
    LineWidth       = 2;
    leg             = strcat(Cluster_Discription(Ind),{' '},'(',num2str(GuniqDetected(Ind)),'/', ...
                                                                num2str(GuniqInGenome(Ind)),')');         
    %
    if (match(SVnames(i),'ALL_PROTEINS')) 
      gColor         = plt.Gray;
      LineWidth      = 6;
      leg            = strcat(Cluster_Discription(Ind),{' '},num2str(GuniqDetected(Ind)));         
    end
    %
    legendInfo{i} = char(leg);
    %
    plot(X, Y, 'LineWidth',LineWidth, 'Color', gColor);
    hold on;
    % queuees
    %
  end
  %    
end

legend1 = legend(legendInfo,'Location','northeastoutside'); 

% Create ylabel
ylabel('Probability density function');

% Create xlabel
xlabel('Protein concentration (nM, log_{10} units)');
    
% Set the remaining axes properties
set(gca,'FontSize',12);

%
title(char('Protein concentration in stage 42 by genesets'), 'FontWeight', 'normal');  

set(gca,'FontSize',14);




%
close all
fclose('all');


%[SAllSets, ISort]   = sort(AllSets);
%SDynamicity         = Dynamicity(ISort);

%also can you sort by proportion of gene set which is detected (last column divided by 
%the cardinality of gene set) as we discussed, for Phospho data, no tby dynamicity 

SortThings           = GuniqDetected./GinSet;

%[SDynamicity, ISort]= sort(Dynamicity, 'descend');

[Res, ISort]         = sort(SortThings, 'descend');

SDynamicity          = round(Dynamicity(ISort),5);
SAllSets             = AllSets(ISort);
SCluster_Discription = Cluster_Discription(ISort);
SGinSet              = GinSet(ISort); 
SGdetected           = Gdetected(ISort); 
SGinGenome           = GinGenome(ISort); 
SGuniqInGenome       = GuniqInGenome(ISort); 
SGuniqDetected       = GuniqDetected(ISort);
SGunDetected         = GunDetected(ISort);
SMedianConcentr      = round(GMedianConcentr(ISort),1);  
SMedianVscore        = round(GMedianVscore(ISort),2);
STableRecord         = TableRecord(ISort);
SGDMean              = GDMean(ISort);
SGDVarn              = GDVarn(ISort);

Title                = Input_Title;
CFileFolder          = mfilename('fullpath');
Directory            = strrep(CFileFolder,'ClusteringRunAllRed','ClustByGeneSet');

ext                  = '.txt';
%//file              = char(strcat(Directory,'/',Title,'_Clusters')); 
file                 = char(strcat(Directory,'/','Cluster_Dynamicity')); 
case_file            = strcat(file, ext);

fid                  = fopen(case_file, 'w');
OK                   = 1;

%
if (fid < 0)
   OK = 0;
   exit;
end
%

for k = 1 : size(SAllSets,1)
  fprintf(fid, '%s\t\t%3.5f\n', char(SAllSets(k)),SDynamicity(k)); 
end   

%
% Print the Global HTML header
%

ext               = '.html';
%//file           = char(strcat(Directory,'/',Title,'_Clusters')); 
file              = char(strcat(Directory,'/','Clustering_Main_WebPage')); 
case_file         = strcat(file, ext);

fid_Main_WebPage  = fopen(case_file, 'w');
OK                = 1;

%
if (fid_Main_WebPage < 0)
   OK = 0;
   exit;
end
%

fprintf(fid_Main_WebPage, '%s\n', char('<head><title>Xenopus at Systems Biology, HMS</title>')); 
fprintf(fid_Main_WebPage, '%s\n', char('<head><title>Xenopus at Systems Biology, HMS</title>'));
fprintf(fid_Main_WebPage, '%s\n', char('<meta name="robots" content="noindex">'));
fprintf(fid_Main_WebPage, '%s\n', char('<meta name="googlebot" content="noindex">'));
fprintf(fid_Main_WebPage, '%s\n', char('</head>'));
fprintf(fid_Main_WebPage, '%s\n', char('<body>'));

fprintf(fid_Main_WebPage, '%s\n', char('<center><pre><font face="Garamond">'));
fprintf(fid_Main_WebPage, '%s\n', char('<img halign=middle border=0 src="http://xenopus.hms.harvard.edu/Oocyte/ClusteringLogo.png" height=400><br>'));
fprintf(fid_Main_WebPage, '%s\n', char('This webpage provides supplementary materials for Clustering and browser for <b><i>"To Be Decided"</i></b>, <a href="http://www.cell.com">Peshkin et al., To Be Decided, 2019</a>.'));


%fprintf(fid_Main_WebPage, '%s\n', char('</pre></center>'));      


fprintf(fid_Main_WebPage, '%s\n', char('<head>'));
fprintf(fid_Main_WebPage, '%s\n', char('<style>'));
fprintf(fid_Main_WebPage, '%s\n', char('table {'));
fprintf(fid_Main_WebPage, '%s\n', char('font-family: arial, sans-serif;'));
fprintf(fid_Main_WebPage, '%s\n', char('border-collapse: collapse;'));
fprintf(fid_Main_WebPage, '%s\n', char('width: 80%;'));
fprintf(fid_Main_WebPage, '%s\n', char('}'));

fprintf(fid_Main_WebPage, '\n');
fprintf(fid_Main_WebPage, '\n');

fprintf(fid_Main_WebPage, '%s\n', char('td, th {'));
fprintf(fid_Main_WebPage, '%s\n', char('border: 1px solid #dddddd;'));
fprintf(fid_Main_WebPage, '%s\n', char('text-align: left;'));
fprintf(fid_Main_WebPage, '%s\n', char('padding: 8px;'));
fprintf(fid_Main_WebPage, '%s\n', char('}'));

fprintf(fid_Main_WebPage, '%s\n', char('tr:nth-child(even) {'));
fprintf(fid_Main_WebPage, '%s\n', char('background-color: #dddddd;'));
fprintf(fid_Main_WebPage, '%s\n', char('}'));
fprintf(fid_Main_WebPage, '%s\n', char('</style>'));
fprintf(fid_Main_WebPage, '%s\n', char('</head>'));
fprintf(fid_Main_WebPage, '%s\n', char('<body>'));

fprintf(fid_Main_WebPage, '\n');
fprintf(fid_Main_WebPage, '\n');


fprintf(fid_Main_WebPage, '%s\n', char('<h2>Table of gene sets features</h2>'));
fprintf(fid_Main_WebPage, '%s\n', char('<table>'));

fprintf(fid_Main_WebPage, '\n');
fprintf(fid_Main_WebPage, '\n');

fprintf(fid_Main_WebPage, '%s\n', char('<tr>'));
%
%

fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('GeneSet/Clusters'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="footnote1" title="Gene Set identifier, which links to per-set page, presenting the dynamics and gene membership, as well as the link to respective gene function as detailed in GeneCards entry."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('Median Dynamicity'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="#footnote-1" title="The median dynamicity (i.e. deviation from a flat line) of temporal patterns in respective set."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('In Set'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="#footnote-1" title="A number of gene symbols in the set definition (based on human genes) as a link to complete list of gene symbols."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('Total In Genome'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="#footnote-1" title="The total number of genes in the genome which correspond to the gene symbols in a given set. N.B. possibly several genomic loci(genes) map to the same human gene symbol, resulting in a higher number here than in &quot;In Set&quot; column."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('Total Detected'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="#footnote-1" title="The total number of objects (proteins or phospho peptideds) detected in the experimental measurement, which map into a given gene set by symbol."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('Total unDetected'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="#footnote-1" title="The total number of symbols in the gene set not identified in the experiment because these were ether not found in the genome or found but not detected."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('Unique in Genome'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="#footnote-1" title="The number of unique symbols which are found in the genome. This is theoretical maximum number of &quot;Unique Detected&quot; if we detected all genes found in the genome."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('Unique Detected'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="#footnote-1" title="The number of unique gene symbols from the set detected in the experiment."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('Median concentration (nM)'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="#footnote-1" title="Median concentration of unique gene symbols from the set detected in the experiment."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('<th>'));
fprintf(fid_Main_WebPage, '%s\n', char('Tissue Specificity'));
fprintf(fid_Main_WebPage, '%s\n', char('<a href="#footnote-1" title="Tissue specificity as discribed by the median V-score of unique gene symbols from the set detected in the experiment."><sup><img src="http://xenopus.hms.harvard.edu/ClustByGeneSet_Dev_Phos/HelpImage.png" width="8%" hight="8%"/></sup></a>'));
fprintf(fid_Main_WebPage, '%s\n', char('</th>'));
%
fprintf(fid_Main_WebPage, '%s\n', char('</tr>'));

fprintf(fid_Main_WebPage, '\n');
fprintf(fid_Main_WebPage, '\n');

%
%

for k = 1 : size(SAllSets,1)
  %
  %
  
  if (STableRecord(k) == true)
    %
    %
    fprintf(fid_Main_WebPage, '%s\n', char('<tr>'));
    %
    % Clusters
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',                    ...
    strcat('<a href="http://xenopus.hms.harvard.edu/Clustering/ClustByGeneSet/', ...
           char(SAllSets(k)),'/',char(SAllSets(k)),'.html">',                ...
           '<strong style="font-size: 12px;">',                              ...  
           SCluster_Discription(k),                                          ...
           '</strong><br />',                                                ...         
           '</a>'),                                                          ...
           '</td>')));
            
    %
    % Median Dynamicity
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SDynamicity(k)),'</a>','</td>')));

    %//fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>', ...
    %//strcat('<a href="http://xenopus.hms.harvard.edu/Clustering/ClustByGeneSet/', ...
    %//       char(SAllSets(k)),'/',char(SAllSets(k)),'.html">', num2str(SDynamicity(k)),'</a>'), ...
    %//       '</td>')));
     
    %
    % In Set
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>', ...
    strcat('<a href="http://xenopus.hms.harvard.edu/Clustering/ClustByGeneSet/',char(SAllSets(k)),'/',char(SAllSets(k)),'_Set.txt">', ... 
           num2str(SGinSet(k)),'</a>'), '</td>')));
    %
    % Total in Genome
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SGinGenome(k)),'</a>','</td>')));
    %
    % Total Detected
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>', ...
    strcat('<a href="http://xenopus.hms.harvard.edu/Clustering/ClustByGeneSet/', ...
           char(SAllSets(k)),'/',char(SAllSets(k)),'_html.html">', char(num2str(SGdetected(k))),'</a>'), ...
           '</td>')));
    %
    % Total unDetected
    %
    if (SGunDetected(k) > 0)
      fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>', ...
      strcat('<a href="http://xenopus.hms.harvard.edu/Clustering/ClustByGeneSet/', ...
             char(SAllSets(k)),'/',char(SAllSets(k)),'_unDetected.txt">', char(num2str(SGunDetected(k))),'</a>'), ...
             '</td>')));
    else
      fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SGunDetected(k)),'</a>','</td>')));      
    end
    %
    % Unique in Genome
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SGuniqInGenome(k)),'</a>','</td>')));
    %
    % Unique in Detected
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SGuniqDetected(k)),'</a>','</td>')));
    %
    % Median concentration (nM)
    %
    if (SMedianConcentr(k) > 0)
      fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>', ...
      strcat('<a href="http://xenopus.hms.harvard.edu/Clustering/ClustByGeneSet/', ...
             char(SAllSets(k)),'/',char(SAllSets(k)),'_Concentration.html">', char(num2str(SMedianConcentr(k))),'</a>'), ...
             '</td>')));
    else
      fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SMedianConcentr(k)),'</a>','</td>')));
    end  
    %
    % Median V-score
    %
    if (SMedianVscore(k) > 0)
     fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>', ...
     strcat('<a href="http://xenopus.hms.harvard.edu/Clustering/ClustByGeneSet/', ...
           char(SAllSets(k)),'/',char(SAllSets(k)),'_VScore.html">', char(num2str(SMedianVscore(k))),'</a>'), ...
           '</td>')));
    else
     fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SMedianVscore(k)),'</a>','</td>')));
    end  
    %
    fprintf(fid_Main_WebPage, '%s\n', char('</tr>'));
    %
    fprintf(fid_Main_WebPage, '\n');
    fprintf(fid_Main_WebPage, '\n');
    %
    %
    %
  else
    %
    %
    fprintf(fid_Main_WebPage, '%s\n', char('<tr>'));
    %
    % Clusters
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',char(SCluster_Discription(k)),'</a>','</td>')));
    %
    % Median Dynamicity
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SDynamicity(k)),'</a>','</td>')));

    %//fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>', ...
    %//strcat('<a href="http://xenopus.hms.harvard.edu/Clustering/ClustByGeneSet/', ...
    %//       char(SAllSets(k)),'/',char(SAllSets(k)),'.html">', num2str(SDynamicity(k)),'</a>'), ...
    %//       '</td>')));

    %
    % In Set
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>', ...
    strcat('<a href="http://xenopus.hms.harvard.edu/Clustering/geneset/',char(SAllSets(k)),'.txt">', ... 
           num2str(SGinSet(k)),'</a>'), '</td>')));     
    %
    % Total in Genome
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SGinGenome(k)),'</a>','</td>')));
    %
    % Total Detected
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SGdetected(k)),'</a>','</td>')));
    %
    % Total unDetected
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SGunDetected(k)),'</a>','</td>'))); 
    %
    % Unique in Genome
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SGuniqInGenome(k)),'</a>','</td>')));
    %
    % Unique in Detected
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SGuniqDetected(k)),'</a>','</td>')));
    %
    % Median concentration (nM)
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SMedianConcentr(k)),'</a>','</td>')));
    %
    % Median V-score
    %
    fprintf(fid_Main_WebPage, '%s\n', char(strcat('<td>',num2str(SMedianVscore(k)),'</a>','</td>')));
    %
    fprintf(fid_Main_WebPage, '%s\n', char('</tr>'));
    %
    fprintf(fid_Main_WebPage, '\n');
    fprintf(fid_Main_WebPage, '\n');
    %
    %      
  end
  %
  %   
end   
%
fprintf(fid_Main_WebPage, '\n');
fprintf(fid_Main_WebPage, '\n');
%
fprintf(fid_Main_WebPage, '%s\n', char('</table>'));
fprintf(fid_Main_WebPage, '%s\n', char('</body>'));
fprintf(fid_Main_WebPage, '<h3><a href="http://xenopus.hms.harvard.edu/Clustering/ClustByGeneSet/Index_All_Zip.html" >Zip files for Clustering Data</a </h3>');
fprintf(fid_Main_WebPage, '<h3><a href="http://xenopus.hms.harvard.edu/Clustering/SourceCode/ClustByGeneSet_Final.zip">Zip files of MATLAB source code for clustering data</a </h3>');
fprintf(fid_Main_WebPage, '</body>');


fclose(fid);


disp('end');
  
