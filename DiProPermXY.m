function [stat,pval,zscore, L, U] = DiProPermXY(mdata1,mdata2,paramstruct) ;
% DiProPermSM, DIrection PROjection PERMutation based
% mean hypothesis test, intended for High Dimension, 
% Low Sample Size settings
%   Steve Marron's matlab function
%     Idea is to provide a 2 sample Mean Hypothesis Test,
%     based on a specified "direction" (e.g. DWD) in
%     High Dimension Low Sample Size contexts.
%     Method starts with an "important direction 
%     separating the two populations", then projects both
%     data sets onto the direction vector, and next 
%     computes a statistic comparing the projected populations.  
%     Statistical significance is assessed, using a 
%     permutation approach, where the combined data is 
%     randomly relabelled into classes (and the direction
%     is recomputed for each permutation).
%
% Inputs:
%     mdata1 - d x n1 data set 1
%     mdata2 - d x n2 data set 2
%         i.e. data vectors are columns, with same number of rows
%
%     paramstruct - a Matlab structure of input parameters
%                      Use: "help struct" and "help datatypes" to
%                           learn about these.
%                      Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%                    Version for easy copying and modification:
%     paramstruct = struct('',, ...
%                          '',, ...
%                          '',) ;
%
%    fields            values
%
%    idir             index for direction vector
%                         1  (default) DWD direction vector
%                         2  Mean Difference (aka centroid) direction vector
%                         3  Maximal Data Piling direction vector
%                         4  Fisher Linear Discrimination direction vector
%                                (based on generalized inverse)
%                         5  Support Vector Machine direction vector
%                         6  DWDLarge direction vector - a new DWD algorithm for solving large scale generalized problems                 
%
%    iper             index for permutation
%                         1  (default) all permutations
%                         2  balanced permutations
%    alpha            use 1-alpha confidence interval for z-score, (default) alpha=0.05         
%    istat            index for test statistic
%                         1  2 sample t-statistics
%                                 Warning:  will generate error for idir = 3
%                         2  (default) Difference of sample means
%                         3  Difference of sample medians
%                         4  Difference of sample medians, 
%                                 divided by pooled MADs
%                                 Warning:  will generate error for idir = 3
%                         5  Area Under the Curve (AUC), from Receiver 
%                                 Operating Characteristic (ROC) Analysis
%                                 Note:  not recommended for easily separable data
%                                 Warning:  will generate error for idir = 3
%                         6  paired sampling t-statistic
%                                 Note:  this gives scatterplot of
%                                               projections & pairwise differences
%                                 Warning:  this will generate error 
%                                               for different sample sizes
%                                 Warning:  will generate error for idir = 3
%                         7  Complement of Cluster Index (CCI)
%                                 This is 1 - Cluster Index 
%                                 (from 2-means clustering)
%                                 It is sample size adjusted squared mean 
%                                 difference divided by sum of squares about
%                                 overall mean, thus in range [0,1]
%                                 Warning:  will generate error for idir = 3
%                                     since all statistics take on value of 1
%
%    mctl             t x 2 matrix of control parameters,
%                         allowing the simultaneous running of t of these 
%                         tests, using the same permutations.
%                         For each row, the 1st entry is idir and the 
%                         2nd is istat.  Default is empty.  When non-empty, 
%                         this overrides any input values of idir & istat.
%                         When vaxh is empty, this creates a separate output 
%                         plot for each row, in Figure Windows 1,...,t.
%                         Numerical outputs are t dimensional column vectors.
%
%    vaxh              vector of 2 axis handles, useful for putting these two plots
%                      into chosen subplots.  E.g. can do:
%                          axh1 = subplot(2,2,2) ;
%                          axh2 = subplot(2,2,4) ;
%                          vaxh = [axh1 axh2] ;
%                      When mctl is not empty, this must be a t x 2 matrix 
%                          of axis handles (or left at the default of empty)
%
%    ipval            index for p-value computation:
%                         1  (default) Display both p-value (empirical quantile)
%                                and Z-score
%                         2  Show only p-value (empirical quantile)
%                         3  Show only Z-score (useful for comparisons, when
%                                            p-values are all 0)
%
%    ishowperm        0  -  Just show two panel (projections & p-value) display
%                     1  -  (default) Also show two panels of permuted 
%                                   projections (first 2, for 4 total panels)
%                           Note:  Permutes colors as well, so subdensities
%                                       show permuted differences, but heights
%                                       of points and symbols show original
%                                       class labels
%                           Note:  This is over-ridden by non-empty vaxh
%                                       or by nsim = 0 
%
%    nsim             Number of simulated relabellings to use
%                           (default is 1000)
%                     0  -  only do projection, and computation of statistic
%                               with no permutation test for signficance.
%                               (only gives the first graphical output,
%                                returns pval & zscore as empty)
%
%    nreport          How often to report permutation step taken to screen
%                     default 100 (report after 100 steps)
%                     (only has effect when iscreenwrite = 1)
%
%    seed          Seed for random number generator, use to allow
%                     identical repetition of random relabellings.
%                     Empty (default) for using the current value
%                     Note:  This uses rng version of "seed" idea.
%                         If have made call using older versions,
%                         then need to use command:  rng('default')
%                     Note:  this required modification of a number 
%                         of the DWD & SDPT3 subroutines, which were 
%                         found using a Windows Explorer search
%                         for "rand"
%
%    icolor           0  fully black and white version (everywhere)
%                     1  (default)  color version (Red for Class 1, Blue for Class 2)
%                     2x3 color matrix, top row for Class 1, bottom row for Class 2
%                           Note:  larger color matrices are deliberately not 
%                                  supported, since color provides important
%                                  cues as to what the two classes are.
%                                  To combine multiple classes (having the
%                                  same color), consider a formulation like:
%                                       unique(icolor,'rows')
%                                
%    statstrcol       Color for strings where statistics are reported
%                             Default is dark green [0 0.6 0]
%                             (this gets over-ridden to 'k' when icolor = 0)
%
%    markerstr        Can be any of:
%                         A single string with symbol to use for marker,
%                             e.g. 'o', '.', '+', 'x'
%                             (see "help plot" for a full list)
%                         A character array (2 x 1), of these symbols,
%                             One for each data set, created using:  strvcat
%                                  or using:  ['+';'o']   (default)
%                         An (n1 + n2) x 1 character array, 
%                             One for each data point
%                     Note:  When this is 'o', and ishowperm = 1,
%                            For projections of permutations 1 & 2,
%                            Use '+' for 1st permuted group and '^' for 2nd
%
%    ibigdot          0  (default)  use Matlab default for dot sizes
%                     1  force large dot size in prints (useful since some
%                              postscript graphics leave dots too small)
%                              (Caution: shows up as small in Matlab view)
%                              Only has effect on plot of simulated t-stats 
%
%    datovlaymax      maximum (on [0,1] scale, with 0 at bottom, 1 at top of plot)
%                     of vertical range for overlaid data.  Default = 0.6
%                     This applies to both the projected data, and the simulated t-stats
%
%    datovlaymin      minimum (on [0,1] scale, with 0 at bottom, 1 at top of plot)
%                     of vertical range for overlaid data.  Default = 0.5
%                     This applies to both the projected data, and the simulated t-stats
%
%    legendcellstr    cell array of strings for data set labels (2 of them),
%                     useful for (colored) classes, create this using
%                     cellstr, or {{string1 string2}}
%                     Also can indicate symbols, by just adding (at least 
%                             for +,x.o) into the text
%
%    title1str         string with title for left hand plot (showing projected data)
%                           default is empty string, '', to use:
%                               ['Projections on ' dirstr ' Direction'] 
%                           When given, and mctl is non-empty,
%                               will append [' - ' dirstr 'Dir.']
%                           For no title, use ' '
%
%    title2str         string with title for right hand plot (showing simulations)
%                           default is empty string, '', to use:
%                               [num2str(nsim) statstr ', from random relab''s'] 
%                           For no title, use ' '
%
%    titlefontsize    font size for title
%                                    (only has effect when the titlestr is nonempty)
%                           default is empty [], for Matlab default
%
%    xlabel1str        string with x axis label for left hand plot
%                           default is empty string, '', for no xlabel
%
%    xlabel2str        string with x axis label for left hand plot
%                           default is empty string, '', for no xlabel
%
%    ylabel1str        string with y axis label for left hand plot
%                           default is empty string, '', for no ylabel
%
%    ylabel2str        string with y axis label for left hand plot
%                           default is empty string, '', for no ylabel
%
%    labelfontsize     font size for axis labels
%                                    (only has effect when plot is made here,
%                                     and when a label str is nonempty)
%                           default is empty [], for Matlab default
%
%    DWDpar            DWDpenalty factor
%                         (will adjust by median pairwise dist.
%                              100 is default)
%                         (inactive, unless idir = 1)
%
%    SVMpar            SVMpenalty factor
%                         (SVMpar = C, when >= 0,
%                          SVMpar = penalty factor, when < 0
%                              (will adjust by median pairwise dist.)
%                         SVMpar = 1000 is default
%                         (inactive, unless idir = 5)
%
%    savestr          string controlling saving of output,
%                         either a full path, or a file prefix to
%                         save in matlab's current directory
%                         Will add .ps, and save as either
%                             color postscript (icolor ~= 0)
%                         or
%                             black&white postscript (when icolor = 0)
%                         Unspecified:  results only appear on screen
%                         When mctl is nonempty, and vaxh is empty,
%                             will create t postscript files, and append text:
%                                   'd#s#'
%                             where 1st # represents idir for each direction, 
%                             and 2nd # represents istat for each statistic, 
%                                 i.e. for each row of mctl
%                         When both mctl and vaxh are nonempty,
%                             will create single postscript file
%
%    iscreenwrite     0  (default)  no screen writes
%                     1  write to screen to show progress
%     
% Outputs:
%     Graphics in current Figure, 
%         showing projected data in left hand plot
%         and population of simulated pvalues in right hand plot
%     When savestr exists,
%        Postscript files saved in 'savestr'.ps
%                 (color postscript for icolor ~= 0)
%                 (B & W postscript for icolor = 0)
%     
%      stat - value of statistic, computed for data projected 
%                  onto direction vector
%
%      pval - pvalue, based on simulated quantiles.
%                 summarizing results of permutation test
%
%     zscore - Z-score summary (in permutation population),
%                 useful for comparisons, when pval = 0
%


% Assumes path can find personal functions:
%    DWD2XQ.m
%    SVM1SM.m
%    projplot1SM.m
%    cprobSM.m
%    axisSM.m
%    kdeSM.m
%    nmfSM.m
%    bwsjpiSM.m
%    lbinrSM.m
%    vec2matSM.m
%    bwrfphSM.m
%    bwosSM.m
%    rootfSM
%    vec2matSM.m
%    bwrotSM.m
%    bwsnrSM.m
%    iqrSM.m
%    cquantSM.m

%    Copyright (c) J. S. Marron 2005-2015



%  First set all parameters to defaults
%
idir = 1 ;
iper = 1;
per={'All' 'Balanced'};
alpha=0.05;
istat = 2 ;
mctl = [] ;
vaxh = [] ;
ipval = 1 ;
ishowperm = 1 ;
nsim = 1000 ;
nreport = 100 ;
seed = [] ;
icolor = 1 ;
statstrcol = [0 0.6 0] ;
markerstr = ['+';'o'] ;
ibigdot = 0 ;
idatovlay = 1 ;
datovlaymax = 0.6 ;
datovlaymin = 0.5 ;
legendcellstr = {} ;
title1str = '' ;
title2str = '' ;
titlefontsize = [] ;
xlabel1str = '' ;
xlabel2str = '' ;
ylabel1str = '' ;
ylabel2str = '' ;
labelfontsize = [] ;
DWDpar = 100 ;
SVMpar = 1000 ;
savestr = [] ;
iscreenwrite = 0 ;



%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 2 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'idir') ;    %  then change to input value
    idir = getfield(paramstruct,'idir') ; 
  end ;
  if isfield(paramstruct,'iper') ;    %  then change to input value
    iper = getfield(paramstruct,'iper') ; 
  end ;
  if isfield(paramstruct,'alpha') ;    %  then change to input value
    alpha = getfield(paramstruct,'alpha') ; 
  end ;
  if isfield(paramstruct,'istat') ;    %  then change to input value
    istat = getfield(paramstruct,'istat') ; 
  end ;

  if isfield(paramstruct,'mctl') ;    %  then change to input value
    mctl = getfield(paramstruct,'mctl') ; 
  end ;

  if isfield(paramstruct,'vaxh') ;    %  then change to input value
    vaxh = getfield(paramstruct,'vaxh') ; 
    if ~isempty(vaxh) ;
      if ~(sum(sum(ishandle(vaxh))) == (size(vaxh,1) * 2)) | ...
             ~(size(vaxh,2) == 2) ;
              %  Check that all entries are handles,
              %  and that vaxh has 2 columns
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        disp('!!!   Warning from DiProPermSM.m:     !!!') ;
        disp('!!!   Invalid vaxh,                   !!!') ;
        disp('!!!   using default of empty vaxh      !!!') ;
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        vaxh = [] ;
      end ;
    end ;
  end ;

  if isfield(paramstruct,'ipval') ;    %  then change to input value
    ipval = getfield(paramstruct,'ipval') ; 
  end ;

  if isfield(paramstruct,'ishowperm') ;    %  then change to input value
    ishowperm = getfield(paramstruct,'ishowperm') ; 
  end ;

  if isfield(paramstruct,'nsim') ;    %  then change to input value
    nsim = getfield(paramstruct,'nsim') ; 
  end ;

  if isfield(paramstruct,'nreport') ;    %  then change to input value
    nreport = getfield(paramstruct,'nreport') ; 
  end ;

  if isfield(paramstruct,'seed') ;    %  then change to input value
    seed = getfield(paramstruct,'seed') ; 
  end ;

  if isfield(paramstruct,'icolor') ;    %  then change to input value
    icolor = getfield(paramstruct,'icolor') ; 
  end ;

  if isfield(paramstruct,'statstrcol') ;    %  then change to input value
    statstrcol = getfield(paramstruct,'statstrcol') ; 
  end ;

  if isfield(paramstruct,'markerstr') ;    %  then change to input value
    markerstr = getfield(paramstruct,'markerstr') ; 
  end ;

  if isfield(paramstruct,'ibigdot') ;    %  then change to input value
    ibigdot = getfield(paramstruct,'ibigdot') ; 
  end ;

  if isfield(paramstruct,'idatovlay') ;    %  then change to input value
    idatovlay = getfield(paramstruct,'idatovlay') ; 
  end ;

  if isfield(paramstruct,'datovlaymax') ;    %  then change to input value
    datovlaymax = getfield(paramstruct,'datovlaymax') ; 
  end ;

  if isfield(paramstruct,'datovlaymin') ;    %  then change to input value
    datovlaymin = getfield(paramstruct,'datovlaymin') ; 
  end ;

  if isfield(paramstruct,'legendcellstr') ;    %  then change to input value
    legendcellstr = getfield(paramstruct,'legendcellstr') ; 
  end ;

  if isfield(paramstruct,'title1str') ;    %  then change to input value
    title1str = getfield(paramstruct,'title1str') ; 
  end ;

  if isfield(paramstruct,'title2str') ;    %  then change to input value
    title2str = getfield(paramstruct,'title2str') ; 
  end ;

  if isfield(paramstruct,'titlefontsize') ;    %  then change to input value
    titlefontsize = getfield(paramstruct,'titlefontsize') ; 
  end ;

  if isfield(paramstruct,'xlabel1str') ;    %  then change to input value
    xlabel1str = getfield(paramstruct,'xlabel1str') ; 
  end ;

  if isfield(paramstruct,'xlabel2str') ;    %  then change to input value
    xlabel2str = getfield(paramstruct,'xlabel2str') ; 
  end ;

  if isfield(paramstruct,'ylabel1str') ;    %  then change to input value
    ylabel1str = getfield(paramstruct,'ylabel1str') ; 
  end ;

  if isfield(paramstruct,'ylabel2str') ;    %  then change to input value
    ylabel2str = getfield(paramstruct,'ylabel2str') ; 
  end ;

  if isfield(paramstruct,'labelfontsize') ;    %  then change to input value
    labelfontsize = getfield(paramstruct,'labelfontsize') ; 
  end ;

  if isfield(paramstruct,'DWDpar') ;    %  then change to input value
    DWDpar = getfield(paramstruct,'DWDpar') ; 
  end ;

  if isfield(paramstruct,'SVMpar') ;    %  then change to input value
    SVMpar = getfield(paramstruct,'SVMpar') ; 
  end ;

  if isfield(paramstruct,'savestr') ;    %  then use input value
    savestr = getfield(paramstruct,'savestr') ; 
    if ~(ischar(savestr) | isempty(savestr)) ;    %  then invalid input, so give warning
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from DiProPermSM.m:    !!!') ;
      disp('!!!   Invalid savestr,               !!!') ;
      disp('!!!   using default of no save       !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      savestr = [] ;
    end ;
  end ;

  if isfield(paramstruct,'iscreenwrite') ;    %  then change to input value
    iscreenwrite = getfield(paramstruct,'iscreenwrite') ; 
  end ;


end ;    %  of resetting of input parameters


%  Initiate parameters
%
d = size(mdata1,1) ;
if ~(d == size(mdata2,1)) ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from DiProPermSM.m:        !!!') ;
  disp('!!!   mdata1 and mdata2 must have      !!!') ;
  disp('!!!   same number of rows              !!!') ;
  disp('!!!   Terminating Execution            !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  stat = [] ;
  pval = [] ;
  zscore = [] ;
  return ;
end ;
n1 = size(mdata1,2) ;
n2 = size(mdata2,2) ;
n = n1 + n2 ;
r0 = round(n1*n2/(n1+n2));
ad=[sqrt(1-(n1+n2)/(4*n1*n2-n1-n2)) sqrt(1-(n1+n2)/(4*n1*n2-2*n1-2*n2))];%% [all balanced]
mdata = [mdata1 mdata2] ;

if ~isempty(seed) ;
  rng(seed) ;
end ;

if nsim == 0 ;
  nax = 1 ;
else ;
  nax = 2 ;
end ;

if icolor == 0 ;    %  fully black and white version (everywhere)
  mcolor = zeros(n1 + n2,3) ;
  mlegcol = zeros(2,3) ;
  statstrcol = 'k' ;
      % override input choice
elseif icolor == 1 ;    %  (default)  color version (Red for Class 1, Blue for Class 2)
  mcolor = [ones(n1,1) * [1 0 0] ; ...
            ones(n2,1) * [0 0 1]] ;
  mlegcol = [[1 0 0]; [0 0 1]] ;
elseif  size(icolor,1) == 2  &  size(icolor,2) == 3  ;
                  %  have 2 rows two icolor, as required
  mcolor = [ones(n1,1) * icolor(1,:) ; ...
            ones(n2,1) * icolor(2,:)] ;
  mlegcol = icolor ;
else ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Warning from DiProPermSM:        !!!') ;
  disp('!!!   Input icolor is invalid          !!!') ;
  disp(['!!!   Need 2 rows,    input has ' num2str(size(icolor,1))]) ;
  disp(['!!!   Need 3 columns, input has ' num2str(size(icolor,2))]) ;
  disp('!!!   Replacing icolor with default    !!!') ;
  disp('!!!       and proceeding               !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  mcolor = [ones(n1,1) * [1 0 0] ; ...
            ones(n2,1) * [0 0 1]] ;
  mlegcol = [[1 0 0]; [0 0 1]] ;
  statstrcol = [0 0.6 0] ;
end ;

if size(markerstr,1) == 1 ;    %  then can use as single symbol
  vmarkerstr = markerstr ;
elseif size(markerstr,1) == 2 ;    %  then need to expand out
  vmarkerstr = [] ;
  for i = 1:n1 ;
    vmarkerstr = strvcat(vmarkerstr,markerstr(1,1)) ;
  end ;
  for i = 1:n2 ;
    vmarkerstr = strvcat(vmarkerstr,markerstr(2,1)) ;
  end ;
else ;    %  pass given markerstr through to graphics
  vmarkerstr = markerstr ;
end ;


if isempty(mctl) ;
  vidir = idir ;
  vistat = istat ;
else ;
  vidir = mctl(:,1) ;
  vistat = mctl(:,2) ;

  if ~isempty(vaxh) ;    %  both mctl and vaxh non-empty, so check sizes
    if ~((size(mctl,1) == size(vaxh,1)) & (size(mctl,2) == size(vaxh,2))) ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Error from DiProPermSM.m:         !!!') ;
      disp('!!!   mctl & vaxh must have same size   !!!') ;
      disp('!!!       (when nonempty)               !!!') ;
      disp('!!!   Terminating Execution             !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      stat = [] ;
      pval = [] ;
      zscore = [] ;
      return ;
    end ;
  end ;

end ;
nt = length(vidir) ;


if  ~isempty(vaxh)  | ...
    nsim == 0       ;

  if ishowperm == 1 ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from DiProPermSM:           !!!') ;
    disp('!!!       ishowperm = 1                   !!!') ;
    disp('!!!   cannot be used with nonempty vaxh   !!!') ;
    disp('!!!   Resetting to                        !!!') ;
    disp('!!!       ishowperm = 0                   !!!') ;
    disp('!!!   and proceeding                      !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    ishowperm = 0 ;
  end ;

end ;



%  Find full data directions, and statistics
%
stat = [] ;
for it = 1:nt ;

  iidir = vidir(it) ;
  iistat = vistat(it) ;


  %  Find Data 1 vs. Data 2 direction vector
  %
  if iidir == 1 ;    %  DWD direction vector

    dirstr = 'DWD' ;
    vdir = DWD2XQ(mdata1,mdata2,1,[],DWDpar) ;
        %  DWD direction vector, pointing from 2nd group towards first

  elseif iidir == 2 ;    %  Mean Difference (aka centroid) direction vector

    dirstr = 'MD' ;
    vdir = mean(mdata1,2) - mean(mdata2,2) ;
        %  MD direction vector, pointing from 2nd group towards first

  elseif iidir == 3 ;    %  Maximal Data Piling direction vector

    dirstr = 'MDP' ;
    vdir = MaxDatPilJA(mdata1,mdata2) ;
        %  Maximal Data Piling Vector

  elseif iidir == 4 ;    %  Fisher Linear Discrimination direction vector

    dirstr = 'FLD' ;
    vmean1 = mean(mdata1,2) ;
    vmean2 = mean(mdata2,2) ;
    mresid1 = mdata1 - vec2matSM(vmean1,size(mdata1,2)) ;
    mresid2 = mdata2 - vec2matSM(vmean2,size(mdata2,2)) ;
    mresid = [mresid1 mresid2] ;
    mcov = cov(mresid') ;
        %  Get covariance matrix, transpose, since want 
        %               "coordinates as variables"
        %  This gives "pooled within class covariance"
    mcovinv = pinv(mcov) ;
        %  pseudo-inverse
    vdir = mcovinv * (vmean1 - vmean2) ;
        %  Fisher Linear Discriminant Vector

  elseif iidir == 5 ;    %  Support Vector Machine direction vector

    dirstr = 'SVM' ;
    vdir = SVM1SM(mdata1,mdata2,SVMpar) ;
        %  SVM direction vector, pointing from 2nd group towards first
  elseif iidir == 6 ;    %  Support Vector Machine direction vector
    dirstr = 'DWD Large' ;
    DWDLarge_X = [mdata1 mdata2];
    DWDLarge_y = [ones(1, size(mdata1, 2)) -ones(1, size(mdata2, 2))]';
    options.method = 1; 
    [C,ddist] = penaltyParameter(DWDLarge_X,DWDLarge_y,1);
    vdir = genDWDweighted(DWDLarge_X,DWDLarge_y,C,1,options);
    vdir = vdir./norm(vdir);
%  DWD new algorithm direction vector, pointing from 2nd group towards first
  else ;

    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Error from DiProPermSM.m:     !!!') ;
    disp('!!!       idir is invalid           !!!') ;
    disp('!!!   Terminating execution         !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    stat = [] ;
    pval = [] ;
    zscore = [] ;
        %  set to empty, to avoid error.
    return ;

  end ;
  vdir = vdir / sqrt(sum(vdir.^2)) ;
      %  divide by length, to make it a unit vector


  %  Compute Statistic
  %
  if iistat == 1 ;    %  2 sample t-statistics
    if iidir == 3 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Error from DiProPermSM:   !!!') ;
      disp('!!!   Cannot compute t-stat,    !!!') ;
      disp('!!!        (istat = 1)          !!!') ;
      disp('!!!   for Max Data Piling Dir   !!!') ;
      disp('!!!        (idir = 3)           !!!') ;
      disp('!!!   Terminating Execution     !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      return ;
    end ;
    [h,pval,ci,stats] = ttest2(mdata1' * vdir, ...
                               mdata2' * vdir) ;
    statt = getfield(stats,'tstat') ;
    statstr = 't-stat' ;
  elseif iistat == 2 ;    %  Mean Difference
    statt = mean(mdata1' * vdir) - mean(mdata2' * vdir) ;
    statstr = 'Mean-Diff' ;
  elseif iistat == 3 ;    %  Median Difference
    statt = median(mdata1' * vdir) - median(mdata2' * vdir) ;
    statstr = 'Med-Diff' ;
  elseif iistat == 4 ;    %  MAD rescaled Median Difference
    if iidir == 3 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Error from DiProPermSM:        !!!') ;
      disp('!!!   Cannot compute MAD rescaling,  !!!') ;
      disp('!!!        (istat = 4)               !!!') ;
      disp('!!!   for Max Data Piling Dir        !!!') ;
      disp('!!!        (idir = 3)                !!!') ;
      disp('!!!   Terminating Execution          !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      return ;
    end ;
    pmdata1 = mdata1' * vdir ;
    pmdata2 = mdata2' * vdir ;
    statt = (median(pmdata1) - median(pmdata2)) / ...
                mean([madSM(pmdata1); madSM(pmdata2)]) ;
    statstr = 'Med/MAD' ;
  elseif iistat == 5 ;    %  AUC, from ROC
    if iidir == 3 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Error from DiProPermSM:   !!!') ;
      disp('!!!   Cannot compute AUC,       !!!') ;
      disp('!!!        (istat = 5)          !!!') ;
      disp('!!!   for Max Data Piling Dir   !!!') ;
      disp('!!!        (idir = 3)           !!!') ;
      disp('!!!   Terminating Execution     !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      return ;
    end ;
    [temp, statt] = ROCcurveSM(mdata1' * vdir,mdata2' * vdir,0) ;
    statstr = 'AUC-ROC' ;
  elseif iistat == 6 ;    %  Paired sampling t-statistic
    if n1 ~= n2 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Error from DiProPermSM:         !!!') ;
      disp('!!!   Cannot compute paired t-stat,   !!!') ;
      disp('!!!        (istat = 6)                !!!') ;
      disp('!!!   for different sample sizes      !!!') ;
      disp('!!!   Terminating Execution           !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      return ;
    end ;
    if iidir == 3 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Error from DiProPermSM:         !!!') ;
      disp('!!!   Cannot compute paired t-stat,   !!!') ;
      disp('!!!        (istat = 6)                !!!') ;
      disp('!!!   for Max Data Piling Dir         !!!') ;
      disp('!!!        (idir = 3)                 !!!') ;
      disp('!!!   Terminating Execution           !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      return ;
    end ;
    pmdata1 = mdata1' * vdir ;
    pmdata2 = mdata2' * vdir ;
    pmdata = pmdata1 - pmdata2 ;
    [h,pval,ci,stats] = ttest(pmdata) ;
    statt = getfield(stats,'tstat') ;
    statstr = 'Paired-t' ;
  elseif iistat == 7 ;    %  Complement Cluster Index
    if iidir == 3 ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Error from DiProPermSM:         !!!') ;
      disp('!!!   CCI is not useful               !!!') ;
      disp('!!!        (istat = 7)                !!!') ;
      disp('!!!   for Max Data Piling Dir         !!!') ;
      disp('!!!        (idir = 3)                 !!!') ;
      disp('!!!   Since all statistics take       !!!') ;
      disp('!!!       on value of 1               !!!') ;
      disp('!!!   Terminating Execution           !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      return ;
    end ;
    pmdata1 = mdata1' * vdir ;
    pmdata2 = mdata2' * vdir ;
    statt = (n1 + n2) * (mean(pmdata1) - mean(pmdata2))^2 / 4 ;
    oamean = mean([pmdata1; pmdata2]) ;
    statt = statt / (sum((pmdata1 - oamean).^2) + sum((pmdata2 - oamean).^2)) ;
    statstr = 'CCI' ;
  end ;
 stat = [stat; statt] ;

      %  save for output purposes, and for later comparisons


  %  Make Hypo test output plot
  %
  if isempty(vaxh) ;
    if ~isempty(mctl) ;
      figure(it) ;
      clf ;
    end ;
    if ishowperm == 1 ;
      subplot(2,nax,1) ;
    else ;
      subplot(1,nax,1) ;
    end ;
  else ;
    axes(vaxh(it,1)) ;
  end ;

  if isempty(title1str) ;
    projtitstr = ['Projections on ' dirstr ' Direction']  ;
  else ;
    if isempty(mctl) ;
      projtitstr = title1str ;
    else ;
      projtitstr = [title1str ' - ' dirstr ' Dir.'] ;
    end ;
  end ;


  paramstructPP1 = struct('icolor',mcolor, ...
                          'isubpopkde',1, ...
                          'markerstr',vmarkerstr, ...
                          'titlestr',projtitstr, ...
                          'titlefontsize',titlefontsize, ...
                          'xlabelstr',xlabel1str, ...
                          'ylabelstr',ylabel1str, ...
                          'labelfontsize',labelfontsize, ...
                          'datovlaymin',datovlaymin, ...
                          'datovlaymax',datovlaymax, ...
                          'iscreenwrite',iscreenwrite) ;
  if ~(isempty(legendcellstr)) ;
    paramstructPP1 = setfield(paramstructPP1,'legendcellstr',legendcellstr) ;
    paramstructPP1 = setfield(paramstructPP1,'mlegendcolor',mlegcol) ;
  end ;

  projplot1SM(mdata,vdir,paramstructPP1) ;


  %  overlay combined hypo test stat
  %
  vax = axis ;
  hold on ;
    text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
         vax(3) + 0.92 * (vax(4) - vax(3)), ...
         [statstr ' = ' num2str(statt)], ...
         'Color',statstrcol) ;
  hold off ;


end ;    %  of 1st it loop



if ~(nsim == 0) ;    %  Then do permutation test

  %  Recompute t-stats over random relabellings
  %
  mstat = zeros(nsim,nt) ;
  if ishowperm == 1 ;
    Perm12Cell = cell(9,nt) ;
        %  create 9 x nt cell array for storing info 
        %          from 1st & 2nd Permutations
  end ;
  for isim = 1:nsim ;


    if (isim / nreport) == floor(isim / nreport) ;
      if iscreenwrite == 1 ;
        disp(['    Working on sim ' num2str(isim) ' of ' num2str(nsim)]) ;
      end ;
    end ;
 flagss1sim = [ones(1,n1), zeros(1,n2)] ;
if iper==1   
    flagss1sim = [ones(1,n1), zeros(1,n2)] ;
    vunif = rand(1,n) ;
    [temp,indperm] = sort(vunif) ;
        %  indices of random permutation
elseif iper==2
     index1=randsample(n1,r0);
     index2=randsample(n2,r0);
     indperm=1:n;
     indperm(index1)=index2+n1;
     indperm(n1+index2)=index1;
else
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Error from DiProPermSM:         !!!') ;
      disp('!!!       iper must be 1 or 2         !!!') ;
      disp('!!!    selection for all/balanced     !!!') ;
      disp('!!!           permutation             !!!') ;
      disp('!!!                                   !!!') ;
      disp('!!!                                   !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;    
end
    flagss1sim = flagss1sim(indperm) ;
        %  random permutation of flagss1sim
    comboss1flag = logical(flagss1sim) ;
      
    
    

    %  Find permuted data directions, and statistics
    %
    for it = 1:nt ;

      iidir = vidir(it) ;
      iistat = vistat(it) ;

      %  Find Data 1 vs. Data 2 direction vector
      %
      if iidir == 1 ;    %  DWD direction vector

        dirstr = 'DWD' ;
        vdirc = DWD2XQ(mdata(:,comboss1flag), ...
                         mdata(:,~comboss1flag),1,[],DWDpar) ;
            %  DWD direction vector, pointing from 2nd group towards first

      elseif iidir == 2 ;    %  Mean Difference (aka centroid) direction vector

        dirstr = 'MD' ;
        vdirc = mean(mdata(:,comboss1flag),2) - mean(mdata(:,~comboss1flag),2) ;
            %  MD direction vector, pointing from 2nd group towards first

      elseif iidir == 3 ;    %  Maximal Data Piling direction vector

        dirstr = 'MDP' ;
        vdirc = MaxDatPilJA(mdata(:,comboss1flag), ...
                               mdata(:,~comboss1flag)) ;
            %  Maximal Data Piling Vector

      elseif iidir == 4 ;    %  Fisher Linear Discrimination direction vector

        dirstr = 'FLD' ;
        vmean1 = mean(mdata(:,comboss1flag),2) ;
        vmean2 = mean(mdata(:,~comboss1flag),2) ;
        mresid1 = mdata(:,comboss1flag) - vec2matSM(vmean1,size(mdata(:,comboss1flag),2)) ;
        mresid2 = mdata(:,~comboss1flag) - vec2matSM(vmean2,size(mdata(:,~comboss1flag),2)) ;
        mresid = [mresid1 mresid2] ;
        mcov = cov(mresid') ;
            %  Get covariance matrix, transpose, since want 
            %               "coordinates as variables"
            %  This gives "pooled within class covariance"
        mcovinv = pinv(mcov) ;
            %  pseudo-inverse
        vdirc = mcovinv * (vmean1 - vmean2) ;
            %  Fisher Linear Discriminant Vector

      elseif iidir == 5 ;    %  Support Vector Machine direction vector

        dirstr = 'SVM' ;
        vdirc = SVM1SM(mdata(:,comboss1flag), ...
                         mdata(:,~comboss1flag),SVMpar) ;
            %  SVM direction vector, pointing from 2nd group towards first
      elseif iidir == 6 ;    %  DWDLarge  algorithm direction vector

        dirstr = 'DWDLarge' ;
        DWDLarge_Xp = [mdata(:,comboss1flag) mdata(:,~comboss1flag)];
        DWDLarge_yp = [ones(1, size(mdata(:,comboss1flag), 2)) -ones(1, size(mdata(:,~comboss1flag), 2))]';
        options.method = 1; 
        [Cp,ddistp] = penaltyParameter(DWDLarge_Xp,DWDLarge_yp,1);
        vdirc = genDWDweighted(DWDLarge_Xp,DWDLarge_yp,Cp,1,options);
      end ;
      vdirc = vdirc / sqrt(sum(vdirc.^2)) ;
          %  divide by length, to make it a unit vector


      %  Compute Statistic
      %
      if iistat == 1 ;    %  2 sample t-statistics
        if iidir == 3 ;
          disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
          disp('!!!   Error from DiProPermSM:   !!!') ;
          disp('!!!   Cannot compute t-stat,    !!!') ;
          disp('!!!        (istat = 1)          !!!') ;
          disp('!!!   for Max Data Piling Dir   !!!') ;
          disp('!!!        (idir = 3)           !!!') ;
          disp('!!!   Terminating Execution     !!!') ;
          disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
          return ;
        end ;
        [h,pval,ci,stats] = ttest2(mdata(:,comboss1flag)' * vdirc, ...
                                   mdata(:,~comboss1flag)' * vdirc) ;
        statsim = getfield(stats,'tstat') ;
      elseif iistat == 2 ;    %  Mean Difference
        statsim = mean(mdata(:,comboss1flag)' * vdirc) - ...
                        mean(mdata(:,~comboss1flag)' * vdirc) ;
      elseif iistat == 3 ;    %  Median Difference
        statsim = median(mdata(:,comboss1flag)' * vdirc) - ...
                        median(mdata(:,~comboss1flag)' * vdirc) ;
      elseif iistat == 4 ;    %  MAD rescaled Median Difference
        if iidir == 3 ;
          disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
          disp('!!!   Error from DiProPermSM:        !!!') ;
          disp('!!!   Cannot compute MAD rescaling,  !!!') ;
          disp('!!!        (istat = 4)               !!!') ;
          disp('!!!   for Max Data Piling Dir        !!!') ;
          disp('!!!        (idir = 3)                !!!') ;
          disp('!!!   Terminating Execution          !!!') ;
          disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
          return ;
        end ;
        pmdata1 = mdata(:,comboss1flag)' * vdirc ;
        pmdata2 = mdata(:,~comboss1flag)' * vdirc ;
        statsim = (median(pmdata1) - median(pmdata2)) / ...
                    mean([madSM(pmdata1); madSM(pmdata2)]) ;
      elseif iistat == 5 ;    %  AUC, from ROC
        if iidir == 3 ;
          disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
          disp('!!!   Error from DiProPermSM:   !!!') ;
          disp('!!!   Cannot compute AUC,       !!!') ;
          disp('!!!        (istat = 5)          !!!') ;
          disp('!!!   for Max Data Piling Dir   !!!') ;
          disp('!!!        (idir = 3)           !!!') ;
          disp('!!!   Terminating Execution     !!!') ;
          disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
          return ;
        end ;
        [temp, statsim] = ROCcurveSM(mdata(:,comboss1flag)' * vdirc, ...
                                      mdata(:,~comboss1flag)' * vdirc,0) ;
      elseif iistat == 6 ;    %  Paired sampling t-statistic
        if n1 ~= n2 ;
          disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
          disp('!!!   Error from DiProPermSM:         !!!') ;
          disp('!!!   Cannot compute paired t-stat,   !!!') ;
          disp('!!!        (istat = 6)                !!!') ;
          disp('!!!   for different sample sizes      !!!') ;
          disp('!!!   Terminating Execution           !!!') ;
          disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
          return ;
        end ;
        pmdata1 = mdata(:,comboss1flag)' * vdirc ;
        pmdata2 = mdata(:,~comboss1flag)' * vdirc ;
        pmdata = pmdata1 - pmdata2 ;
        [h,pval,ci,stats] = ttest(pmdata) ;
        statsim = getfield(stats,'tstat') ;
      elseif iistat == 7 ;    %  Complement Cluster Index
        pmdata1 = mdata(:,comboss1flag)' * vdirc ;
        pmdata2 = mdata(:,~comboss1flag)' * vdirc ;
        statsim = (n1 + n2) * (mean(pmdata1) - mean(pmdata2))^2 / 4 ;
        oamean = mean([pmdata1; pmdata2]) ;
        statsim = statsim / (sum((pmdata1 - oamean).^2) + sum((pmdata2 - oamean).^2)) ;
      end ;


      mstat(isim,it) = statsim ;
          %  save for plotting, and p-value combination


      if ishowperm == 1 ;
        if isim == 1 ;
          vdirperm1 = vdirc ;
          statperm1 = statsim ;
          Perm12Cell{1,it} = vdirperm1 ;
          Perm12Cell{2,it} = statperm1 ;
          if strcmp(vmarkerstr,'o') ;
            vmarkerstr1 = [] ;
            for ipm = 1:n ;
              if comboss1flag(ipm) ;
                vmarkerstr1 = strvcat(vmarkerstr1,'+') ;
              else ;
                vmarkerstr1 = strvcat(vmarkerstr1,'^') ;
              end ;
            end ;
            Perm12Cell{3,it} = vmarkerstr1 ;
          end ;
          mcolorsim = mcolor(indperm,:) ;
          Perm12Cell{4,it} = mcolorsim ;
        elseif isim == 2 ;
          vdirperm2 = vdirc ;
          statperm2 = statsim ;
          Perm12Cell{5,it} = vdirperm2 ;
          Perm12Cell{6,it} = statperm2 ;
          if strcmp(vmarkerstr,'o') ;
            vmarkerstr2 = [] ;
            for ipm = 1:n ;
              if comboss1flag(ipm) ;
                vmarkerstr2 = strvcat(vmarkerstr2,'+') ;
              else ;
                vmarkerstr2 = strvcat(vmarkerstr2,'^') ;
              end ;
            end ;
            Perm12Cell{7,it} = vmarkerstr2 ;
          end ;

          mcolorsim = mcolor(indperm,:) ;
          Perm12Cell{8,it} = mcolorsim ;
        end ;
        Perm12Cell{9,it} = dirstr ;
      end ;


    end ;    %  of it loop inside sim loop



  end ;    %  of isim loop through simulated permutations



  %  Loop through tests once more, to get p-values & generate output graphics
  %
  pval = [] ;
  zscore = [] ;
  zscore_L = [] ;
  zscore_U = [] ;
  for it = 1:nt ;

    iistat = vistat(it) ;
    if iistat == 1 ;    %  2 sample t-statistics
      statstr = 't-stat' ;
    elseif iistat == 2 ;    %  Mean Difference
      statstr = 'Mean-Diff' ;
    elseif iistat == 3 ;    %  Median Difference
      statstr = 'Med-Diff' ;
    elseif iistat == 4 ;    %  MAD rescaled Median Difference
      statstr = 'Med/MAD' ;
    elseif iistat == 5 ;    %  AUC, from ROC
      statstr = 'AUC-ROC' ;
    elseif iistat == 6 ;    %  Paired Difference t-statistic
      statstr = 'Paired-t' ;
    elseif iistat == 7 ;    %  Complement Cluster Index
      statstr = 'CCI' ;
    end ;


    %  Compute p-values
    %
    pval = [pval; 1 - cprobSM(mstat(:,it),stat(it))] ;
        %  empirical p-value
    simmean = mean(mstat(:,it)) ;
    simsd = std(mstat(:,it)) ;
    zscore = [zscore; (stat(it) - simmean) / simsd] ;

    %%%%%%%%%%%%%%%%%% added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z_ci=[];
    for i = 1:1000
    pdf_ci=datasample(mstat(:,it),length(mstat(:,it)));
    simmean_ci = mean(pdf_ci) ;
    simsd_ci = std(pdf_ci) ;
    z_ci=[z_ci (stat(it) - simmean_ci) / simsd_ci];
    end
    zscore_L_ci=prctile(z_ci,alpha*50);
    zscore_U_ci=prctile(z_ci,100-alpha*50);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  empirical p-value
    simmean = mean(mstat(:,it)) ;
    simsd = std(mstat(:,it)) ;
 if sum(mstat(:,it)==mstat(1,it))==length(mstat(:,it));
     zscore = [zscore; 0] ;
     %%%%%%%%%%%%%%%%%%% added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      zscore_L = [zscore_L; 0] ;
      zscore_U = [zscore_U; 0] ;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 else
    zscore = [zscore; (stat(it) - simmean) / simsd] ;
    %%%%%%%%%%%%%%%%%%% added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     zscore_L = [zscore_L; zscore_L_ci] ;
     zscore_U = [zscore_U; zscore_U_ci] ;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end
 zscore=ad(iper)*zscore;
 zscore_L=ad(iper)*zscore_L;
 zscore_U=ad(iper)*zscore_U;
    
    %  add results to graphics
    %
    if isempty(vaxh) ;
      if ~isempty(mctl) ;
        figure(it) ;
      end ;
      if ishowperm == 1 ;
        subplot(2,2,2) ;
      else ;
        subplot(1,2,2) ;
      end ;
    else ;
      axes(vaxh(it,2)) ;
    end ;

    vax = axisSM([mstat(:,it); stat(it)]) ;
    if isempty(title2str) ;
      kdetitstr = [num2str(nsim) ' ' statstr 's, from ' per{iper} ' Permuted Labels'] ;
    else ;
      kdetitstr = title2str ;
    end ;

    kdeparamstruct = struct('vxgrid',vax, ...
                            'linecolor','k', ...
                            'dolcolor','k', ...
                            'ibigdot',ibigdot, ...
                            'titlestr',kdetitstr, ...
                            'titlefontsize',titlefontsize, ...
                            'xlabelstr',xlabel2str, ...
                            'ylabelstr',ylabel2str, ...
                            'labelfontsize',labelfontsize, ...
                            'datovlaymin',datovlaymin, ...
                            'datovlaymax',datovlaymax, ...
                            'iscreenwrite',1) ;
    kdeSM(mstat(:,it),kdeparamstruct) ;
    vax = axis ;
    vax(3) = 0 ;
        %  put lower end of vertical axis at 0
    axis(vax) ;
    hold on ;
      plot([stat(it); stat(it)],[vax(3); vax(4)],'Color',statstrcol) ;
      text(vax(1) + 0.3 * (vax(2) - vax(1)), ...
           vax(3) + 0.9 * (vax(4) - vax(3)), ...
           [statstr ' = ' num2str(stat(it))],'Color',statstrcol) ;
      if ipval == 1 ;
        text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
             vax(3) + 0.8 * (vax(4) - vax(3)), ...
             ['pval = ' num2str(pval(it))],'Color','k') ;
        text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
             vax(3) + 0.7 * (vax(4) - vax(3)), ...
             ['Z-score = ' num2str(zscore(it))],'Color','k') ;
        text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
              vax(3) + 0.6 * (vax(4) - vax(3)), ...
             ['Z-score CI =' '[' num2str(zscore_L(it)) ',' num2str(zscore_U(it)) ']'],'Color','k') ;

      elseif ipval == 2 ;
        text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
             vax(3) + 0.8 * (vax(4) - vax(3)), ...
             ['pval = ' num2str(pval(it))],'Color','k') ;
      elseif ipval == 3 ;
        text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
             vax(3) + 0.8 * (vax(4) - vax(3)), ...
             ['Z-score = ' num2str(zscore(it))],'Color','k') ;
         text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
              vax(3) + 0.7 * (vax(4) - vax(3)), ...
             ['Z-score CI = ' '[' num2str(zscore_L(it)) ',' num2str(zscore_U(it)) ']'],'Color','k') ;

      else ;
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        disp('!!!   Warning from DiProPermSM.m:         !!!') ;
        disp('!!!   Invalid ipval,                      !!!') ;
        disp('!!!   using default of both displays      !!!') ;
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
             vax(3) + 0.8 * (vax(4) - vax(3)), ...
             ['pval = ' num2str(pval(it))],'Color','k') ;
        text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
             vax(3) + 0.7 * (vax(4) - vax(3)), ...
             ['Z-score = ' num2str(zscore(it))],'Color','k') ;
        text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
              vax(3) + 0.6 * (vax(4) - vax(3)), ...
             ['Z-score CI = ' '[' num2str(zscore_L(it)) ',' num2str(zscore_U(it)) ']'],'Color','k') ;
   
      end ;
  L=zscore_L(it);U=zscore_U(it);    
    hold off ;


    if ishowperm == 1 ;

      dirstr = Perm12Cell{9,it} ;

      subplot(2,2,3)
      vdirperm1 = Perm12Cell{1,it} ;
      statperm1 = Perm12Cell{2,it} ;
      if strcmp(vmarkerstr,'o') ;
        vmarkerstr1 = Perm12Cell{3,it} ;
      else ;
        vmarkerstr1 = vmarkerstr ;
      end ;
      mcolorsim = Perm12Cell{4,it} ;
      paramstructPP1 = struct('icolor',mcolorsim, ...
                              'isubpopkde',1, ...
                              'markerstr',vmarkerstr1, ...
                              'titlestr',['1st Permuted ' dirstr ' Projections'], ...
                              'titlefontsize',titlefontsize, ...
                              'xlabelstr',xlabel1str, ...
                              'ylabelstr',ylabel1str, ...
                              'labelfontsize',labelfontsize, ...
                              'datovlaymin',datovlaymin, ...
                              'datovlaymax',datovlaymax, ...
                              'iscreenwrite',iscreenwrite) ;
      projplot1SM(mdata,vdirperm1,paramstructPP1) ;
      vax = axis ;
      hold on ;
        text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
             vax(3) + 0.92 * (vax(4) - vax(3)), ...
             [statstr ' = ' num2str(statperm1)], ...
             'Color','k') ;
      hold off ;

      subplot(2,2,4)
      vdirperm2 = Perm12Cell{5,it} ;
      statperm2 = Perm12Cell{6,it} ;
      if strcmp(vmarkerstr,'o') ;
        vmarkerstr2 = Perm12Cell{7,it} ;
      else ;
        vmarkerstr2 = vmarkerstr ;
      end ;
      mcolorsim = Perm12Cell{8,it} ;
      paramstructPP1 = struct('icolor',mcolorsim, ...
                              'isubpopkde',1, ...
                              'markerstr',vmarkerstr2, ...
                              'titlestr',['2nd Permuted ' dirstr ' Projections'], ...
                              'titlefontsize',titlefontsize, ...
                              'xlabelstr',xlabel1str, ...
                              'ylabelstr',ylabel1str, ...
                              'labelfontsize',labelfontsize, ...
                              'datovlaymin',datovlaymin, ...
                              'datovlaymax',datovlaymax, ...
                              'iscreenwrite',iscreenwrite) ;
      projplot1SM(mdata,vdirperm2,paramstructPP1) ;
      vax = axis ;
      hold on ;
        text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
             vax(3) + 0.92 * (vax(4) - vax(3)), ...
             [statstr ' = ' num2str(statperm2)], ...
             'Color','k') ;
      hold off ;


    end ; 


  end ;    %  of it loop for p-value calculation, and graphics generation



else ;    %  for nsim = 0

  pval = [] ;
  zscore = [] ;
      %  set to empty, to avoid error.

end ;



%  Save graphical output(s) (if needed)
%
if ~isempty(savestr) ;   %  then create postscript file(s)

  if  ~isempty(vaxh)  |  (nt == 1)  ;    %  do single page print

    orient landscape ;

    if  (size(icolor,2) > 1.5)  |  (icolor ~= 0)  ;     %  then make color postscript
      print('-dpsc',savestr) ;
    else ;                %  then make black and white
      print('-dps',savestr) ;
    end ;

  else ;

    for it = 1:nt ;

      figure(it) ;

      if ~isempty(mctl) ;
        savestradd = ['d' num2str(vidir(it)) 's' num2str(vistat(it))] ;
      else ;
        savestradd = [] ;
      end ;

      orient landscape ;

      if  (size(icolor,2) > 1.5)  |  (icolor ~= 0)  ;     %  then make color postscript
        print('-dpsc',[savestr savestradd]) ;
      else ;                %  then make black and white
        print('-dps',[savestr savestradd]) ;
      end ;

    end ;    %  of it loop for postscript file generation

  end ;

end ;




