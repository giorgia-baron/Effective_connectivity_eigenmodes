function jttrend(x,varargin)
% JTTREND: Perform the Jonckheere-Terpstra test on trend.
% There are situations in which treatments are ordered in some
% way, for example the increasing dosages of a drug. In these
% cases a test with the more specific alternative hypothesis that
% the population medians are ordered in a particular direction
% may be required. For example, the alternative hypothesis
% could be as follows: population median1 <= population
% median2 <= population median3. This is a one-tail test, and
% reversing the inequalities gives an analagous test in the
% opposite tail. Here, the Jonckheere-Terpstra test can be
% used.
% Bewick V., Cheek L., Ball J. Statistics review 10: further nonparametric
% methods. Critical Care 2004, 8: 196-199
%
% Assumptions:
% - Data must be at least ordinal
% - Groups must be selected in a meaningful order i.e. ordered
% If you do not choose to enter your own group scores then scores are
% allocated uniformly (1 ... n) in order of selection of the n groups.
%
% Syntax: 	jttrend(x,score)
%      
%     Inputs:
%           X - Nx2 data matrix 
%           SCORE - order of selection of the groups
%     Outputs:
%           - Jonckheere-Terpstra statistics and p-value
%
%   Example:
% Mice were inoculated with cell lines, CMT 64 to 181, which had been
% selected for their increasing metastatic potential. The number of lung
% metastases found in each mouse after inoculation are quoted below:
%
%                                 Sample
%                   ---------------------------------
%                      64   167  170  175  181
%                   ---------------------------------
%                      0    0    2    0    2
%                      0    0    3    3    4
%                      1    5    6    5    6
%                      1    7    9    6    6
%                      2    8    10   10   6
%                      2    11   11   19   7
%                      4    13   11   56   18
%                      9    23   12   100  39    
%                           25   21   132  60
%                           97
%                   ---------------------------------
%
%       Data matrix must be:
%    d=[0 0 1 1 2 2 4 9 0 0 5 7 8 11 13 23 25 97 2 3 6 9 10 11 11 12 21 ...
%       0 3 5 6 10 19 56 100 132 2 4 6 6 6 7 18 39 60];
%    g=[ones(1,8) 2.*ones(1,10) 3.*ones(1,9) 4.*ones(1,9) 5.*ones(1,9)];
%    x=[d' g'];
%
%           Calling on Matlab the function: jttrend(x)
% (in this case, the groups are automated scored from 1 to 5)
%
%           Answer is:
%
% JONCKHEERE-TERPSTRSA TEST FOR NON PARAMETRIC TREND ANALYSIS
% --------------------------------------------------------------------------------
%     Comparison    Nx    Ny    Uxy 
%     __________    __    __    ____
% 
%     '1-2'          8    10      63
%     '1-3'          8     9    65.5
%     '1-4'          8     9      61
%     '1-5'          8     9    63.5
%     '2-3'         10     9      41
%     '2-4'         10     9    49.5
%     '2-5'         10     9    41.5
%     '3-4'          9     9    45.5
%     '3-5'          9     9      39
%     '4-5'          9     9    35.5
% 
% --------------------------------------------------------------------------------
%  
% JONCKHEERE-TERPSTRSA STATISTICS
% --------------------------------------------------------------------------------
%     Uxy_sum      JT      one_tailed_p_values
%     _______    ______    ___________________
% 
%     505        2.0116    0.022129           
% We have shown a statistically significant trend for increasing number of
% metastases across these malignant cell lines in this order.
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) Jonckheere-Terpstra test: A nonparametric Test for Trend
% http://www.mathworks.com/matlabcentral/fileexchange/22159
%Input Error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','nonempty','ncols',2}));
addOptional(p,'score',[],@(x) isempty(x) || (all(isnumeric(x(:))) && isrow(x) && all(isreal(x(:))) && all(isfinite(x(:))) && ~all(isnan(x(:))) && all(x(:)>0) && all(fix(x(:))==x(:))));
parse(p,x,varargin{:});
assert(all(x(:,2) == fix(x(:,2))),'Warning: all elements of column 2 of input matrix must be whole numbers')
score=p.Results.score;
clear p
k=max(x(:,2)); %number of groups
if isempty(score) %check score
   score=1:1:k;
end
ni=crosstab(x(:,2)); %elements for each group
N=sum(ni); %total observation
%Build-up the matrix of observation
X=NaN(max(ni),k);
for I=1:k
    X(1:ni(score(I)),I)=x(x(:,2)==score(I),1);
end
%vector and variable preallocation
A=cell(0.5*k*(k-1),4); G=1;
tr=repmat('-',1,80);% set divisor
disp('JONCKHEERE-TERPSTRSA TEST FOR NON PARAMETRIC TREND ANALYSIS')
disp(tr)
for I=1:k-1
    Uxy=zeros(1,ni(score(I)));
    for J=I+1:k
        for F=1:ni(score(I)) %for each element of the i-esim group...
            %...check how many elements of j-esim group is major or equal
            %then (in this case multiple them for 1/2)
            Uxy(F)=length(X(X(:,J)>X(F,I)))+0.5*length(X(X(:,J)==X(F,I)));
        end
        A{G,1}=sprintf('%i-%i',I,J); A{G,2}=ni(score(I)); A{G,3}=ni(score(J)); A{G,4}=sum(Uxy);
        G=G+1;
    end
end
disp(cell2table(A,'VariableNames',{'Comparison','Nx','Ny','Uxy'}))
Ut=sum(cellfun(@double,A(:,4))); N2=N^2; ni2=ni.^2;
clear A
%Compute the Jonckheere-Terpstra statistics
num=Ut-(N2-sum(ni2))/4;
denom=sqrt((N2*(2*N+3)-sum(ni2.*(2.*ni+3)))/72);
JT=abs(num/denom);
%Compare the JT stats with a standard Normal Distribution
p=1-0.5*erfc(-JT/realsqrt(2)); %p-value
disp(tr); disp(' ')
disp('JONCKHEERE-TERPSTRSA STATISTICS')
disp(tr)
disp(table(Ut,JT,p,'VariableNames',{'Uxy_sum','JT','one_tailed_p_values'}))