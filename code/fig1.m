% fig3TestFun - produces Figure 3. of the paper
%
% Efficiency of the V-fold model selection for localized bases
%
%   Copyright (c) 2017 Fabien Navarro and Adrien Saumard

clear all
close all

n = 4096;
signal_names = {'Angles','Parabolas','Corner'};       
              
for in=1:length(signal_names)
  figure(in);
  set(gcf,'Name',signal_names{in},'NumberTitle','off');
  [X,~,s_star] = makedata(n,signal_names{in});
  plot(X,s_star,'k','LineWidth',2)
  xlim([0 1]);ylim([0 1])
end
tilefigs