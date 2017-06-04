% fig2 - produces Figure 2. of the paper
%
% Efficiency of the V-fold model selection for localized bases
%
% Copyright (c) 2017 Fabien Navarro and Adrien Saumard

clear all
close all

signal_name = {'Angles','Parabolas','Corner'};
Pn = @(x0,x)(mean((x-x0).^2,2));

n = 4096;
J = log2(n);
j0 = 0;
D = 2.^(1:J-1);
filter = 'Daubechies';
qmf = MakeONFilter(filter,8);

for in=1:length(signal_name)
  rand('seed',0)
  randn('seed',0)
  hat_s_m = zeros(length(D),n);
  wc_hat_s_m = zeros(1,n);
  [X,Y,s_star] = makedata(n,signal_name{in});
  wc = FWT_PO(Y,j0,qmf);
    
  n1 = n/2;
  wc_odd_lin  = zeros(1,n1);
  wc_even_lin = zeros(1,n1);
  Y_odd = Y(1:2:end);
  Y_even = Y(2:2:end);
      
  wc_odd  = FWT_PO(Y_odd,j0,qmf);
  wc_even = FWT_PO(Y_even,j0,qmf);
  hat_s_odd  = zeros(length(D),n1);
  hat_s_even = zeros(length(D),n1);
      
  for ii =1:length(D)
    wc_hat_s_m(1:2^ii) = wc(1:2^ii);
    hat_s_m(ii,:) = IWT_PO(wc_hat_s_m,j0,qmf);
      
    wc_odd_lin(1:2^(ii))  = wc_odd(1:2^(ii));
    wc_even_lin(1:2^(ii)) = wc_even(1:2^(ii));
    hat_s_odd(ii,:)  = IWT_PO(wc_odd_lin,j0,qmf);
    hat_s_even(ii,:) = IWT_PO(wc_even_lin,j0,qmf);
  end
  err = Pn(hat_s_m,repmat(s_star,length(D),1));
  [~,m_star] = min(err);
  hat_s_m_star = hat_s_m(m_star,:);
        
  % 2-fold CV
  mat_Y = repmat(Y,length(D),1);
  bar_s_odd = 0.5.*(hat_s_odd(:,1:n1-1) + hat_s_odd(:,2:n1));
  bar_s_odd(:,n1) = hat_s_odd(:,1);
  bar_s_even = 0.5.*(hat_s_even(:,1:n1-1) + hat_s_even(:,2:n1));
  bar_s_even(:,n1) = hat_s_even(:,1);
      
  mat_Y_odd = repmat(Y_odd,length(D),1);
  mat_Y_even = repmat(Y_even,length(D),1);
  Crit_CV = sum( (mat_Y_odd-bar_s_even).^2+...
                 (mat_Y_even-bar_s_odd).^2,2 );
  [~,hat_m_2FCV] = min(Crit_CV);
  hat_s_m_2FCV = hat_s_m(hat_m_2FCV,:);
      
  % 2-CVpen
  Crit_pen2F = 1/n*sum((hat_s_m-mat_Y).^2,2)+...
               1/(2*n)*(-sum((mat_Y_even- hat_s_even).^2,2)...
                        +sum((mat_Y_odd - bar_s_even).^2,2)...
                        -sum((mat_Y_odd - hat_s_odd).^2,2)...
                        +sum((mat_Y_even- bar_s_odd).^2,2) );              
  [~,hat_m_pen2F] = min(Crit_pen2F);
  hat_s_m_pen2F = hat_s_m(hat_m_pen2F,:);
   
  figure;handle_axis = gca;set(handle_axis,'FontSize', 20);
  set(gcf,'Name',[signal_name{in}],'NumberTitle','off');
  plot(X,Y,'k','LineWidth',2);
  xlim([0 1]);ylim([-0.1 1.1]);
  xlabel('$X$','Interpreter','latex');
  ylabel('$Y$','Interpreter','latex');
    
  figure;handle_axis = gca;set(handle_axis,'FontSize', 20);
  plot(X,s_star,'k--','LineWidth',2);hold on
  plot(X,hat_s_m_2FCV,'k','LineWidth',2);
  plot(X,hat_s_m_pen2F,'k-.','LineWidth',2);
  legend('$s_*$','$\widehat{s}_{\widehat{m}_\mathrm{2FCV}}$',...
         '$\widehat{s}_{\widehat{m}_\mathrm{pen2F}}$',...
         'Location','NorthEast','Orientation','Vertical');
  legend boxoff;
  h0 = legend;set(h0, 'interpreter', 'latex','FontSize',30)
  xlim([0 1]);ylim([0 1]);
    
  figure;handle_axis = gca;set(handle_axis,'FontSize', 20);
  rCrit_CV = rescale(Crit_CV,min(err),max(err));
  rCrit_pen2F = rescale(Crit_pen2F,min(err),max(err));
  loglog(D,err,'k-o','LineWidth',1.2);grid;hold on;axis tight
  h = loglog(D,rCrit_CV,'k-o','LineWidth',1.2);
  set(h, 'color', [0.5 0.5 0.5])
  hp = loglog(D,rCrit_pen2F,'k-o','LineWidth',1.2);
  set(hp, 'color', [0.75 0.75 0.75])   
  loglog(D(hat_m_2FCV),rCrit_CV(hat_m_2FCV),'o','LineWidth',2,...
         'MarkerEdgeColor','k','MarkerFaceColor',[.5 .5 .5],...
         'MarkerSize',24);
  loglog(D(hat_m_pen2F),rCrit_pen2F(hat_m_pen2F),'d','LineWidth',2,...
         'MarkerEdgeColor','k','MarkerFaceColor',[.75 .75 .75],...
         'MarkerSize',18);
  loglog(D(m_star),err(m_star),'k*','LineWidth',2,'MarkerSize',12)    
  xlabel('$D_m$', 'interpreter', 'latex');
  legend('$\ell(s_{\ast},\widehat{s}_{m})$',...
         '$\mathrm{crit_{\mathrm{2FCV}}}(m)$',...
         '$\mathrm{crit_{\mathrm{pen2F}}}(m)$',...
        ['$D_{\widehat{m}_\mathrm{2FCV}}=$' num2str(2^hat_m_2FCV)],...
        ['$D_{\widehat{m}_\mathrm{pen2F}}=$' num2str(2^hat_m_pen2F)],...
        ['$D_{m_{\ast}}=$' num2str(2^m_star)],...
         'Location','NorthEast','Orientation','Vertical');
  h1 = legend;set(h1, 'interpreter', 'latex','FontSize',22)
end

tilefigs
