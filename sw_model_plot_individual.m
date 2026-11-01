% Thu 11 Jun 19:03:05 CEST 2026
% Karl Kastner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% plot results of an individual model run
%%
function sw_model_plot_individual(fid,rad,t,yy,sd,sd_t,sg,stat_t, R,cva,pflag,T_avg,psl,ps)
	fflag = pflag;
	tscale = 'lin';
	meta = sw_metadata();

	cmap = [0,0,0;
                meta.colour.red;
                meta.colour.blue;
                0.4,0,0.4];

	bb = rad.extract1(yy);
	bb = single(bb);

	tdx = (t>=T_avg);

	label = rad.opt.output.base_str;
	%xscale = sg.tab.fc.star;
	normalized = true;
	if (normalized)
		xscale = sd.fc.hat;
		bscale = 0.5./sg.tab.b.thresh;
	else
		xscale = 1;
		bscale = 1;
	% TODO taba
	end

	x   = rad.x;
	lc  = 1./sd.fc.bar;
	% plot range of pattern
	Lp = 4*lc;

	% select a representative subrange of the pattern that pars-pro-toto
	% for the whole, i.e. is closest to the mean spectrun
	nx  = size(bb,2);
	nx1 = round(nx*Lp/rad.L(1));
	nx2 = floor(nx/nx1);
	b_  = reshape(bb(end,1:nx1*nx2),nx1,nx2);
	S_ = abs(fft(b_)).^2;
	% do not subtract the mean, the mean biomass matters here
	Sbar_ = mean(S_,2);
	hd2  = sum((sqrt(S_)-sqrt(Sbar_)).^2);
	[~,mdx] = min(hd2);
	fdx = nx1(1)*(mdx-1)+(1:nx1(1));
	x0 = x(fdx(1));

	% plot pattern
	splitfigure([3,2],[fid,1],fflag,num2str(fid),100);
	cla();
	l = length(sg.tab.xse);
	fdx_ = find((sg.tab.xse(:,2) >= x(fdx(1))) & (sg.tab.xse(:,1) <= x(fdx(end))));
	area(xscale*(flat([ (sg.tab.xse(fdx_,1)*[1,1])';
                            (sg.tab.xse(fdx_,2)*[1,1])' ...
                          ])-x0) ...
		    , flat((ones(length(fdx_),1)*[0,1,1,0])'), ...
		'FaceColor',meta.colour.lightgreen,'edgecolor','none');
	hold on
	plot(xscale*(x(fdx)-x0),bscale*bb(end,fdx),'color',meta.colour.blue,'linewidth',meta.linewidth)
	hold on
	xlim([0,Lp/lc])
	ylim([0,1.1*max(bscale*bb(end,fdx))]);
	if (normalized)
		ylabel(['Biomass density $b/b_t$']);
		xlabel(['Position $x/\lambda_c$']);
	else
		ylabel(['Biomass density $b/',rad.unit.b,'$']);
		xlabel(['Position $x/',rad.unit.L,'$']);
	end
	if (pflag)
		aspect = 3/4;
		pdfprint(100*fid+1,sprintf('img/sw-%s-biomass-along-transect-R-%0.2g-cva-%0.2g-psl-%g.pdf',label,R,cva,psl),ps);
	else
		title(rad.filename);
	end

	% biomass moments over time
	splitfigure([3,3],[fid,3],fflag,num2str(fid),100);
	yyaxis left
	plot(t,stat_t.b.mean);
	ylabel('Mean Biomass $\bar b$');
	yyaxis right
	plot(t,stat_t.b.sd)
	ylabel('std(b)');
	xlabel('Time t');
	set(gca,'xscale',tscale)

	% stripe and gap width distributions
	% for gaps and stripes
	type  = {'swd','gwd'}
	wtype = {'Stripe','Gap'};
	ylimscale = [1,1.5]; 
	for pdx=1:2
	% without and with fitted parametric distributions
	for qdx=1:2
	
	splitfigure([3,4],[fid,4+pdx],fflag,num2str(fid),100);
	cla();
	wscale = 1/lc;
	hp=histplot(wscale*inner2outer(sg.tab.w(pdx).xc),sg.tab.w(pdx).pdf.ksdens/wscale,'edgecolor',meta.colour.lightgrey,'facecolor',[1,1,1]*2/3);
	hold on
	f_C = {'normal','no2','qgauss'};
	l_C = {'Normal','NO_2','Q-Gauss'};
	if (qdx > 1)
	for fdx=1:length(f_C)
		plot(wscale*sg.tab.w(pdx).xc,sg.tab.w(pdx).pdf.(f_C{fdx})/wscale,'linewidth',meta.linewidth);
	end % for fdx
	end
	set(gca,'colororder',cmap);
	xlim([0,pdx*0.8])
	if (~pflag)
		legend([upper(label(1)),label(2:end)],l_C{:});
	end
	xlabel([wtype{pdx},' width $w_c/\lambda_c$']);
	ylabel('$\hspace{-1.2em}$Probability density $f\cdot{\lambda_c}$');
	set(gca,'xtick',0:0.2:0.8);
	axis square
	if (pflag)
		name = sprintf('img/sw-%s-%s-R-%0.2g-cva-%0.2g-psl-%g-with-fit-%g.pdf',label,type{pdx},R,cva,psl,qdx);
		pdfprint(fid*100+4+pdx,name,ps);
	end % if pflag
	end % for qdx
	end % for pdx

if (0)
	splitfigure([3,4],[fid,4+2],fflag,num2str(fid),100);
	cla();
	histplot(stat_a.xi,stat_a.hhg(:,end));
	xlabel('Gap width $g$');
	ylabel('Probibility density $f$');
end

	% plot swd legend
	if (1 == fid)
	for qdx=1:2
	figure(1e6);
	clf();
	plot(NaN,'k-','color',meta.colour.lightgrey,'linewidth',10);
	hold on
	if (1==qdx)
	plot(NaN(2),'-','linewidth',meta.linewidth);
	else
	plot(NaN(3),'-','linewidth',meta.linewidth);
	end
	axis off
	%legend(label,l_C{:});
	set(gca,'colororder',cmap);
	legend([upper(label(1)),label(2:end)],l_C{:});
	if (pflag)
		pdfprint(1e6,sprintf('img/sw-%s-swd-legend-%d.pdf',label,qdx),ps);
	end % if pflag
	end % if qdx
	end

	splitfigure([3,4],[fid,4+3],fflag,num2str(fid),100);
	cla();
try
	fx = sd.fourier_axis();
	fc = sd.fc.bar;
	fdx=(fx>=0 & fx <= 4*fc);
	plot(fx(fdx)/fc,fc*sd.S.bar(fdx),'color',meta.colour.blue,'linewidth',meta.linewidth);
	switch(label)
	case {'tlidi'}
	text(1.8,0.9*fc*sd.Sc.bar,sprintf('$\\lambda_c$=%2d\\,[L]',round(1./fc)),'interpreter','latex');
	case {'rietkerk'}
	text(1.8,0.9*fc*sd.Sc.bar,sprintf('$\\lambda_c$=%2d\\,m',round(1./fc)),'interpreter','latex');
	end
	set(gca,'xtick',0:4);
	set(gca,'ytick',0:0.2:2);
	xlim([0,3.5])
	ylim([0,1.1*sd.Sc.bar*fc]);
	axis square
	xlabel('Wavenumber $k/k_c$');
	ylabel('$\hspace{-1em}$Spectral density $S_{x}^+/\lambda_c$')
	if (pflag)
		name = sprintf('img/sw-%s-spectral-density-R-%0.2g-cva-%0.2g-psl-%g.pdf',label,R,cva,psl);
		pdfprint(fid*100+4+3,name,ps);
	end % if pflag
catch
end

	splitfigure([3,4],[fid,2*4+1],fflag,num2str(fid),100);
	yyaxis left
	cla
	plot(t,stat_t.w(1).mean);
	ylabel('Mean stripe width $\bar w$');
	yyaxis right
	cla
	plot(t,stat_t.w(1).sd)
	ylabel('std(w)');
	xlabel('Time $t$');
	set(gca,'xscale',tscale)
	
	splitfigure([3,4],[fid,2*4+2],fflag,num2str(fid),100);
	yyaxis left
	plot(t,[stat_t.w(1).sk])
	ylabel('Skewness(w)')
	yyaxis right
	plot(t,[stat_t.w(1).ku])
	ylabel('Kurtosis(w)')
	xlabel('Time $t$');
	set(gca,'xscale',tscale)

	% regularity over time
	splitfigure([3,4],[fid,2*4+3],fflag,num2str(fid),100);
	yyaxis left
	cla();
	try
	plot(t,[stat_t.Sc.bar.*stat_t.fc.bar,stat_t.Sc.bar.*stat_t.fc.bar, ...
			stat_t.Sc.hat.*stat_t.fc.hat])
	ylabel('Regularity $S_{xc}^+/\lambda_c$');
	xlabel('Time $t$');
	set(gca,'xscale',tscale)
	yyaxis right
	cla();
	plot(t,1./stat_t.fc.bar)
	catch e
		e
	end 

	splitfigure([3,4],[fid,2*4+4],fflag,num2str(fid),100);
	subplot(3,5,15)
	plot(t,stat_t.celerity)
	ylabel('Celerity $c$');
	xlabel('Time $t$');
	set(gca,'xscale',tscale)
	
end % sw_tlidi_plot_1

