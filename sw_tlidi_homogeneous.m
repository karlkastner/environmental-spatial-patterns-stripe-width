% Mon 13 Jul 12:22:47 CEST 2026
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
%%function sw_tlidi_homogeneous_convergence(pflag)
%%
%% generate a pattern with the Tlidi model in homogenous coordinate and plot
%% its line spectrum and degenerate SWD
function [sd_t,sg_t,stat_t] = sw_tlidi_homogeneous(fid,pflag)
	if (nargin()<1||isempty(fid))
		fid = 1;
	end % if nargin()<1
	if (nargin()<2)
		pflag = 0;
	end % if nargin()<2
	fflag = pflag;
	plotscale = 5;
	init_plot();
	meta = sw_metadata();
	maxNumCompThreads(meta.maxNumCompThreads);
	
	% spatial extent
	L   = 1024;

	% spatial resolution
	dx = 1;

	% duration of model run (in simulated time, not cpu runtime)
	T = 4e3*365;

	% start time for avraging pattern properties
	T_avg = 0.99*T;

	% mortality rate, determining relative water availability
	mu = 0.7;
	
	% degree of exogenous heterogeneity
	% spatial coefficient of variation
	cva = 0;

	% correlation length
	psl = 0;

	% SWD histogram bin edges
	xi = (0:1:100);
	
	cv = [];
	tlimit = NaN(length(dx),2);
	for idx=1:length(dx)
		dt = 0.1*dx(idx);
		
		% base model configuration
		param = sw_tlidi_param(mu,cva,psl,T,L,dx(idx));
		param.opt.rng = 0;

		param.initial_condition = 'obj.random_state()';

		% create output folder
		mkdir(param.opt.output.path_str);

		rad = Tlidi(param);
		% simulate
		[t,bb,out] = rad.run();
		bb = double(bb);
	
		b = cvec(bb(end,:));	
		[bmin,mdx]=min(b);
		b=circshift(b,-mdx);
		
		[sd,sd_t,sg,sg_t,stat_t] = sw_model_analyze(rad,t,bb,xi,T_avg);
		
		% water-availability proxy
		R = (1-param.pmu.mu/param.pmu.ef);
		% quick fix: the sd esimate fails for the line spectrum, so override with unsmoothed value
		sd.fc.bar = sd.fc.hat;
		sw_model_plot_individual(fid,rad,t,bb,sd,sd_t,sg,stat_t,R,cva,pflag,T_avg,psl,plotscale);
		
		fx = sd_t(end).fourier_axis;
		fc = sd_t(end).fc.hat;
		S  = sd_t(end).S.hat;
		Sc = sd_t(end).Sc.hat;
		xscale = fc;

		% determine time where sd of stripe width first drops below dx
		sd    = stat_t.w(1).sd;
		sd(1) = inf;
		fdx         = find(sd<dx,1,'first');
		if (~isempty(fdx))
		tlimit(idx,1) = t(fdx);
		end
		wrange      = arrayfun(@(x) range(x.tab.w(1).w), sg_t);
		fdx         = find(wrange<dx,1,'first');
		if (~isempty(fdx))
		tlimit(idx,2) = t(fdx);
		end 

		% biomass pattern along transect
		splitfigure([2, 4], [fid+1, 1], fflag);
		cla();
		plot(rad.x,b)
		hold on
		xlabel('Distance $x$/m');
		ylabel('Biomass density $b/[b]$');
		
		splitfigure([2, 4], [fid+1, 2], fflag);
		plot(t,mean(bb,2))
		xlabel('Time $t$/d');
		ylabel('Biomass denisty \bar b')
		
		splitfigure([2, 4], [fid+1, 3], fflag);
		histplot(xscale*sg_t(end).tab.w(1).xe, sg_t(end).tab.w(1).pdf.hist/xscale,'facecolor',meta.colour.blue,'edgecolor',meta.colour.blue)
		xlabel('Stripe width $w / \lambda_c$');
		ylabel('Probability density $f \cdot \lambda_c$');
		xlim([0,0.6])
		set(gca,'xtick',0:0.2:1);
		axis square
		hmax = max(sg_t(end).tab.w(1).pdf.hist/xscale);
  	        ylim([0,1.1*hmax]);
		text(0.31,0.8*hmax,sprintf('$\\displaystyle\\frac{\\bar w}{\\lambda_c}\\!=$%0.2g', sg_t(end).tab.w(1).mean*fc),'interpreter','latex');
		
		splitfigure([2, 4], [fid+1, 4], fflag);
		stem(fx/fc,S*fc,'linewidth',1,'marker','none','color',meta.colour.blue);
		xlim([0,3.5]);
		yl=([0,max(Sc*fc)*1.1]);
		ylim(yl);
		ylabel('Spectral density $S_x^+/\lambda_c$');
		xlabel('Wavenumber $k/k_c$');
		xlim([0,3.5]);
		set(gca,'xtick',0:4);
		text(1.75,0.8*Sc*fc,sprintf('$\\lambda_c$=%0.3g{m}', 1./fc),'interpreter','latex');
		axis square
		
		splitfigure([2, 4], [fid+1, 5], fflag);
		d = diff(bb)./diff(t);
		d = rms(d,2);
		plot(mid(t),d);
		ylabel('Rate of change ||db/dt||');
		
		splitfigure([2, 4], [fid+1, 6], fflag);
		semilogx(t,stat_t.Sc.hat.*stat_t.fc.hat);
		xlabel('Time t/d');
		ylabel('Regularity $S_c/\lambda_c$')
		
		sk=[];
		ku=[];
		fullpop=true;
		for jdx=1:size(bb,1)
		l = analyze_stripe_width_simple(bb(jdx,:),dx(idx));
		cv(jdx,idx) = std(l,fullpop)./mean(l);
		sk(jdx,1) = skewness(l,fullpop);
		ku(jdx,1) = kurtosis(l,fullpop);
		end
		
		splitfigure([2, 4], [fid+1, 7], fflag);
		yyaxis left
		cla();
		plot(t,cv)
		yyaxis right
		cla();
		plot(t,sk)
		hold on
		plot(t,ku)
		legend('cv','sk','ku')
		xlabel('Time t/d');
		
	
	end % for idx (dx)

	printf('Regularity %g\n',fc*Sc);
	printf('Number of stripes %g\n',fc*L);
	printf('Limit sd time %g years\n',tlimit(1)/365);
	printf('Limit range time %g years\n',tlimit(2)/365);
	%printf('Stripe width limits %g %g\n',wlimit);
	
	if (pflag)
		base_str = rad.opt.output.base_str;
		pdfprint(10*(fid+1)+3,['img/sw-',base_str,'-homogenous-spectral-density.pdf'],plotscale);
		pdfprint(10*(fid+1)+4,['img/sw-',base_str,'-homogenous-stripe-width-pdf.pdf'],plotscale);
	end
	
	
	splitfigure([2,2],[fid+2,1],fflag);
	semilogy(t,cv./dx.^2,'linewidth',meta.linewidth);
	ylabel('Coefficient of variation $cv_w/(\Delta{x})^2$')
	xlabel(['Time $t / ',rad.unit.T,'$']);
	lh=legend(arrayfun(@num2str,dx,'uniformoutput',false));
	title(lh,['$\Delta{x}/',rad.unit.L,'$'],'interpreter','latex');
	set(gca,'colororder',[0,0,0;meta.colour.red;meta.colour.blue]);
	axis square
	set(gca,'ytick',10.^(-15:0))
	
	splitfigure([2,2],[fid+2,2],fflag);
	semilogy(t,cv,'linewidth',meta.linewidth);
	set(gca,'colororder',[0,0,0;meta.colour.red;meta.colour.blue])
	ylabel('Coefficient of variation $cv_w$')
	xlabel(['Time $t / ',rad.unit.T,'$']);
	lh=legend(arrayfun(@num2str,dx,'uniformoutput',false));
	title(lh,['$\Delta x / ',rad.unit.L,'$'],'interpreter','latex')
	axis square
	set(gca,'ytick',10.^(-15:0))
	set(gca,'ylim',10.^-[3.3,0]);
	xlim([0,3e4]);
	
	if (pflag)
		pdfprint(10*(fid+2)+2,'img/sw-tlidi-swd-vs-t-homogeneous.pdf',plotscale);
	end
end % sw_tlidi_homogenous_convergence

