% 2026-07-13 23:16:11.988460421 +0200
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
%% function sw_rk1d_homogeneous_convergence(pflag)
%%
%% generate a pattern with the Rietkerk model in homogenous coordinate and plot
%% its line spectrum and degenerate SWD
function [sd_t,sg_t,stat_t] = sw_rk1d_homogeneous_convergence(fid,pflag)
	if (nargin()<1||isempty(fid))
		fid = 1;
	end % if nargin()<1
	if (nargin()<2)
		pflag = 0;
	end % if nargin()<2
	fflag = pflag;
	plotscale    = 5;
	init_plot();
        meta = sw_metadata();
	maxNumCompThreads(meta.maxNumCompThreads);
	
	% spatial extent
	L     = 1024;

	% spatial resolution
	dx    = 1;

	% duration of model run (in simulated time, not cpu runtime)
	T     = 200*365;
	
	% start time for avraging pattern properties
	T_avg = 0.99*T;

	% precipitation rate
	R = 0.8;

	% degree of exogenous heterogeneity
	% spatial coefficient of variation of the infiltration coefficient
	cva = 0;

	% correlation length
	psl = 0;

	% SWD histogram bin edges
	xi = (0:1:100);

	% declare output vars
	tlimit = NaN(length(dx),2);
	for idx=1:length(dx)
		% base model configuration
		param = sw_rk1d_param(R,cva,psl,T,L,dx);

		% create output folder
		mkdir(param.opt.output.path_str);

		% run specific parameters
	        % data type during computtion file
		% high number of digits here to check convergence for dx -> 0
		param.opt.compute_class = @double;
		% data type of output file
		param.opt.output.class  = @double;
	
		rad = Rietkerk(param);

		% simulate
		[t,y,out_rad]	  = rad.run();
		y = double(y);
		[b,w,h] = rad.extract1(y);
		bb = rad.extract1(y);
		b=bb(end,:);
		[bmin,mdx]=min(b);
		b=circshift(b,-mdx);
	
		rad.aux.dt = 0;
		[sd,sd_t,sg,sg_t,stat_t] = sw_model_analyze(rad,t,y,xi,T_avg);
	
		cvw = stat_t.w(1).sd./stat_t.w(1).mean;
		
		% quick fix: the sd esimate fails for the line spectrum, so override with unsmoothed value
		sd.fc.bar = sd.fc.hat;
		sw_model_plot_individual(fid,rad,t,y,sd,sd_t,sg,stat_t,R,cva,pflag,T_avg,psl,plotscale);
	
		fx = sd_t(end).fourier_axis;
		fc = sd_t(end).fc.hat;
		S  = sd_t(end).S.hat;
		Sc = sd_t(end).Sc.hat;
		xscale = fc;

		% determine time where sd of stripe width first drops below dx
		sd = stat_t.w(1).sd;
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
		ylabel('Biomass density $b$/(g/m$^2$)');
		
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
	end % for idx

	printf('Regularity %g\n',fc*Sc);
	printf('Number of stripes %g\n',fc*L);
	printf('Limit sd time %g years\n',tlimit(1)/365);
	printf('Limit range time %g years\n',tlimit(2)/365);

	splitfigure([2,4],[fid+1,5],fflag);
	semilogx(t/365,stat_t.Sc.hat.*stat_t.fc.hat);
	xlabel('Time $t$/year');
	ylabel('Regularity $S_c/\lambda_c$')
		  
	splitfigure([2,4],[fid+1,6],fflag);
	semilogy(t/365,cvw,'linewidth',meta.linewidth);
	ylabel('Coefficient of variation $cv_w$');
	xlabel('Time $t$ / y');
	leg_C = arrayfun(@(x) num2str(x),dx,'uniformoutput',false);
	if (length(dx)>1)
		lh=legend(leg_C)
		title(lh,'$\Delta{x}$ / m','interpreter','latex');
	end
	set(gca,'colororder',[0,0,0;meta.colour.red;meta.colour.blue])
	xlim([0,500])
	set(gca,'ytick',10.^(-15:0))
	axis square
	yyaxis right
	semilogy(t/365,cvw./dx.^2,'linewidth',meta.linewidth);
	ylabel('cv(w)/dx$^2$')
	
	
	if (pflag)
		base_str = rad.opt.output.base_str;
		pdfprint(10*(fid+1)+3,['img/sw-',base_str,'-homogenous-spectral-density.pdf'],plotscale);
		pdfprint(10*(fid+1)+4,['img/sw-',base_str,'-homogenous-stripe-width-pdf.pdf'],plotscale);
		pdfprint(10*(fid+1)+5,'img/sw-rk-cvw-vs-t-homogeneous.pdf',plotscale);
	end
end % sw_rk1d_homogeneous_convergence

