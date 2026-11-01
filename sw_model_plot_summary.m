% 2026-08-03 16:08:23.468421720 +0200
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
%% function sw_model_heterogeneous_plot(base_str,fid,pflag,ps)
%% 2D-surface plots of control variables against the two control variables
%%
%% input:
%%	base_str : model identifier string, either "rietkerk" or "tlidi"
%%	fid      : figure index
%%	pflag    : (boolean) set to true for writing plots to files
%%	ps	 : (scalar double) scale of plots
function sw_model_plot_summary(path_str,fid,pflag)
	if (nargin()<2)
		fid = 1;
	end
	if (nargin()<3)
		pflag = false;
	end
	% note, colormaps in matlab are figure-wide, for subplot individual
	% colormaps set fflag to 1 to plot each surface into an own figure
	fflag = pflag;	
	init_plot();
	meta = sw_metadata();

	% smooting window width
	ns = 3;
	
	pfun = @contourf;
	
	% load output
	load([path_str,filesep(),'summary.mat'],'o');
	base_str = o.model;
	%load([path_str,filesep(),base_str,'.mat'],'o');

	switch (base_str)
	case {'rietkerk'}
	plotscale = 4.5;
	otherwise
	% because the label has a fraction the vertical size is larger, zoom needs to be adapted
	plotscale = 4.1;
	end

	% select the x-label, corresponding to the respective model
	switch (base_str)
	case {'rietkerk'}
		R          = o.R;
		xlabel_str = '\hspace{1.3em}Precipitation $R$/(mm/d)';
	case {'tlidi'}
		R          = o.water_availability;
		xlabel_str = '\hspace{2.5em}Water\,availability\,$\left(1\!-\!\frac{\mu}{a}\right)$';
	otherwise
		error('sw_model_plot_summary');
	end
	cva = o.cva;
	psl = o.psl;
	
	% for each correlation length
	for psldx=1:length(psl)
	
	splitfigure([2, 5], [fid+psldx-1, 1], fflag);
	cla();
	level.tlidi = 0:0.05:0.8;
	level.rietkerk    = 0:0.10:0.8;
	pfun(R,cva,smooth(o.w(1).mean(:,:,psldx).*o.fc.bar(:,:,psldx),ns),level.(base_str));
	view(0,90);
	ch=colorbar();
	axis tight
	axis square
	axis xy;
	ylabel('$\hspace{-0.8em}$Exogenous heterogeneity $\mathrm{cv}(e)$');
	xlabel(xlabel_str);
	title('Stripe width $\bar w/\lambda_c$')
	if (strcmp(base_str,'tlidi'))
		caxis([0.15,0.55]);
		colormap(parula(8))
		set(ch,'xtick',0.2:0.1:0.5);
	else
		caxis([0,0.8])
		colormap(parula(8))
	end
	
	splitfigure([2, 5], [fid+psldx-1, 2], fflag);
	level.tlidi = 0:0.05:0.7;
	level.rietkerk    = 0:0.1:0.7;
%	if (strcmp(base_str,'tlidi'))
%		level = 0:0.05:0.7;
%	else
%		level = 0:0.1:0.7;
%		
%	end
	pfun(R,cva,smooth(o.w(1).sd(:,:,psldx)./o.w(1).mean(:,:,psldx),ns),level.(base_str));
	view(0,90);
	colorbar();
	axis tight;
	ylabel('$\hspace{-0.8em}$Exogenous heterogeneity $\mathrm{cv}(e)$');
	xlabel(xlabel_str);
	axis square
	axis xy;
	title('\hspace{3em}Coefficient\,of\,variation\,$\mathrm{cv}{(w)}$');
	if (strcmp(base_str,'tlidi'))
		caxis([0, 0.55])
		colormap(parula(11))
	else
		caxis([0,0.7])
		colormap(parula(7))
	end
	
	
	splitfigure([2, 5], [fid+psldx-1, 3], fflag);
	cla;
	level.tlidi = -1.4:0.20:1.6;
	level.rietkerk    = -1.0:0.25:1.5;
	caxis_.tlidi = [-1,0];
	caxis_.rietkerk    = [-1,1.5];
	pfun(R,cva,smooth(o.w(1).sk(:,:,psldx),ns),level.(base_str));
	view(0,90);
	ch=colorbar();
	axis tight
	axis square
	axis xy;
	xlabel(xlabel_str);
	ylabel('$\hspace{-0.8em}$Exogenous heterogeneity $\mathrm{cv}(e)$');
	title('Skewness $\mathrm{sk}(w)$');
	if (strcmp(base_str,'tlidi'))
		colormap(parula(10))
		set(ch,'xtick',-1:0.2:0);
	else
		colormap(parula(10))
	end
	caxis(caxis_.(base_str))

	splitfigure([2, 5], [fid+psldx-1, 4], fflag);
	cla();
	% bimodality
	if (1)
		ns_ = ns;
        	bimodality_c = smooth(o.w(1).bimodality_c,ns_);
        	bimodality_small = smooth(o.w(1).bimodality_small,ns_);
        	bimodality_large = smooth(o.w(1).bimodality_large,ns_);
        	bimodality_rel   = 1 - bimodality_c./bimodality_small;
	        % ignore too small lobes                                                
        	 bimodality_rel(bimodality_small < 0.05*bimodality_large)  = 0;
	end
	%r = o.w(1).bimodality_rel;
	r = bimodality_rel;
	level = 0:0.1:0.5;
	ns_ = 1;
	pfun(R,cva,r,level);
	%pfun(R,cva,r,level);
	view(0,90);
	axis square
	axis xy;
	ch=colorbar();
	xlabel(xlabel_str);
	ylabel('$\hspace{-0.8em}$Exogenous heterogeneity $\mathrm{cv}(e)$');
	title('Bimodality $\mathrm{d}(w)$');
	if (strcmp(base_str,'tlidi'))
		caxis([0, 0.4])
		colormap(parula(8))
	else
		colormap(parula(10))
		caxis([0 0.5])
	end

	splitfigure([2, 5], [fid+psldx-1, 5], fflag);
	level = -1.75:0.25:2;
	pfun(R,cva,smooth(log2(o.Sc.bar(:,:,psldx).*o.fc.bar(:,:,psldx)),ns),level);
	view(0,90);
	ch = colorbar();
	axis tight
	axis square
	axis xy;
	xt=(-2:0.5:5);
	set(ch,'xtick',xt,'xticklabel',roundn(2.^xt,2));
	axis square
	axis xy;
	xlabel(xlabel_str);
	ylabel('$\hspace{-0.8em}$Exogenous heterogeneity $\mathrm{cv}(e)$');
	title('Regularity $S_{xc}^{+}/\lambda_c$')
	if (strcmp(base_str,'tlidi'))
		caxis([-1,1.25]);
		colormap(parula(9))	
	else
		caxis([-1.75 1])
		colormap(parula(11))
	end
	
	splitfigure([2, 5], [fid+psldx-1, 6], fflag);
	cla();
	level = struct();
	level.tlidi = 4.75:0.06125:8;
	level.rietkerk    = 4.75:0.25:7;
	pfun(R,cva,smooth(log2(1./o.fc.bar(:,:,psldx)),ns),level.(base_str));
	view(0,90);
	ch=colorbar();
	axis tight
	axis square
	axis xy;
	axis square
	axis xy;
	xlabel(xlabel_str);
	ylabel('$\hspace{-0.8em}$Exogenous heterogeneity $\mathrm{cv}(e)$');
	title('Wavelength $\lambda_c$/m');
	if (strcmp(base_str,'tlidi'))
		xt=(0:0.125:8);
		set(ch,'xtick',xt,'xticklabel',round(2.^xt));
		caxis([5.75+0.125+0.0625,6.5+0.06125])
		colormap(parula(10))
	else
		xt=(0:0.5:8);
		set(ch,'xtick',xt,'xticklabel',round(2.^xt));
		caxis([5.25,7]);
		colormap(parula(7))
	end

	splitfigure([2, 5], [fid+psldx-1, 7], fflag);
	cla();
	sig = median(sign(o.celerity),'all');
	pfun(R,cva,sig*smooth(o.celerity.*o.fc.bar,ns));
	view(0,90);
	ch=colorbar();
	axis xy;
	axis tight
	axis square
	xlabel(xlabel_str);
	ylabel('$\hspace{-0.8em}$Exogenous heterogeneity $\mathrm{cv}(e)$');
	switch (base_str)
	case {'rietkerk'}
		title('Celerity $c/(\lambda_c \, \mathrm{d})$');
	case {'tlidi'}
		title('Celerity $c/(\lambda_c \, \mathrm{[T]})$');
	end
	colormap(parula(8))

	% cover
	splitfigure([2, 5], [fid+psldx-1, 8], fflag);
	cla();
	pfun(R,cva,smooth(o.cover,ns),0:0.05:1);
	view(0,90);
	ch=colorbar();
	axis xy;
	axis tight
	axis square
	xlabel(xlabel_str);
	title('Vegetation Cover $p_c$');
	ylabel('$\hspace{-0.8em}$Exogenous heterogeneity $\mathrm{cv}(e)$');
	caxis_.rietkerk = [0.2,0.85];
	caxis_.tlidi    = [0.2,0.65];
	caxis(caxis_.(base_str));
	switch (base_str)
	case {'rietkerk'}
		colormap(parula(13))	
	case{'tlidi'}
		colormap(parula(9))
	end


	% mean biomass (0)
	splitfigure([2, 5], [fid+psldx-1, 9], fflag);
	cla();
	level.tlidi       = 0:0.05:1;
	level.rietkerk    = 0:1:20;
	pfun(R,cva,smooth(o.b_mean,ns),level.(base_str));
	view(0,90);
	ch=colorbar();
	axis xy;
	axis tight
	axis square
	xlabel(xlabel_str);
	title('\hspace{1em}Mean biomass density $\bar b$');
	ylabel('$\hspace{-0.8em}$Exogenous heterogeneity $\mathrm{cv}(e)$');
	caxis_.rietkerk = [5,16];
	caxis_.tlidi = [0.15,0.55]
	nc.rietkerk = 11;
	nc.tlidi = 8;
	caxis(caxis_.(base_str));
	colormap(parula(nc.(base_str)))	

	if (pflag)
		pdfprint((fid+psldx-1)*10+1,sprintf('img/sw-%s-psl-%g-stripe-width-vs-R-and-cva.pdf',base_str,psl(psldx)),plotscale);
		pdfprint((fid+psldx-1)*10+2,sprintf('img/sw-%s-psl-%g-cvw-vs-R-and-cva.pdf',base_str,psl(psldx)),plotscale);
		pdfprint((fid+psldx-1)*10+3,sprintf('img/sw-%s-psl-%g-skw-vs-R-and-cva.pdf',base_str,psl(psldx)),plotscale);
		pdfprint((fid+psldx-1)*10+4,sprintf('img/sw-%s-psl-%g-bimodality-vs-R-and-cva.pdf',base_str,psl(psldx)),plotscale);
		pdfprint((fid+psldx-1)*10+5,sprintf('img/sw-%s-psl-%g-regularity-vs-R-and-cva.pdf',base_str,psl(psldx)),plotscale);
		pdfprint((fid+psldx-1)*10+6,sprintf('img/sw-%s-psl-%g-wavelength-vs-R-and-cva.pdf',base_str,psl(psldx)),plotscale);
		pdfprint((fid+psldx-1)*10+7,sprintf('img/sw-%s-psl-%g-celerity-vs-R-and-cva.pdf',base_str,psl(psldx)),plotscale);
		pdfprint((fid+psldx-1)*10+8,sprintf('img/sw-%s-psl-%g-cover-vs-R-and-cva.pdf',base_str,psl(psldx)),plotscale);
		pdfprint((fid+psldx-1)*10+9,sprintf('img/sw-%s-psl-%g-biomass-vs-R-and-cva.pdf',base_str,psl(psldx)),plotscale);
	end % pflag
	
	end % psldx
	
	% hellinger distance over all
	figure(fid+1);
	f_C = {'normal','no2','qgauss'};
	%s=sprintf('Noisy\rOscillator');
	dx_leg = 0.2;
	leg_C = {'Normal','Noisy   ','Oscillator','Q-Gauss'};
	cla;
	hd_a = [];
	col = [0.05,0.05,0.05;
	       meta.colour.blue;
	       meta.colour.red];
	p = 0.8;
	col = p*col + (1-p)*ones(3);
	for idx=1:length(f_C)
		hd_a(idx) = nanrms(abs(flat(o.hd.(f_C{idx}))));
		bar(idx,hd_a(idx),'facecolor',col(idx,:));
		hold on;
		y = hd_a(idx);
		text(idx-0.33,y+0.01,sprintf('%0.3f',y));
	end
	ylabel('$\hspace{-1.5em}$RMS Hellinger Distance');
	xlim([0.5,3.5]);
	ylim([0,1.2*max(hd_a)]);
	ylim([0,0.175]);
	set(gca,'xtick',[1,2-dx_leg,2+dx_leg,3],'xticklabel',leg_C,'xticklabelrot',45);
	%set(gca,'xtick',1:length(f_C),'xticklabel',leg_C,'xticklabelrot',45)
	axis square
	% expected random error
	if (0)
		hd_random = rms(sqrt(1./(8*o.n(:,:,pdx).*o.dxh(:,:,pdx)) + 1./(96*o.wsd(:,:,pdx).^2)),'all');
	end
	
	if (pflag)
		plotscale_ = 3.5;
		pdfprint(fid+1,sprintf('img/sw-%s-hellinger-distance-rms.pdf',o.model),plotscale_);
	end
	
	function x=smooth(x,nf)
		w=normpdf([-1.5:1.5:1.5],0,1);
		w=w'*w;
		w=w/sum(w,'all');

		for idx=1:nf
			x = [x(:,1),x,x(:,end)];
			x = [x(1,:);x;x(end,:)];
		if (1)
			x = conv2(x,w,'valid');
		else
			m = 3;
			x = ordfilt2(x,(3^2+1)/2,ones(m));
			x = x(2:end-1,2:end-1);
		end
		end
	end % smooth
	
end % sw_model_heterogeneous_plot

