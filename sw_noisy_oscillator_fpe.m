% Tue  2 Jun 09:44:32 CEST 2026
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
%% sw_noisy_oscillator_fpe(pflag)
%%
%% plot the stripe width distribution of the noisy oscillator, determined
%% by solving the Fokker-Planck equation
function sw_noisy_oscillator_fpe(fid,pflag)
	if (nargin()<1)
		fid = 1;
	end
	if (nargin()<2)
		pflag = 0;
	end
	fflag = pflag;
	plotscale = 4;
	init_plot();
	meta = sw_metadata();
	set(groot,'defaultaxescolororder',[meta.colour.blue; meta.colour.red])
	n = 12;
	x=(0:n-1)'/(n-1);
	hsv=([0.6*ones(n,1), 1-x, x.^0.5]);
	cmap_surface=flipud(hsv2rgb(hsv));
	
	% charactersitic wavelength
	lc     = 100;
	% fraction of ground covered by vegetation
	pc     = 0.3;
	% corresponding phase 
	phi_c  = 2*pi*pc;
	% characteristic stripe wdith
	wc     = pc*lc;
	% spatial extent
	Lx     = 3*wc;
	% standard deviation of exogeneout heterogeneity
	s      = 2.6;
	% hysteresis for thresholding biomass
	x_hyst = 0.26*wc;
	% plot extent
	xlim_  = [0,0.6*lc];
	% set test=true for comparison with monte-carlo simulation
	test   = false;
	% discretized phase space
	nphi = 1e3;
	x   = innerspace(0,Lx,1e3)';
	phi = linspace(-phi_c,3*phi_c,nphi);
	philim_ = [-0.067*2*pi,1.001*phi_c];

	meta = sw_metadata();
	
	% allocate memory
	clear ax1 ax2 cb1 cb2 hs1 hs2
	% first loop run without short stripes then with short stripes 
	for pdx=1:2
		
		switch (pdx)
		case {1}
			% swd of the noisy osicallator, short stripes ignored
			%[f,F,out] = stripe_width_1_pdf(x/lc,wc/lc,s/sqrt(lc));
			[f,F,out] = stripe_width_1_pdf_fpe(x/lc,phi,wc/lc,s/sqrt(lc));
			p = out.p';
			if (test)
				% analytic solution
				f(:,2) = stripe_width_1_pdf(x,wc,s*sqrt(wc));
			end
		case {2}
			% swd of the noisy osicallator short stripes considered
			[f,F,out] = stripe_width_2_pdf_fpe(x/lc,wc/lc,s/sqrt(lc),x_hyst/lc,Lx/lc,meta.n_fpe);
			p = out.p;
			if (test)
				% TODO use simulation
				%n  = 1e3*[1,1];
				%x_hyst_ = wc;
				%f_ = stripe_width_2_pdf_fpe(x,wc/lc,s/sqrt(lc),x_hyst_/lc,Lx/lc,n_fit);
			end
		case {3}
			error('here')
		end
		
		% plot probability p of process to have the phase value phi at distance x
		splitfigure([2,2],[fid+pdx-1,1],fflag);
		cla();
		hs1(pdx) = surface(lc*out.x,out.phi,p,'edgecolor','none');
		ax1(pdx) = gca;
		axis xy;
		ylim(philim_)
		xlabel('Distance x = Stripe width w');
		ylabel('Phase $\varphi$');
		cb1(pdx)=colorbar();
		title(cb1(pdx),'$p$','interpreter','latex');
		colormap(cmap_surface);
		caxis([0,2])
		axis square
		cb1(pdx).Ticks = (0:0.5:2);
		xlim(xlim_);
		
		% normalized plot
		splitfigure([2,2],[fid+pdx-1,2],fflag);
		cla();
		% p is integrated along phi but differentiated along x
		hs2(pdx) = surface(out.x,out.phi/(2*pi),p*(2*pi/lc^0),'edgecolor','none');
		ax2(pdx) = gca;
		axis xy
		ylim(philim_/(2*pi));
		xlabel('Stripe width $w/\lambda_c$');
		ylabel('Phase $\varphi/(2 \pi)$');
		cb2(pdx)=colorbar();
		title(cb2(pdx),'$\frac{p}{2 \pi}$','interpreter','latex');
		colormap(cmap_surface);
		caxis([0,1.5]*1e3/1e2);
		axis square
		cb2(pdx).Ticks = (0:0.5:2)*1e3/1e2;
		xlim(xlim_/lc);

		% survivial function and probability density
		splitfigure([2,2],[fid+pdx-1,3],fflag);
		yyaxis left
		cla();
		plot(x,1-F,'linewidth',meta.linewidth);
		ylabel('Survival $1-F=\int{p}\,\mathrm{d}\,\varphi$')
		ylim([-0.0,1.05]);
		yyaxis right;
		cla();
		plot(x,f/lc,'linewidth',meta.linewidth);
		ylabel('Density $f = \frac{\partial{F}}{\partial{x}}$')
		%ylabel('Density $-\lambda_c \, \frac{\partial F}{\partial x}$')
		xlabel('Stripe width $w$');
		ylim([-0.0,1.05]*max(f(:,1))/lc);
		axis square
		xlim(xlim_);
		
		% normalized
		splitfigure([2,2],[fid+pdx-1,4],fflag);
		yyaxis left
		cla
		plot(x/lc,1-F,'linewidth',meta.linewidth);
		ylabel('Survival $1-F = \int p \, \mathrm{d} \, \varphi$')
		ylim([-0.0,1.05]);
		yyaxis right;
		cla();
		plot(x/lc,f(:,1),'linewidth',meta.linewidth);
		hold on
		if (1 == pdx)
		if(test)
			plot(x/lc,f(:,2)*lc,'--b','linewidth',meta.linewidth)
			plot(x/lc,f_*lc,':g','linewidth',meta.linewidth)
		end
		else
			plot(out.x,[out.f12],'--','color',meta.colour.red,'linewidth',meta.linewidth)
		end
		ylabel('Probability density $\lambda_c \, \frac{\partial{F}}{\partial{x}}$')
		xlabel('Stripe width $w/\lambda_c$');
		ylim([0,3.75]);
		axis square
		xlim(xlim_/lc);
		
		if (pflag)
			% recent matlab version cannot plot the content as raster and axis as vector
			% so plot them separately and recombine them manually
			pdfprint(10*(fid+pdx-1)+2,sprintf('img/stripe-width-fokker-planck-%g.pdf',pdx),plotscale);
			hs2(pdx).Visible = 'off';
			ax2(pdx).Box = 'on';
			pdfprint(10*(fid+pdx-1)+2,sprintf('img/stripe-width-fokker-planck-%g-axes.pdf',pdx),plotscale);
			hs2(pdx).Visible = 'on';
			cb2(pdx).Visible = 'off';
			ax2(pdx).Box = 'on';
			ax2(pdx).XTick = [];
			ax2(pdx).YTick = [];
			ax2(pdx).XLabel.String = '';
			ax2(pdx).YLabel.String = '';
			pdfprint(10*(fid+pdx-1)+2,sprintf('img/stripe-width-fokker-planck-%g-content.png',pdx),plotscale,[],'png');
		
			pdfprint(10*(fid+pdx-1)+4,sprintf('img/stripe-width-distribution-with-survival-function-%g.pdf',pdx),plotscale);
		end % if pflag
	end % for pdx
end % sw_noisy_oscillator_fpe

