% Fri 27 Mar 14:02:34 CETh2026
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
function sw_experiment_correlated_heterogeneity(fid,pflag)
	if (nargin()<1)
		fid = 1;
	end
	if (nargin()<2)
		pflag = 0;
	end
	fflag = pflag;
	% plotscale
	plotscale = 5;
	meta = sw_metadata();

	init_plot();
	
	% spatial resolution
	dx = 0.1;
	% characteristic wavelength 
	%lc = 100;
	% correlation length
	Lc = 3;
	% spatial extent
	L=1e3*Lc;
	% number of grid cells
	n=round(L/dx);
	
	% number of repetitions for monte-carlo simulation
	m=1e3; % e5
	
	
	% x-axis coordinates
	x=(1:n)'*dx;
	% index
	id = (1:n)';
	
	
	% spatially (symetrically) correlated noise
	[ec,e0,rho,R,T,x_] = exp_spatial_generate(L,[n,m],Lc);
	
	% variance of the process
	s2    = mean(var(ec,[],2));
	% the mean is zero, so this should be constant
	s2(2) = mean(mean(ec.*ec));
	printf('Variance of correllated process: %g %g\n',s2);
	scale          = 2*Lc;
	scale_discrete = dx*sum_ii(rho,n)/n;
	printf('Ratio of asymptotic variance %g\n',scale);
	printf('Ratio of asymptotic variance, discrete %g\n',scale_discrete);
	
	% allocate memnory
	W2 = zeros(n,2);
	W2_analytic = zeros(n,3);
	
	% brownian motion with uncorrelated noise
	W0 = sqrt(dx)*cumsum(e0);
	% variance of the simulated processes
	% (increases with distance)
	W2(:,1) = mean(W0.*W0,2);
	
	% brownian motion with correlated noise
	Wc      = sqrt(dx)*cumsum(ec);
	% variance of the simulated processes
	W2(:,2) = mean(Wc.*Wc,2);
	
	% scaled correlated process
	%Wcs = scale*Wc;
	% variance of the simulated processes
	%W2(:,3) = mean(Wcs.^2,2);
	
	% analytic solution - discrete
	%y0=(1/(1-4*r^2))^(1/2)
	ir_0 = (1+rho^2)/(1-rho^2);
	ir_sum = ir_0*(1+rho)/(1-rho);
	ir2_sum = ir_0.^2*(1+rho^2)/(1-rho^2);
	% analytic variance, uncorrelated process
	W2_analytic(:,1) = dx*(1:n)';
	% analytic variance, discrete correlated process
	W2_analytic(:,2) = cvec(dx^2*sum_ii(rho,id));
	% analytic variance, continuous correlated process
	W2_analytic(:,3) = 2*Lc*(x - Lc*(1-exp(-x/Lc)));
	
	% plot ratio of variance between scaled correlated and uncorrelated process
	splitfigure([2,2],[fid,1],fflag);
	cla()
	plot(x/Lc,W2_analytic(:,2)./(scale*W2_analytic(:,1)),'linewidth',meta.linewidth,'color',meta.colour.blue)
	% ylabel('$E[W_{Lc}^2]/\left(2 L_c E[W_0^2]\right)$')
	ylabel('$\hspace{-2em}\mathrm{E}\!\left[(\varphi_{e,L_c})^2\right]\!/\!\left(2 \, L_c  \mathrm{E}\!\left[(\varphi_{e,0})^2\right]\right)$')
	if (~pflag)
		hold on;
		plot(x/Lc,[W2_analytic(:,3)./(scale*W2_analytic(:,1)),W2(:,2)./(scale*W2(:,1))]);
		legend('analytic discrete','analytic continuous','simulated')
	end
	hline(1,'color','k','handlevisibility','off');
	xlim([0,10])
	xlabel('Distance $x / L_c$');
	ylim([0,1.05])
	axis square
	
	
	% plot ratio of variance between unscaled correlated and uncorrelated process
	splitfigure([2,2],[fid,2],fflag);
	cla()
	plot(x/Lc,W2(:,2)./W2(:,1),'linewidth',meta.linewidth,'color',meta.colour.blue)
	ylabel('$E[W_{L_c}^2]/E[W_0^2]$')
	if (~pflag)
		hold on;
		plot(x/Lc,W2_analytic(:,2:3)./W2_analytic(:,1));
		legend('simulated','analytic discrete','analytic continuous')
	end
	hline(scale,'color','k','handlevisibility','off');
	xlim([0,10])
	xlabel('Distance $x / L_c$');
	ylim([0,scale*1.05])
	axis square
	
	% plot variance directly, i.e. not as a ratio
	splitfigure([2,2],[fid,3],fflag);
	cla();
	plot(x/Lc,[scale*W2(:,1),W2(:,2)]);
	xlabel('Distance $x / L_c$');
	ylabel('Variance $E[W^2]$');
	xlim([0,10])
	legend('uncorrelated','correlated, scaled');
	
	% plot standard deviation directly, i.e. not as a ratio
	splitfigure([2,2],[fid,4],fflag);
	cla();
	plot(x/Lc,sqrt([scale*W2(:,1),W2(:,2)]));
	xlim([0,10])
	xlabel('Distance $x / L_c$');
	ylabel('Standard deviation $\sqrt{E[W^2]}$');
	legend('uncorrelated','correlated, scaled');
	
	if (pflag)
		pdfprint(10*fid+1,'img/correlated-noise-variance-ratio.pdf',plotscale);
	end
end

