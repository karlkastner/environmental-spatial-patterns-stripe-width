% Thu 11 Jun 19:03:05 CEST 2026
% Karl Kastner (Berlin)
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
%% tabulate results of individual model runs
%% and determine parameter sensitivity by linear regression
%%
%% input :
%%	base_str : model that generated patterns "rietkerk" or "tlidi"
%%	R        : n1 x 1 water availability (proxy for Tlidi)
%%	cva	 : n2 x 1 degree of exogenous spatial heterogeneity
%%	psl	 : correlation length of heterogeneoity
%%	sd_a     : n1 x n2 spectral density of respetive model configuration
%%	sg_a     : n1 x n2 stripe and gap width statistics of respetive model configuration
function sw_model_tabulate(path_str,base_str,R,cva,psl,sd_a,sg_a,stat_a) 
	o       = struct();
	o.model = base_str;
	
	switch(base_str)
	case {'rietkerk'}
		o.R = R;
	case {'tlidi'}
	    	% proxy value water availability (critical biomass)
		mu = R;
		ef = 1;
	        water_availability  = 1-mu/ef;
		o.water_availability = water_availability;
	otherwise
		error('here');
	end
	o.cva = cva;
	o.psl = psl;
	
	% failsave, select only model runs that completed and were analyzed
	fdx=((cellfun(@(x) isfield(x.tab,'w'),sg_a)));
	
	o.w                     = struct();
	o.b_mean                = NaN(size(sg_a));
	o.cover                 = NaN(size(sg_a));
	o.fc.bar                = NaN(size(sg_a));
	o.Sc.bar                = NaN(size(sg_a));
	o.celerity              = NaN(size(sg_a));
	o.n_step                = NaN(size(sg_a));
	% for stripes (jdx=1) and gaps (jdx=2)
	for jdx=1:2
		o.w(jdx).mean      = NaN(size(sg_a));
		o.w(jdx).sd        = NaN(size(sg_a));
		o.w(jdx).sk        = NaN(size(sg_a));
		o.w(jdx).ku        = NaN(size(sg_a));
		o.w(jdx).mode      = NaN(size(sg_a));
		o.w(jdx).bimodality       = NaN(size(sg_a));
		o.w(jdx).bimodality_rel   = NaN(size(sg_a));
		o.w(jdx).bimodality_c     = NaN(size(sg_a));
		o.w(jdx).bimodality_small = NaN(size(sg_a));
		o.w(jdx).bimodality_large = NaN(size(sg_a));
		o.w(jdx).mean(fdx) = cellfun(@(x) x.tab.w(jdx).mean,sg_a(fdx));
		o.w(jdx).sd(fdx)   = cellfun(@(x) x.tab.w(jdx).sd,sg_a(fdx));
		o.w(jdx).sk(fdx)   = cellfun(@(x) x.tab.w(jdx).sk,sg_a(fdx));
		o.w(jdx).ku(fdx)   = cellfun(@(x) x.tab.w(jdx).ku,sg_a(fdx));
		% primary mode
		o.w(jdx).mode(fdx)             = cellfun(@(x) x.tab.w(jdx).bimodality.out.x.r,sg_a(fdx));
		o.w(jdx).bimodality(fdx)       = cellfun(@(x) x.tab.w(jdx).bimodality.d,sg_a(fdx));
		o.w(jdx).bimodality_c(fdx)     = cellfun(@(x) x.tab.w(jdx).bimodality.out.y.c,sg_a(fdx));
		o.w(jdx).bimodality_small(fdx) = cellfun(@(x) min(x.tab.w(jdx).bimodality.out.y.l,x.tab.w(jdx).bimodality.out.y.r),sg_a(fdx));
		o.w(jdx).bimodality_large(fdx) = cellfun(@(x) max(x.tab.w(jdx).bimodality.out.y.l,x.tab.w(jdx).bimodality.out.y.r),sg_a(fdx));
		o.w(jdx).bimodality_rel(fdx)   = 1 - o.w(jdx).bimodality_c(fdx)./o.w(jdx).bimodality_small(fdx);
		% the estimated density has local maxima (modes) due to the random variation,
		% small modes are therefore discarded here
		o.w(jdx).bimodality_rel(o.w(jdx).bimodality_small < 0.05*o.w(jdx).bimodality_large)  = 0;
	end
	%o.w(1).bimodality_rel(o.w(1).bimodality_small == 0) = 0;
	o.b_mean(fdx)                = cellfun(@(x) x.b.mean(end), stat_a(fdx));
	o.celerity(fdx)              = cellfun(@(x) x.tab.celerity,sg_a(fdx));
	o.cover(fdx)                 = cellfun(@(x) x.tab.cover, sg_a(fdx));
	o.fc.bar(fdx)                = cellfun(@(x) x.fc.bar,sd_a(fdx));
	o.n_step(fdx)                = cellfun(@(x) x.tab.n_step,sg_a(fdx));
	o.runtime(fdx)               = cellfun(@(x) x.tab.runtime,sg_a(fdx));
	o.Sc.bar(fdx)                = cellfun(@(x) x.Sc.bar,sd_a(fdx));
	o.hd = struct();
	f_C = fieldnames(sg_a{1}.tab.w(1).hd);
	for idx=1:length(f_C)
		o.hd.(f_C{idx})      = NaN(size(sg_a));
		o.hd.(f_C{idx})(fdx) = cellfun(@(x) x(1).tab.w(1).hd.(f_C{idx}),sg_a(fdx));
	end
	save([path_str,filesep(),'summary.mat'],'o');

	pdf  = struct();
	pdf.hist = cell(size(sg_a));
	pdf.ksdens = cell(size(sg_a));
	pdf.hist(fdx)   = cellfun(@(x) x.tab.w(1).pdf.hist,sg_a(fdx),'uniformoutput',false);
	pdf.ksdens(fdx) = cellfun(@(x) x.tab.w(1).pdf.ksdens,sg_a(fdx),'uniformoutput',false);
	save([path_str,filesep(),'pdf.mat'],'pdf');

	% expand dimensions for tabulating
	switch(base_str)
	case {'rietkerk'}
		RR = repmat(rvec(R),[length(cva),1,length(psl)]);
		n2 = length(R);
	case {'tlidi'}
		RR = repmat(rvec(water_availability),[length(cva),1,length(psl)]);
		n2 = length(water_availability);
	otherwise
		error('here');
	end
	cc = repmat(cvec(cva),[1,n2,length(psl)]);
	
	% write table with control and response variables to file
	otab = table;
	switch(base_str)
	case {'tlidi'}
		otab.water_availability = RR(:);
	case {'rietkerk'}
		otab.precipitation      = RR(:);
	otherwise
		error('here');
	end
	otab.cva           = cc(:);
	otab.fc            = o.fc.bar(:);
	otab.regularity    = flat((o.Sc.bar.*o.fc.bar));
	otab.celerity      = o.celerity(:);
	otab.sw_mean       = o.w(1).mean(:);
	otab.sw_sd         = o.w(1).sd(:);
	otab.sw_sk         = o.w(1).sk(:);
	otab.sw_bimodality = o.w(1).bimodality(:); 
	writetable(otab,['output/',base_str,'-summary.csv']);
	
	pp = repmat(cvec(cva),[1,n2,length(psl)]);
	pp = shiftdim(psl,-1);
	pp = repmat(pp,[length(cva),n2,1]);
	C  = {RR,cc}; %,log(pp)};
	sign_c        = median(sign(o.celerity),'all');
	
	% tabulate linear sensitivity
	lin_tab  = table();
	n        = numel(RR);
	mean_R   = mean(R);
	mean_cva = mean(cva);
	Lhs_lin = [ones(n,1),flat((RR - mean(R))/mean_R),flat((cc-mean_cva)/mean_cva)];
	rhs_lin = [  flat(o.b_mean) ...
		   , flat(o.cover) ...
		   , flat((1./o.fc.bar)) ...
		   , flat((o.Sc.bar.*o.fc.bar)) ...
		   , flat(sign_c*o.celerity.*o.fc.bar) ...
		   , flat(o.w(1).mean.*o.fc.bar) ...
		   , flat(o.w(1).sd./o.w(1).mean) ...
		   , flat(o.w(1).sk) ...
		   , flat(o.w(1).bimodality_rel) ... 
		];
	lin_tab{:,:} = Lhs_lin \ rhs_lin;
	lin_tab.Properties.VariableNames = { 'Biomass' ...
					   , 'Cover' ...
					   , 'Wavelength' ...
					   , 'Regularity' ...
					   , 'Celerity' ...
					   , 'sw_mean' ...
					   , 'sw_cv' ...
					   , 'sw_sk' ...
					   , 'sw_bimodality' ...
					   };
	log_tab.Properties.RowNames = {'mean','R-sensitivity','cva-sensitivity'};
	lin_tab{:,:}=roundn(lin_tab{:,:},2);
	disp(lin_tab)
	writetable(lin_tab, ['output/',base_str,'-lin-sensitivity.csv'] ...
			  , 'WriteRowNames',true);
	
	% tabulate geometry sensitivity
	log_tab = table();
	Lhs_log = [ones(n,1), flat(log(RR)), flat(log(cc))];
	rhs_log = [  flat(log(o.b_mean)) ...
		   , flat(log(o.cover)) ...
	           , flat(-log(o.fc.bar)) ...
		   , flat(log(o.Sc.bar.*o.fc.bar)) ...
		   , flat(log(sign_c*o.celerity.*o.fc.bar)) ...
		   , flat(log(o.w(1).mean.*o.fc.bar)) ...
		   , flat(log(o.w(1).sd./o.w(1).mean)) ...
		  ];
	log_tab{:,:} = Lhs_log \ rhs_log;
	% geometric mean
	log_tab{1,:} = exp(mean(rhs_log));
	log_tab.Properties.VariableNames = { 'Biomass' ...
					   , 'Cover' ...
					   , 'Wavelength' ...
					   , 'Regularity' ...
					   , 'Celerity' ...
					   , 'sw_mean' ...
					   , 'sw_cv' ...
					   };
	log_tab.Properties.RowNames = {'mean','R-sensitivity','cva-sensitivity'};
	log_tab{:,:}=roundn(log_tab{:,:},2);
	disp(log_tab)
	writetable(log_tab, ['output/',base_str,'-log-sensitivity.csv'] ...
			   , 'WriteRowNames',true);
end % sw_model_heterogeneous_tabulate

