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
%% function [sd_a,sd_t,sg_a,sg_t,stat_t] = sw_model_analyze(rad,t,b,xi,T_avg,meta)
%%
%% input:
%%
%%	rad : RAD model object
%%	t   : time of output
%%	yy  : rad state at time t (biomass and water)
%%	xi  : histogram bins, ksdensity ordinate
%%     T_avg : time to average samples over
%%	meta : processing settings
%%
%% output:
%%
%% sd_a   : spectral density, averaged over the last T_avg time steps
%% st_t   : spectral density, every time step
%% sg_a   : stripe width distribution, averaged over the last T_avg time steps
%% sg_t   : stripe width distribution for every time step
%% stat_t : statistics for each time step in b
function [sd_a,sd_t,sg_a,sg_t,stat_t] = sw_model_analyze(rad,t,yy,xi,T_avg,meta)
	if (nargin()<6)
		meta = sw_metadata();
	end
	tdx = (t>T_avg);
	yy = double(yy);

	bb = rad.extract1(yy);
	
	% initialize output objects
	sd_a   = Spectral_Density();
	sd_t   = Spectral_Density();
	sg_a   = Stripes_and_Gaps();
	sg_t   = Stripes_and_Gaps();
	stat_t = struct();
	stat_t.b    = struct();
	stat_t.w    = struct();

	if (isempty(bb))
		disp('Empty imput')
	else
	sg_a.opt.n_fpe   = meta.n_fpe;

	% load preprocessed file name, if existing
	name = rad.filename();
	name = [name(1:end-4),'-analyzed.mat'];
	if (exist(name,'file'))
		try
			load(name,'sd_a','sd_t','sg_a','sg_t','stat_t');
		catch e
			e
		end
	else

	% since the advection-diffusion discretization matrices are not stored, they have to be reinitialized
	% matrices are required for the celerity computation
	rad.init();

	fx = rad.fourier_axis;
	nt = length(t);
	
	% reserve filename for parallel processing by writing an empty file
	save(name,'sd_a','sd_t','sg_a','sg_t','stat_t');
	
	
	% for each time step in the output record
	for idx=1:length(t)
		b = cvec(bb(idx,:));
		
		sd_t(idx,1) = Spectral_Density();
		sd_t(idx,1).init(b,rad.L(1));
		sd_t(idx,1).estimate();
		
		sg_t(idx,1) = Stripes_and_Gaps();
		sg_t(idx,1).init(b,rad.L(1));
		sg_t(idx,1).opt.xe = inner2outer(xi);
		sg_t(idx,1).opt.bandwidth = meta.bandwidth;
		sg_t(idx,1).analyze();
		% sg_t(idx).fit_parametric_densities();
	
		stat_t.celerity(idx,1) = rad.celerity(t,yy(idx,:));
	end % for idx

	% spectral density estimate (only for final time step)
	%sd_t.init(b,rad.L(1));
	%sd_t.estimate();
	% init the last pattern, just to set length
	sg_a.L  = rad.L(1);
	sg_a.nx = rad.nx(1);
	
	% store final start and end points
	sg_a.tab.xse = sg_t(end).tab.xse;

	% extract values into stat_t
	stat_t.b.mean = arrayfun(@(x) x.tab.b.mean, sg_t);
	stat_t.b.sd   = arrayfun(@(x) x.tab.b.sd, sg_t);
	stat_t.b.cv   = arrayfun(@(x) x.tab.b.cv, sg_t);
	stat_t.b.min  = arrayfun(@(x) x.tab.b.min, sg_t);
	stat_t.b.max  = arrayfun(@(x) x.tab.b.max, sg_t);
	stat_t.cover  = arrayfun(@(x) x.tab.cover, sg_t);

	stat_t.S.hat     = cell2mat(arrayfun(@(x) x.S.hat,sd_t,'uniformoutput',false)')';
	stat_t.S.bar     = cell2mat(arrayfun(@(x) x.S.bar,sd_t,'uniformoutput',false)')';
	stat_t.Sc.hat    = arrayfun(@(x) x.Sc.hat,sd_t);
	stat_t.Sc.bar    = arrayfun(@(x) x.Sc.bar,sd_t);
	stat_t.fc.hat    = arrayfun(@(x) x.fc.hat,sd_t);
	stat_t.fc.bar    = arrayfun(@(x) x.fc.bar,sd_t);

	% average over last T_avg time steps
	sg_a.tab.b.thresh = sg_t(end).tab.b.thresh;
	sg_a.tab.b.mean   = mean(stat_t.b.mean(tdx)); 
	sg_a.tab.b.sd     = mean(stat_t.b.sd(tdx)); 
	sg_a.tab.b.min    = mean(stat_t.b.min(tdx)); 
	sg_a.tab.b.max    = mean(stat_t.b.max(tdx)); 
	sg_a.tab.b.cv     = sg_a.tab.b.sd./sg_a.tab.b.mean;
	sg_a.tab.cover    = mean(stat_t.cover(tdx)); 
        sg_a.tab.celerity  = mean(stat_t.celerity(tdx));

	% TODO extract from stat_a.Sc from average density
	tdx_ = tdx;
	sd_a.L = rad.L;
	sd_a.nx = rad.nx;
	% average periodogram and density
	sd_a.S.hat     = mean(stat_t.S.hat(tdx_,:),1);
	sd_a.S.bar     = mean(stat_t.S.bar(tdx_,:),1);
	sd_a.Sc.hat    = mean(stat_t.Sc.hat(tdx_));
	sd_a.Sc.bar    = mean(stat_t.Sc.bar(tdx_));
	sd_a.fc.hat    = mean(stat_t.fc.hat(tdx_));
	sd_a.fc.bar    = mean(stat_t.fc.bar(tdx_));

	% average stripe width (jdx=1) and gap width (jdx=2) distribution
	for jdx=1:2
		sg_a.tab.w(jdx).xe   = sg_t(end).tab.w(jdx).xe;
		sg_a.tab.w(jdx).xc   = sg_t(end).tab.w(jdx).xc;
		stat_t.w(jdx).moment = cell2mat(arrayfun(@(x) x.tab.w(jdx).moment,sg_t,'uniformoutput',false));
                stat_t.w(jdx).mean   = arrayfun(@(x) x.tab.w(jdx).mean,sg_t);
                stat_t.w(jdx).sd     = arrayfun(@(x) x.tab.w(jdx).sd,sg_t);
                stat_t.w(jdx).sk     = arrayfun(@(x) x.tab.w(jdx).sk,sg_t);
                stat_t.w(jdx).ku     = arrayfun(@(x) x.tab.w(jdx).ku,sg_t);
        	stat_t.w(jdx).pdf.hist   = cell2mat(arrayfun(@(x) x.tab.w(jdx).pdf.hist,sg_t,'uniformoutput',false)');
        	stat_t.w(jdx).pdf.ksdens = cell2mat(arrayfun(@(x) x.tab.w(jdx).pdf.ksdens,sg_t,'uniformoutput',false)');

		sg_a.tab.w(jdx).moment     = mean(stat_t.w(jdx).moment(tdx,:),1);
                sg_a.tab.w(jdx).mean       = mean(stat_t.w(jdx).mean(tdx));
		sg_a.tab.w(jdx).sd         = sqrt(mean(sg_a.tab.w(jdx).moment(2)));
		sg_a.tab.w(jdx).sk         = sg_a.tab.w(jdx).moment(3)./sg_a.tab.w(jdx).sd.^3;
		sg_a.tab.w(jdx).ku         = sg_a.tab.w(jdx).moment(4)./sg_a.tab.w(jdx).sd.^4;
        	sg_a.tab.w(jdx).pdf.hist   = mean(stat_t.w(jdx).pdf.hist(:,tdx),2);
        	sg_a.tab.w(jdx).pdf.ksdens = mean(stat_t.w(jdx).pdf.ksdens(:,tdx),2);
	end % for jdx

	% fit parametric densities to average stripe width density
	sg_a.fit_parametric_densities();

	% bimodality
	for idx=1:2
		[sg_a.tab.w(idx).bimodality.d, sg_a.tab.w(idx).bimodality.out] = bimodality(sg_a.tab.w(idx).pdf.ksdens,sg_a.tab.w(idx).xc,'absolute');
	end

	% write results to disc
	% write first to temporary file to avoid accidental load of files which are currently written
	tmpname_str = [tempname(dirname(name)),'.mat'];
	save(tmpname_str,'sg_a','sd_a','sd_t','sg_t','stat_t');
	movefile(tmpname_str,name);

    end % if ~exist analysis output file
    end % ~isempty(bb)
end % sw_model_analyze

