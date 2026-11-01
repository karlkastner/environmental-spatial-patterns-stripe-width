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
%% basic parameter settings for Tlidi model runs
%
% function param = sw_tlidi_param(mu,cva,psl,T,L,dx)
function param = sw_tlidi_param(mu,cva,psl,T,L,dx)
	% model parameter	
	param = struct();

	% spatial extent of computational domain (determines spectral resolution)
	param.L    = L;

        % number of grid cells
	param.nx   = param.L(1)./dx;
    
	% simulated time interval
	param.T             = [0,T];

	tscale = 0.055;
	%xscale = 3.7;
	xscale = 2;

	% biomass equilibrium
	param.pmu.bc   = 1;
	param.pmu.dirc = -0.2*xscale;
	param.pmu.dirf =  0.8*xscale;
	% facilitation scale
	param.pmu.ef   = tscale*1;
	% biomass diffusion coefficient 
	param.pmu.ex = {tscale*0.01/xscale^2};
	% competition kernel extent
	param.pmu.lc   = 30*xscale;
	% facilitation kernel extent
	param.pmu.lf = 1*xscale;
	% mortality rate
	param.pmu.mu   = tscale*mu;
	% competition power
	param.pmu.xc   = 3; 
	% facilitation power
	param.pmu.xf   = 3; 

        % standard deviation per unit area of of the spatially heterogeneous parameter
        param.pss.ef   = cva*param.pmu.ef;
      
	if (cva ~= 0)
		% correlation length of the spatially heterogeneous parameter
	        param.psl.ef   = psl;
		% distribution of the spatially heterogeneous parameter        
		param.psdist.ef  = 'geometric-ornstein-uhlenbeck';
	end

	% time step of numerical integration
	dt = 0.1*dx;
	% data type during compution
	param.opt.compute_class = @single;
	% reload values of intermediate time steps
	% TODO this should not be a param, but argument to run
	param.opt.loadfinal     = false;
	% data type of output file
	param.opt.output.class  = @half;
	% time step for storing output
	param.opt.output.dt     = T/100;
	param.opt.rng           = 0;
	% activate time step adaptation
	param.opt.time_integration.adapt_step_length = 1;
	% initial time step
	param.opt.time_integration.dt_initial   = dt;
	param.opt.time_integration.dt_max       = inf; % dt;
	param.opt.time_integration.dt_min       = 0;
	param.opt.time_integration.abstol       = 1e-3;
	param.opt.time_integration.reltol       = 1e-2;
	param.opt.time_integration.dt_max_scale_up   = sqrt(2);
	param.opt.time_integration.dt_min_scale_down = 0;
	param.opt.time_integration.dt_safety_factor = 0.9;
	% integrate in time with trapezoidal scheme
	param.opt.time_integration.q = 0.5;
	% integrate in real space
	%param.opt.time_integration.scheme = 'step_advect_diffuse_implicit_q_fft';
	param.opt.time_integration.scheme = 'step_advect_diffuse_implicit_q';

	% random initial condition
	param.initial_condition           = 'obj.random_state()';

	param.opt.output.path_str = sprintf('mat/tlidi-L-%g-T-%g-dx-%g-sl-%g/',L(1),T,dx,psl);
end

