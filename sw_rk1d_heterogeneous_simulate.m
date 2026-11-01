% Sat 31 Aug 12:49:37 CEST 2024
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
%% function sw_rk1d_simulate(pflag)
%%
%% generate patterns with the RK model in heterogenous environment
%% vary water availability and the degree of heterogeneity between runs
%%
%% input:
%%
%% tabulate_summary : set to true to tabulate output files, in case of parallel
%% processing, run parallel processes with plot_summary = false, and after
%% all simulations have finished, once with plot_summary = true
%
%function swt = sw_rk1d_heterogeneous_simulate(fid,tabulate_summary,plot_individual,pflag)
function swt = sw_rk1d_heterogeneous_simulate(fid,tabulate_summary,plot_individual,pflag)
    if (nargin()<1 || isempty(fid))
        fid = 1;
    end
    if (nargin()<2 || isempty(tabulate_summary))
        tabulate_summary = true;
    end
    if (nargin()<3||isempty(plot_individual))
        plot_individual = true;
    end
    if (nargin()<4)
        pflag = false;
    end
    plotscale = 4.5;
    meta = sw_metadata();

    % use only 1 core, run in parallel by launching several matlab instances
    maxNumCompThreads(meta.maxNumCompThreads);
    
    % spatial extent
    L = 2^14;

    % spatial resolution
    dx       = 1;

    % duration of run
    T = 2*1e2*365;

    % time from which on pattern properties are averaged to estimate the
    % mean asymptotic properties
    T_avg = 0.9*T;
  
    % degree of exogeneous spatial variation (coefficient of variation)    
    % should not include 0, as homogeneous pattern structurally differ
    cva = (0.01:0.01:0.21); % 0.01
    %cva = 0.05;
    %cva = 0.08; 

    % correlation length scale of exogeneous heterogeneity
    psl = 1;
 
    % precipitation rate
    R = (0.70:0.05:1.3);
        
    % values for which the individual pattern is plot
    cvaplot= NaN; %cvaplot = 0.08;
    Rplot = NaN;  %Rplot  = 0.85;
 
    % bins width of stripe width histogram
    dxi  = 0.5*dx;
    % histogram bins (the largest bin is larger than the maximum stripe width
    % to prevent the fit of distributions with heavy tails)
    xi   = ((0.5*dxi):dxi:100)';
 
    % declare output variables
    sd_a   = cell(length(cva),length(R),length(psl));
    sg_a   = cell(length(cva),length(R),length(psl));
    stat_a   = cell(length(cva),length(R),length(psl));

    printf('Running %g simulations.\n',length(cva)*length(R));

    % for each correlation length
    for psldx=1:length(psl)
        for rdx=1:length(R)
            % for each degree of environemntal heterogeneity
            for cdx=1:length(cva)
                disp([psldx,cdx,rdx,psl(psldx),cva(cdx),R(rdx)]);

                % base model configuration
                param = sw_rk1d_param(R(rdx),cva(cdx),psl(psldx),T,L,dx);
        
                % create output folder
                mkdir(param.opt.output.path_str);
                    
                % inialize simulation
                rad = Rietkerk(param);
            
                % run the model with current parameter configuration
                [t,y,out_rad] = rad.run();

                % analyze the swd of the generated pattern
                if (~isempty(y))
                    [sd,sd_t,sg,sg_t,stat_t] = sw_model_analyze(rad,t,y,xi,T_avg);

                    % store result in ram
                    if (tabulate_summary)
                        sg.tab.filename = rad.filename;
                        sg.tab.n_step   = sum(out_rad.n_step);
                        sg.tab.runtime  = out_rad.runtime(end);
                        sd_a{cdx,rdx,psldx} = sd;
                        sg_a{cdx,rdx,psldx} = sg;
			stat_a{cdx,rdx,psldx} = stat_t;
                    end
            
                    % plot the swd for the current simulated pattern
                    if (isempty(fieldnames(sg.tab)))
                        disp('empty analysis')
                    else
                        if (  plot_individual ...
                           && any(   (abs(R(rdx)   - Rplot)   < sqrt(eps)) ...
                                 &   (abs(cva(cdx) - cvaplot) < sqrt(eps)) ...
                                 ) ...
                           )
                           sw_model_plot_individual( ...
                             fid ...
                            ,rad,t,y,sd,sd_t,sg, stat_t, R(rdx), cva(cdx) ...
                            ,pflag,T_avg,psl(psldx),plotscale ...
                            );
	                    fid = fid+1;
                            drawnow();
                            % avoid clogging with windows
                            %close all;
                        end % if plot_individual
                    end % if ~isempty sg.tab
                end % if ~isempty y
            end % for cdx
        end % for rdx
    end % for psldx
    
    % tabulate summary statistic
    if (tabulate_summary)
        sw_model_heterogeneous_tabulate(rad.opt.output.path_str ...
					, rad.opt.output.base_str ...
					, R, cva, psl, sd_a, sg_a ...
				        , stat_a ...
					);
    end % if tabulate_summary
end % sw_rk1d_heterogenous_simulate

