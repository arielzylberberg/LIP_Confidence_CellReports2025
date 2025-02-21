classdef VarCE < handle
    properties
        x % data from each trial
        f % neuron label
        s % times of spikes ('hist' format)
        cn % column labels for x
        
        srt
        
        times
        s_time % time corresponding to s
        time_win
        ntimes
        
        ntr
        nneurons
        ncoh
        ndir
        ntrg
        ntrialspercond
        meancountpercond
        proptrialspercond
        varcountpercond
        uni_conditions
        conditions
        conditions_columns
        
        uni_coh
        uni_dir
        uni_trg
        uni_neuron
        idx_times
        
        spk
        spk_nomean
        
        fi
        varce
        corce
        covce
        fano
        
        
    end
    
    methods
        % constructor
        function obj = VarCE(x,f,s,cn,s_time,times,time_win)
            
            % sanity check
            if size(x,2)~=length(fields(cn))
                error('sizes do not match')
            end
            obj.x = x; % trial data
            obj.f = f; % neuron number/id
            obj.s = s; % spike times
            obj.cn = cn; % column names (from x)
            
            obj.s_time   = s_time;
            obj.times   = times;
            obj.time_win = time_win;
            obj.ntimes  = length(obj.times);
            obj.ntr     = length(obj.f);
            
            obj.uni_coh    = unique(obj.x(:,cn.R_COH));
            obj.uni_dir    = unique(obj.x(:,cn.R_DIR));
            % obj.uni_trg    = unique(obj.x(:,cn.R_TRG));
            obj.uni_neuron = unique(obj.f);
            
            obj.nneurons = length(obj.uni_neuron);
            obj.ncoh     = length(obj.uni_coh);
            % obj.ntrg     = length(obj.uni_trg);
            obj.ndir     = length(obj.uni_dir);
            
        end
        
        function count_spikes(obj)
            
            obj.spk = zeros(obj.ntr,obj.ntimes);
            for j = 1:obj.ntimes
                inds = obj.s_time>(obj.times(j)-obj.time_win/2) & ...
                    obj.s_time<=(obj.times(j)+obj.time_win/2);
                obj.spk(:,j) = nansum(obj.s(:,inds),2);
            end
            
        end
        
        function remove_mean_count_from_each_condition(obj)
            
            cn = obj.cn;
            
            [~,~,neuron] = unique(obj.f);
            conditions = [neuron obj.x(:,[cn.R_DIR cn.R_COH cn.R_TRG])];
            data = obj.spk;
            
            % just substract the mean count, as in the paper
            [Mean,Stdev,uni_conditions,tr_per_cond,idx_cond] = ...
                average_per_condition(data, conditions);
            obj.spk_nomean      = nan(size(obj.spk));
            inds = not(isnan(idx_cond));
            obj.spk_nomean(inds,:) = data(inds,:) - Mean(idx_cond(inds),:);
            
            
            
            obj.ntrialspercond      = tr_per_cond;
            obj.meancountpercond    = Mean;
            obj.conditions          = conditions;
            obj.uni_conditions       = uni_conditions;
            obj.conditions_columns  = {'NEURON','DIR','COH','TRG'};
            obj.proptrialspercond   = tr_per_cond/sum(tr_per_cond);
            obj.varcountpercond     = Stdev.^2;
            
        end
        
        
        function calc_fi(obj,varargin)
            boolCalcFi = true;
            for i=1:length(varargin)
                if isequal(varargin{i},'fi')
                    obj.fi = varargin{i+1};
                    boolCalcFi = false;
                end
            end
            
            if boolCalcFi
                cn = obj.cn;
                forficalc = nan(obj.nneurons,obj.ntimes);
                
                for i=1:obj.nneurons
                    inds = ismember(obj.f,obj.uni_neuron(i));% & obj.idx_times==1;
                    
                    m = nanmean(obj.spk(inds,:));
                    %                 m = obj.meancountpercond(i,:);
                    forficalc(i,:) = nanvar(obj.spk_nomean(inds,:))./m; %no se si esta bien; hay que usar spk??
                    
                end
                obj.fi = nanmin(forficalc,[],2);
            end
        end
        
        function calc_spike_res_autocorr(obj)
            
            
            res_autocorr = nan(obj.ntimes);
            for i=1:obj.ntimes
                for j=1:obj.ntimes
                    inds = not(isnan(obj.spk_nomean(:,i))) & not(isnan(obj.spk_nomean(:,j)));
                    res_autocorr(i,j) = corr(obj.spk_nomean(inds,i),obj.spk_nomean(inds,j));
                end
            end
            
            if 1
                figure
                imagesc(obj.times,obj.times,res_autocorr,[0,1])
                axis square
            end
        end
        
        function calc_varCE(obj,plot_flag)
            
            % varce
            aux_fi      = obj.fi(obj.uni_conditions(:,1));
            temp        = bsxfun(@times,obj.meancountpercond,obj.proptrialspercond.*aux_fi);
            
            obj.varce   = nanvar(obj.spk_nomean) - nansum(temp);
            
            % fano
            temp = bsxfun(@times,obj.meancountpercond,obj.proptrialspercond);
            obj.fano = nanvar(obj.spk_nomean)./nansum(temp);
            
            if plot_flag
                figure
                subplot(2,1,1)
                plot(obj.times,obj.varce)
                xlabel('time')
                ylabel('VarCE')
                
                subplot(2,1,2)
                plot(obj.times,obj.fano)
                xlabel('time')
                ylabel('FanoFactor')
            end
            
        end
        
        function search_fi_for_corCE(obj)
            
            fi = obj.fi;
            p = 1;
            plot_flag = 0;
            while p~=0 % not positive defined
                disp(['fi:',num2str(fi)]);
                obj.calc_fi('fi',fi);
                obj.calc_varCE(plot_flag);
                obj.calc_covCE(plot_flag);
%                 [~,p] = chol(obj.corce);
                [~,p] = chol(obj.covce);
                fi = fi - 0.005;
            end
        end
        
        function calc_covCE(obj,plot_flag)
            % no se si esto esta bien... parece mucho mas simple que la
            % explicacion
            covmat = nancov(obj.spk_nomean);
            covmat(eye(size(covmat))==1) = obj.varce;
            
            corrmat = corr(obj.spk_nomean,'rows','pairwise');
            
            % test if positive defined
            % [~,p] = chol(corrmat);
            
            CorCE = covmat./sqrt(obj.varce'*obj.varce);
            
            obj.corce = CorCE;
            obj.covce = covmat;
            if plot_flag
                figure
                subplot(1,3,1)
                imagesc(CorCE)
                colorbar
                title('CorCE')
                subplot(1,3,2)
                imagesc(covmat)
                title('CovCE')
                subplot(1,3,3)
                imagesc(corrmat)
                title('corr.')
                set(gcf,'Position',[87  322  955  277])
                
                colormap(jet)
            end
            
            
        end
        
        function auto_regression(obj)
            cn = obj.cn;
            % idx = obj.x(:,cn.R_TRG)==1 & obj.x(:,cn.R_COH)==0;
            idx = ones(obj.ntr,1)==1;
            %             idx = obj.x(:,cn.R_NTARGS)==2;
            
            y = obj.spk_nomean(idx,end);
            x = [obj.spk_nomean(idx,1:end-1) obj.x(idx,cn.R_COH) ...
                obj.x(idx,cn.R_TRG) obj.x(idx,cn.R_DIR) ones(sum(idx),1)];
            
            inds = not(isnan(y));
            [ball,dev,stats] = glmfit(x(inds,:),y(inds),'normal','link','identity','constant','off');
            
            % now the same thing, for each time step at a time
            clear b
            for i = 1:obj.ntimes-1
                x = [obj.spk_nomean(idx,i) obj.x(idx,cn.R_COH) ...
                    obj.x(idx,cn.R_TRG) obj.x(idx,cn.R_DIR) ones(sum(idx),1)];
                btemp = glmfit(x(inds,:),y(inds),'normal','link','identity','constant','off');
                b(i) = btemp(1);
            end
            
            
            p = publish_plot(1,2);
            subplot(1,2,1)
            inds = 1:obj.ntimes-1;
            plot(obj.times(inds),ball(inds),'.-')
            hold on
            plot(obj.times(inds),b(inds),'r.-')
            axis tight
            pxlim
            
            xlabel('time (ms)')
            ylabel('\beta')
            
            cn = obj.cn;
            % idx = (obj.x(:,cn.R_TRG)==1 & obj.x(:,cn.R_DIR)==1);
            % idx = ones(obj.ntr,1)==1;
            timezero = find(obj.timesrt == 0);
            
            y = obj.spk_rt_nomean(idx,timezero);
            x = [obj.spk_rt_nomean(idx,1:timezero-1) obj.x(idx,cn.R_COH) ...
                obj.x(idx,cn.R_TRG) obj.x(idx,cn.R_DIR) ones(sum(idx),1)];
            
            inds = not(isnan(y));
            [b,dev,stats] = glmfit(x(inds,:),y(inds),'normal','link','identity','constant','off');
            
            subplot(1,2,2)
            plot(obj.timesrt(1:timezero-1),b(1:timezero-1),'.-')
            axis tight
            pxlim
            xlabel('time to rt (ms)')
            ylabel('\beta')
            set(gca,'yaxislocation','right')
            same_ylim(p.h_ax)
            p.format('FontSize',15);
            
        end
        
    end
end