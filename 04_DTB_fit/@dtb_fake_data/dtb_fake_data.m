classdef dtb_fake_data < handle
    properties
        
        kappa = 15
        B0 = 1
        a = 0
        d = 0
        srho = -1
        ndt_mu = 0.2
        ndt_sigma = 0.01
        ndt_method_flag = 2
        coh0 = 0
        y0 = 0
        y0std = 0
        y0noise
        Brectif = -200 % a very negative val
        USfunc = 'None'
        boost = 1
        urgency_flag = 0


        tau_kernel_noise = 0.04 % seconds- kernel to smooth the noise
        tdelay_acum = 0.15 % seconds - delay between the ev stream and the acum
        tdelay_acum_std = 0.04 % std of tdelay acum

        p_start_with_right = 0.5
        
        std_noise = 1
        seed = nan
        rstream
        
        ntr_per_coh
        ntrials
        
        
        attention_shift_rate = 0; % shift attention exponentially with mean 1/attention_shift_rate
        delta_attention_shift = 0; %relative influence of the attention shift; 1:binary effect ;0: no effect
        attention_shift
        
        ucoh
        t = [0:0.001:6]
        
        noise
        momentary_evidence
        decision_variable
        ext_noise
        convolve_with_IR = true
        
        coh
        
        winner
        isWinner
        decision_time
        be
        dv1_atdt
        dv2_atdt
        dv1_atbe
        dv2_atbe
        ndt
        RT
        correct
        
        conf_extra_m
        conf_extra_s
        
        confidence
        
        Bup
    end
    
    methods
        function obj = dtb_fake_data(varargin)
            for i=1:2:length(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
            obj.ntrials = size(obj.coh,1);
            obj.ucoh = unique(obj.coh);
            
            obj.set_seed();
        end
        
        function set_seed(obj)
            if isnan(obj.seed)
                %aux = rng('shuffle');
                %aux = rng('shuffle', 'twister')
                obj.seed = sum(100*clock);
            end
            % obj.rstream = RandStream('twister','Seed',obj.seed);
            rng(obj.seed,'twister');
            obj.seed
        end
        
        function make_fakedata(obj)
            
            ntimes = length(obj.t);
            
            use_chole = false;
            if (use_chole) %this has a problem with the variance of the resulting noise ??
                ev_noise_cho = obj.std_noise*randn(obj.ntrials*ntimes,2);
                U = [1 obj.srho; obj.srho 1];
                noise = ev_noise_cho*U;
                obj.noise{1} = reshape(noise(:,1),[obj.ntrials,ntimes]);
                obj.noise{2} = reshape(noise(:,2),[obj.ntrials,ntimes]);
            else
                
                nt = ceil(ntimes * 1.2); % longer because of the kernel

                if 1
                    ev_noise = obj.std_noise * randn(obj.ntrials,nt,3); % first one is shared noise
                    
    
                    rho = [sqrt(abs(obj.srho)) sqrt(1-abs(obj.srho)) 0; ...
                        sign(obj.srho)*sqrt(abs(obj.srho)) 0 sqrt(1-abs(obj.srho))];
    
    
                    obj.noise = cell(2,1);
                    
                    obj.noise{1} = zeros(obj.ntrials,nt);
                    obj.noise{2} = zeros(obj.ntrials,nt);
                    
                    
                    for i=1:3
                        for j=1:2
                            obj.noise{j} = obj.noise{j} + squeeze(ev_noise(:,:,i) * rho(j,i));
                        end
                    end

                else % testing
                    mu = [0 0]; % Mean (of the noise)
                    sigma = 1; % Standard deviation for both variables
                    rho = obj.srho; % Correlation coefficient
                    
                    nSamples = obj.ntrials*nt;
                    % Covariance matrix calculation
                    % Since Cov(X,Y) = rho * sigma_X * sigma_Y, and sigma_X = sigma_Y = sigma,
                    % the off-diagonal elements (covariance) are rho * sigma^2
                    
                    covarianceMatrix = [sigma^2, rho*sigma^2; rho*sigma^2, sigma^2];
                    
                    % Generate samples
                    ev_noise = mvnrnd(mu, covarianceMatrix, nSamples);
                    obj.noise = cell(2,1);
                    for i=1:2
                        obj.noise{i} = reshape(ev_noise(:,i), obj.ntrials, nt);
                    end

                end            

                % make the noise auto-correlated by convolving with
                % exponential kernel
                if obj.tau_kernel_noise~=0 
%                     stdev_pre = nanstd(obj.noise{1}(:));
                    kernel = exp(-obj.t/obj.tau_kernel_noise);
                    ind = find(kernel<eps,1);
                    kernel = kernel(1:ind);
                    kernel = kernel/sum(kernel);
                    for i=1:size(obj.noise)
                        aux = conv2(1,kernel,obj.noise{i},'valid');
                        % keep ntimes size:
                        obj.noise{i} = aux(:,1:ntimes,:);
                    end
                    % preserve the variance of the noise
%                     stdev_post = nanstd(obj.noise{1}(:));
%                     for i=1:size(obj.noise)
%                         obj.noise{i} = obj.noise{i} * stdev_pre/stdev_post;
%                     end
                else
                    for i=1:size(obj.noise)
                        % keep ntimes size:
                        obj.noise{i} = obj.noise{i}(:,1:ntimes,:);
                    end
                    
                end
                
                dt = obj.t(2) - obj.t(1);
                

            
            end



            % ev
            %obj.coh = repmat(obj.ucoh(:),obj.ntr_per_coh,1);
            %obj.coh = shuffle(obj.coh);
            evup = obj.kappa * dt * (repmat(obj.coh,1,ntimes));
            evlo = obj.kappa * dt * (repmat(-obj.coh,1,ntimes));
            
            % attention shifts
            %attention_shift_rate = 1; % shift exponentially with mean 1/attention_shift_rate
            
            if obj.attention_shift_rate==0
                
                s = {evlo - obj.kappa * dt * obj.coh0, evup + obj.kappa * dt * obj.coh0};
                s{1}(:,end) = 1000;%to always get a response
                s{2}(:,end) = 1000;%to always get a response
                
                
            elseif obj.boost==1
                
                attention_shift = cumsum(poissrnd(obj.attention_shift_rate * dt,size(evup)),2);
                %so that the initial focus is random:
                initial_focus = double(rand(obj.ntrials,1)>obj.p_start_with_right);
                attention_shift = bsxfun(@plus,attention_shift,initial_focus);
                attention_shift = mod(attention_shift,2);%which one the model is "looking"
                obj.attention_shift = attention_shift;
                
                
                
                attention_shift_up = ones(size(attention_shift)) + ...
                    obj.delta_attention_shift.*attention_shift - obj.delta_attention_shift.*(1-attention_shift);
                
                attention_shift_lo = ones(size(attention_shift)) + ...
                    obj.delta_attention_shift.*(1-attention_shift)-obj.delta_attention_shift.*attention_shift;
                
                
                s = {evup.*attention_shift_up - evlo.*attention_shift_lo + obj.coh0, ...
                     evlo.*attention_shift_lo - evup.*attention_shift_up - obj.coh0};
                s{1}(:,end) = 1000;%to always get a response
                s{2}(:,end) = 1000;%to always get a response
                
            else
                
                attention_shift = cumsum(poissrnd(obj.attention_shift_rate * dt,size(evup)),2);
                %so that the initial focus is random:
                initial_focus = double(rand(obj.ntrials,1)>obj.p_start_with_right);
                attention_shift = bsxfun(@plus,attention_shift,initial_focus);
                attention_shift = mod(attention_shift,2);%which one the model is "looking"
                obj.attention_shift = attention_shift;
                
                
                
                attention_shift_up = attention_shift;
                
                attention_shift_lo = 1-attention_shift;
                
                evup_at_up = evup-evlo*obj.boost + obj.coh0;
                evup_at_lo = evup*obj.boost-evlo + obj.coh0;
                
                evlo_at_up = evlo*obj.boost-evup - obj.coh0;
                evlo_at_lo = evlo-evup*obj.boost - obj.coh0;
                
                s = {evup_at_up.*attention_shift_up + evup_at_lo.*attention_shift_lo, ...
                    evlo_at_up.*attention_shift_up + evlo_at_lo.*attention_shift_lo};
                
%                 s = {evup.*attention_shift_up - evlo.*attention_shift_lo*obj.boost + obj.coh0, ...
%                      evlo.*attention_shift_lo - evup.*attention_shift_up*obj.boost - obj.coh0};
                s{1}(:,end) = 1000;%to always get a response
                s{2}(:,end) = 1000;%to always get a response
                
            end
            
            % one bias
            obj.y0noise = randn(obj.ntrials,2)*obj.y0std;
            s{1}(:,1) = s{1}(:,1) + obj.y0*obj.B0 + obj.y0noise(:,1);
            % s{2}(:,1) = s{2}(:,1) - obj.y0 + obj.y0noise(:,2);
            s{2}(:,1) = s{2}(:,1) + obj.y0*obj.B0 + obj.y0noise(:,2);
            for i=1:2
                s{i} = s{i} + obj.noise{i}*sqrt(dt);
            end
            %
            obj.momentary_evidence = s;
            
        end
        
        function gaze_post_decision(obj, p_switch_to_chosen_one)
            
            if nargin==0
                p_switch_to_chosen_one = 0;
            end
            dt = obj.t(2)-obj.t(1);
            dt_samp = ceil(obj.decision_time/dt);
%             rt_samp = ceil(obj.RT/dt);
            p = rand(obj.ntrials,1)<p_switch_to_chosen_one;
            nt = length(obj.t);
            for i=nt:-1:1
                I = p & i>dt_samp;
                obj.attention_shift(I,i) = obj.winner(I) - 1;
%                 I = p & i>rt_samp;
%                 obj.attention_shift(I,i) = nan;
            end
            
            
        end
            
            
        
        function add_external_noise(obj,sigma_ext,repeat_each_column)
            % creates and adds external noise
            
            dt = obj.t(2)-obj.t(1);
            
            s = obj.momentary_evidence;
            [ntr,nt] = size(s{1});
            for i=1:3
                ext_noise{i} = randn(ntr,ceil(nt/repeat_each_column))*sqrt(dt)*sigma_ext;
            end
            
            for i=1:3 %first is shared noise
                ext_noise{i} = repeat_columns(ext_noise{i},repeat_each_column);
                ext_noise{i} = ext_noise{i}(:,1:nt);
                
%                 s{i} = s{i} + obj.ext_noise{i};
            end
            
%             obj.momentary_evidence = s;
            
            % use the same anti-correlation as for the races
            rho = [sqrt(abs(obj.srho)) sqrt(1-abs(obj.srho)) 0; ...
                sign(obj.srho)*sqrt(abs(obj.srho)) 0 sqrt(1-abs(obj.srho))];
            obj.ext_noise{1} = zeros(obj.ntrials,nt);
            obj.ext_noise{2} = zeros(obj.ntrials,nt);
            for i=1:3
                for j=1:2
                    obj.ext_noise{j} = obj.ext_noise{j} + ...
                        squeeze(ext_noise{i} * rho(j,i));
                end
            end
            
            
            if obj.convolve_with_IR
                IR = load('impulse_resp');
                % interpolate time
                IRr = interp1(IR.t,IR.me,obj.t * 1000);
                IRr = IRr(1:find(IRr>0,1,'last'));
                IRr = IRr(:);
                for i=1:2
                    nn = conv2(obj.ext_noise{i},IRr')/sum(IRr);
                    nn = nn(:,1:length(obj.t));
                    obj.ext_noise{i} = nn;
                end
            end
            
            for i=1:2
                % obj.ext_noise{i} = obj.ext_noise{i};
                s{i} = s{i} + obj.ext_noise{i};
            end
            
            obj.momentary_evidence = s;
        end
        
        function remove_weak_evidence(obj,thres)
            % for robust integration
            for i=1:length(obj.momentary_evidence)
                ev = obj.momentary_evidence{i};
                ev(abs(ev)<thres) = 0;
                obj.momentary_evidence{i} = ev;
            end
            
        end
        
        function diffuse_to_bound(obj)
            
            ev = obj.momentary_evidence;
            % make the temporal delay
            if obj.tdelay_acum~=0
                dt = obj.t(2)-obj.t(1);
                ind = ceil(obj.tdelay_acum/dt + obj.tdelay_acum_std/dt * randn(obj.ntrials,1));
                ind(ind<0) = 0;
                for i=1:length(ev)
                    aux = zeros(size(ev{i}));
                    for j=1:obj.ntrials
                        aux(j, ind+1:end) = ev{i}(j, 1:end-ind);
                    end
                    aux(:,end) = ev{i}(:,end); % keep the high value at the end so that it terminates

                    ev{i} = aux;
                end
            end

            if obj.urgency_flag==1
                [obj.winner,obj.isWinner,obj.decision_time,obj.Bup,cev] = ...
                    dtb_multi_urgency(ev,obj.t,obj.a,obj.d,obj.B0,obj.Brectif,obj.USfunc);
            else
                [obj.winner,obj.isWinner,obj.decision_time,obj.Bup,cev] = ...
                    dtb_multi(ev,obj.t,obj.a,obj.d,obj.B0,obj.Brectif,obj.USfunc);
            end
            obj.decision_variable = cev;
            
            obj.ndt = calc_ndt(obj.ndt_mu,obj.ndt_sigma,obj.ntrials,obj.ndt_method_flag);
            
            obj.RT = obj.decision_time + obj.ndt;
            
            %accuracy
            model_correct = nan(obj.ntrials,1);
            model_correct(obj.coh>0 & obj.winner==2 | obj.coh<0 & obj.winner==1) = 1;% correct
            model_correct(obj.coh>0 & obj.winner==1 | obj.coh<0 & obj.winner==2) = 0;% incorrect
            model_correct(obj.coh==0) = rand(sum(obj.coh==0),1)>0.5;% 
            obj.correct = model_correct;
            
            obj.confidence = [];
            
        end
        
        function calc_balance_of_evidence(obj,conf_extra_m,conf_extra_s)
            
            if nargin<3 || isempty(conf_extra_m)
                ndt_extra = zeros(obj.ntrials,1);
            else
                obj.conf_extra_m = conf_extra_m;
                obj.conf_extra_s = conf_extra_s;
                ndt_extra = calc_ndt(conf_extra_m,conf_extra_s,obj.ntrials,obj.ndt_method_flag);
            end
            
            %balance of evidence
            [obj.be,obj.dv1_atbe,obj.dv2_atbe] = calc_balance_evidence(obj.decision_variable,...
                obj.decision_time + ndt_extra,obj.t,obj.winner,obj.y0noise);
            
        end
        
        function yhat = calc_confidence_logistic_dt(obj)
            % confidence
            % as a function of decision time
            X = [obj.decision_time,ones(size(obj.decision_time))];
            b = glmfit(X,obj.correct,'binomial','link','logit','constant','off');
            yhat = glmval(b,X,'logit','constant','off');
            obj.confidence = yhat;
        end
        
        function yhat = calc_confidence_logistic_dt_be(obj)
            % confidence
            % as a function of decision time
            X = [obj.decision_time,obj.be,ones(size(obj.decision_time))];
            b = glmfit(X,obj.correct,'binomial','link','logit','constant','off');
            yhat = glmval(b,X,'logit','constant','off');
            obj.confidence = yhat;
        end
        
        function yhat = calc_confidence_logistic_dt_be_attention(obj)
            % confidence
            % as a function of decision time
            
            attention = motionenergy.remove_post_decision_samples(obj.attention_shift,obj.t,obj.decision_time);
            attention(obj.winner==2,:) = 1 - attention(obj.winner==2,:);
            prop_att = nanmean(attention,2);
            time_att = prop_att.*obj.decision_time;
            
            X = [obj.decision_time,obj.be,ones(size(obj.decision_time)),time_att,prop_att];
            b = glmfit(X,obj.correct,'binomial','link','logit','constant','off');
            yhat = glmval(b,X,'logit','constant','off');
            obj.confidence = yhat;
        end
        
        function yhat = calc_confidence_logistic_be(obj)
            % confidence
            % as a function of decision time
            X = [obj.be,ones(size(obj.decision_time))];
            b = glmfit(X,obj.correct,'binomial','link','logit','constant','off');
            yhat = glmval(b,X,'logit','constant','off');
            obj.confidence = yhat;
        end
        
        function yhat = calc_confidence_logistic_be_attention(obj)
            % confidence
            % as a function of decision time
            
            attention = motionenergy.remove_post_decision_samples(obj.attention_shift,obj.t,obj.decision_time);
            attention(obj.winner==2,:) = 1 - attention(obj.winner==2,:);
            prop_att = nanmean(attention,2);
            time_att = prop_att.*obj.decision_time;
            
            X = [obj.be,ones(size(obj.decision_time)),time_att,prop_att];
            b = glmfit(X,obj.correct,'binomial','link','logit','constant','off');
            yhat = glmval(b,X,'logit','constant','off');
            obj.confidence = yhat;
        end
        
        function do_basic_plots(obj)
            p = publish_plot(1,3);
            
            %choice
            p.next()
            curva_media(obj.winner==1,obj.coh,[],1);
            xlabel('coh');ylabel('choice');
            
            %RT
            p.next()
            curva_media(obj.RT,obj.coh,obj.correct==1,1);
            hold all
            curva_media(obj.RT,obj.coh,obj.correct==0,1);
            xlabel('coh');ylabel('RT');
            
            %confidence
            if ~isempty(obj.confidence)
                p.next()
                curva_media(obj.confidence,obj.coh,obj.correct==1,1);
                hold all
                curva_media(obj.confidence,obj.coh,obj.correct==0,1);
                xlabel('coh');ylabel('confidence');
            end
            p.format('FontSize',16);
            set(p.h_ax,'xlim',[-1*max(obj.ucoh),max(obj.ucoh)]);
            set(gcf,'Position',[253  925  950  228])
            
            figure()
            hist(obj.RT,200)
            
            %kernels
%             if ~isempty(obj.ext_noise)
            if 0
                p = publish_plot(1,2);
%                 dt = obj.t(2)-obj.t(1);
                inds = obj.winner==1;
                sel(inds,:)  = obj.ext_noise{1}(inds,:);
                nsel(inds,:) = obj.ext_noise{2}(inds,:);
                inds = obj.winner==2;
                sel(inds,:)  = obj.ext_noise{2}(inds,:);
                nsel(inds,:) = obj.ext_noise{1}(inds,:);
                
                rt_sample   = round(obj.RT/((obj.t(2)-obj.t(1))));
                sel         = motionenergy.remove_post_decision_samples(sel,rt_sample);
                nsel        = motionenergy.remove_post_decision_samples(nsel,rt_sample);
                [rt_sel,rt_t_sel]   = eventlockedmatc(sel,obj.t,obj.RT,[1,1]);
                [rt_nsel,rt_t_nsel] = eventlockedmatc(nsel,obj.t,obj.RT,[1,1]);
                
                %noise free
                %             high = obj.confidence>nanmedian(obj.confidence);
                
                %with noise
                if ~isempty(obj.confidence)
                    addconfnoise_flag = false;
                    if addconfnoise_flag
                        a = 50;
                        y = 1./(1+exp(-a*(obj.confidence-nanmedian(obj.confidence))));
                        high = rand(size(y))<y;
                    else
                        high = obj.confidence>nanmedian(obj.confidence);
                    end
                end
                p.next();
                k1 = nanmean(sel);
                k2 = nanmean(nsel);
                %             k1 = conv(k1,ones(round(0.05/dt),1),'same');
                %             k2 = conv(k2,ones(round(0.05/dt),1),'same');
                plot(obj.t,k1);
                hold all
                plot(obj.t,k2);
                
                if ~isempty(obj.confidence)
                    k1 = nanmean(sel(high,:)) - nanmean(sel(~high,:));
                    k2 = nanmean(nsel(high,:)) - nanmean(nsel(~high,:));
                    %             k1 = conv(k1,ones(round(0.05/dt),1),'same');
                    %             k2 = conv(k2,ones(round(0.05/dt),1),'same');
                    plot(obj.t,k1);
                    hold all
                    plot(obj.t,k2);
                end
                xlim([0,prctile(obj.RT,70)]);
                
                p.next();
                plot(rt_t_sel,nanmean(rt_sel,2))
                hold all
                plot(rt_t_nsel,nanmean(rt_nsel,2))
                if ~isempty(obj.confidence)
                    plot(rt_t_sel,nanmean(rt_sel(:,high),2) - nanmean(rt_sel(:,~high),2))
                    plot(rt_t_sel,nanmean(rt_nsel(:,high),2) - nanmean(rt_nsel(:,~high),2))
                end
                xlim([-0.5,0])
                
                same_ylim(p.h_ax)
                
                p.current_ax(1);
                legend('sel','nsel','sel,H-L','nsel,H-L');
                
                set(p.h_ax(2),'yaxislocation','right')
                p.format('FontSize',15);
                
                %% just rightward choice, split by coh 
                p = publish_plot(1,2);
                w = 1;
                ext_noise_diff = obj.ext_noise{1} - obj.ext_noise{2};
                I = obj.winner==w;
                p.next();
                plot(obj.t,nanmean(ext_noise_diff(I,:)))
                hold all
                plot(obj.t,nanmean(ext_noise_diff(I & high==1,:)))
                plot(obj.t,nanmean(ext_noise_diff(I & high==0,:)))
                xlim([0,1])
                legend('all','high','low');
                
                p.next();
                [datamp,timeslocked] = eventlockedmatc(ext_noise_diff,obj.t,obj.RT,[1,1]);
                plot(timeslocked,nanmean(datamp(:,I),2));
                hold all
                plot(timeslocked,nanmean(datamp(:,I & high==1),2));
                plot(timeslocked,nanmean(datamp(:,I & high==0),2));
                xlim([-0.5,0])
                
                set(p.h_ax(2),'yaxislocation','right')
                p.format('FontSize',15);
                set(gcf,'Position',[487  332  938  341])
                same_ylim(p.h_ax)
                
            end
            
            %             set(gcf,'Position',[316  796  859  235])
            
            p = publish_plot(1,2);
            p.next()
            dotsanalysis.choice_against_duration_var(obj.correct,obj.coh,obj.RT,10,1,[])
            p.next()
            if ~isempty(obj.confidence)
                dotsanalysis.choice_against_duration_var(obj.confidence,obj.coh,obj.RT,10,1,obj.correct==1)
            end
            set(gcf,'Position',[248  942  772  269])
            p.format('FontSize',16);
            
        end
        
        function plot_kernels(obj)
            motion_energy_residuals = obj.ext_noise{1}-obj.ext_noise{2};
            
            filt = abs(obj.coh)<=0.064;
            t = obj.t;
            rt = obj.RT(filt);
            prc = prctile(obj.confidence,25);
            conf = obj.confidence(filt)>prc;
            choice = obj.winner(filt)==1;
            plot_kernels(motion_energy_residuals(filt,:),t,rt,choice,conf);

        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% helper functions %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndt = calc_ndt(ndt_mu,ndt_sigma,ntrials,method_flag)

% method_flag = 1;
% method_flag = 1;
switch method_flag
    case 1
        m = ndt_mu;
        s = ndt_sigma;
        v = s.^2;
        
        mm  = log((m^2)/sqrt(v+m^2));
        ss  = sqrt(log(v/(m^2)+1));
        
        ndt = lognrnd(mm,ss,[ntrials,1]);
        
    case 2
        
        m = ndt_mu;
        s = ndt_sigma;
        
        ndt = randn([ntrials,1])*s + m;
        while any(ndt<0)
            ndt(ndt<0) = randn(sum(ndt<0),1)*s + m;
        end
end


end


