%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Calc_chi_eps_chameleon_AP
%
% Modified from average_data_gen1.m . Trying to make a function to compute
% chameleon chi,eps etc. for a specified bin instead of the standard 1m
% bins. This is so I can compute these quantites patch-wise.
%
% 10/31/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

inds=1:100;



for n=1:nmax %step through all depth bins
                    nu=sw_visc(sal(n),temp(n),pres(n));  % calculate viscosity
                    % 	  eval(['avg.EPSILON' prb '(n)=calc_epsilon(cal.S' prb '(1+(min_ind(n)-1)*head.irep.S' prb ...
                    % 										  ':max_ind(n)*head.irep.S' prb '),avg.FALLSPD(n),nfft,nu,' ...
                    % 								   'head.sensor_index.S' prb ');'])
                    %
                    Sdata = getfield(cal,['S' prb]); %find shear probe 1 or 2
                    sensorindex = getfield(head.sensor_index,['S' prb]);
                    irep = getfield(head.irep,['S' prb]);
                    samplerate = head.slow_samp_rate.*irep;
                    % cut-off and filter order for the anti-aliasing Butterworth
                    % filter.  Needed to correct shear spectrum.
                    filtord = 4;
                    fcutoff = head.filter_freq(sensorindex);
                    fallspd = avg.FALLSPD(n);
                    sind = [1+(min_ind(n)-1)*irep:max_ind(n)*irep];
                    [f,ss] = calc_shearspec(Sdata(sind),samplerate,nfft,0.01*fallspd,filtord,fcutoff);
                    eval(['set_filters_' q.script.prefix]);
                    %variables FILTERS K_START,K_STOP,KSTOP_COEF that used in the next statement are calculated
                    %in set_filters script
                    if avg.FALLSPD(n)<200
                        eval(['avg.EPSILON' prb '(n)=calc_epsilon_filt1(ss,f,avg.FALLSPD(n),nfft,nu,' 'head.sensor_index.S' prb...
                            ',filters,k_start,k_stop,kstop_coef);'])
                    else
                        eval(['avg.EPSILON' prb '(n)=NaN;'])
                    end
                    %       eval(['avg.EPSILON' prb '(n)=calc_epsilon_filt_gen(cal.S' prb '(1+(min_ind(n)-1)*head.irep.S' prb ...
                    %                ':max_ind(n)*head.irep.S' prb '),avg.FALLSPD(n),nfft,nu,' 'head.sensor_index.S' prb...
                    %                ',filters,k_start,k_stop,kstop_coef);'])
                end % nmax (looping =over depth bins)
                
            elseif strncmp(active,'CHI',3)
                from=active(4:length(active));
                if isempty(from), from='1';,from1='';
                else from1=from;
                end
                % use the entire suffix if it is not a single number, otherwise
                if length(from)==1
                    if any(strcmp((avail_series),['T' from 'P']))
                        % use TXP if it is available
                        from=['T' from 'P'];
                    elseif any(strcmp((avail_series),['TP' from ]))
                        % use TPX if it is available
                        from=['TP' from ];
                    else
                        % use TP as a last resort
                        from='TP';
                    end
                end
                avg=select_epsilon(avg,epsilon_glitch_factor);
                for n=1:nmax
                    eval(['	avg.CHI' from1 '(n)=calc_chi_AP(cal.' from '(1+(min_ind(n)-1)*head.irep.' from  ...
                        ':max_ind(n)*head.irep.' from '),avg.FALLSPD(n),avg.EPSILON(n),nfft,sw_visc(sal(n),temp(n),pres(n)),' ...
                        'sw_tdif(sal(n),temp(n),pres(n)),' ...
                        ' head.sensor_index.' from ');'])
                end
            end
        end
        
    else % average data that is NOT epsilon or chi
        for n=1:nmax
            eval(['avg.' active '(n)=mean(cal.' active ...
                '(1+(min_ind(n)-1)*head.irep.' active ...
                ':max_ind(n)*head.irep.' active '));']);
        end
    end
end

%%