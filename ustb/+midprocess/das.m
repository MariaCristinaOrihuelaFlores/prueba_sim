classdef das < midprocess
    %DAS   Implementation of USTB DAS general beamformer
    %
    %   authors:    Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
    %               Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %               Stefano Fiorentini <stefano.fiorentini@ntnu.no>
    %
    %   $Last updated: 01/02/2023$
    
    %% Additional properties
    properties
        dimension = dimension.both;         % Dimension enumeration class that specifies whether the process will run only on transmit, receive, both, or none.
        code = code.mex;                    % Code enumeration class that specifies the code to be run (code.matlab, code.mex, ...)
        gpu_id = 0;                         % In case code.matlab_gpu or code.mex_gpu is used, this variable specfifies which CUDA-enabled gpu the code should run on 
        
        spherical_transmit_delay_model = spherical_transmit_delay_model.hybrid; % spherical transmit delay model enumeration for deciding model when the source is in front of the transducer
        pw_margin = 1e-3;                   % The margin of the area around focus in m for the spherical_transmit_delay_model.hybrid
        
        transmit_delay                      % Variable returning the calculated tx part of the receive delay so that it can be plotted
        receive_delay                       % Variable returning the calculated rx part of the receive delay so that it can be plotted

        elapsed_time                        % Variable to store the beamforming time. Used for benchmarking
    end
    
    %% constructor
    methods (Access = public)
        function h=das()
            h.name='USTB DAS General Beamformer';
            h.reference= 'www.ustb.no';
            h.implemented_by={'Stefano Fiorentini <stefano.fiorentini@ntnu.no>', 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>','Ole Marius Hoel Rindal <olemarius@olemarius.net>'};
            h.version='v1.1.0';
        end
    end
    
    %% go method
    methods
        function beamformed_data=go(h)
            
            % short names
            N_pixels = h.scan.N_pixels;
            N_channels = h.channel_data.N_channels;
            N_waves = h.channel_data.N_waves;
            N_frames = h.channel_data.N_frames;
            
            % constants
            sampling_frequency=single(h.channel_data.sampling_frequency);
            initial_time=single(h.channel_data.initial_time);
            modulation_frequency=single(h.channel_data.modulation_frequency);
            w0=2*pi*modulation_frequency;

            % calculate transmit apodization according to 10.1109/TUFFC.2015.007183
            h.transmit_apodization.sequence=h.channel_data.sequence;
            h.transmit_apodization.focus=h.scan;
            tx_apodization=single(h.transmit_apodization.data);
            
            % calculate receive apodization
            h.receive_apodization.probe=h.channel_data.probe;
            h.receive_apodization.focus=h.scan;
            rx_apodization=single(h.receive_apodization.data);
            
            % calculate receive delay
            xm=bsxfun(@minus,h.channel_data.probe.x.',h.scan.x);
            ym=bsxfun(@minus,h.channel_data.probe.y.',h.scan.y);
            zm=bsxfun(@minus,h.channel_data.probe.z.',h.scan.z);
            receive_delay=single(sqrt(xm.^2+ym.^2+zm.^2)/h.channel_data.sound_speed); %#ok<*PROP> 
            h.receive_delay = receive_delay;
            
            % calculate transmit delay
            transmit_delay=zeros(N_pixels,N_waves);
            for n_wave=1:numel(h.channel_data.sequence)
                switch(h.channel_data.sequence(n_wave).wavefront)
                    % point source
                    case uff.wavefront.spherical
                        % check if the point source is at infinity -> for backcompatibility of examples
                        if isinf(h.channel_data.sequence(n_wave).source.distance)
                            transmit_delay(:,n_wave)=h.scan.z*cos(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.x*sin(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.y*sin(h.channel_data.sequence(n_wave).source.elevation);
                        else
                            % distance between source and elements
                            transmit_delay(:,n_wave)=(-1).^(h.scan.z<h.channel_data.sequence(n_wave).source.z).*sqrt((h.channel_data.sequence(n_wave).source.x-h.scan.x).^2+(h.channel_data.sequence(n_wave).source.y-h.scan.y).^2+(h.channel_data.sequence(n_wave).source.z-h.scan.z).^2);
                            
                            % add distance from source to origin
                            if (h.channel_data.sequence(n_wave).source.z<0) %Diverging Wave (DW) transmit
                                transmit_delay(:,n_wave)=transmit_delay(:,n_wave)-abs(h.channel_data.sequence(n_wave).source.distance);
                            else % if virtual source in front of transducer
                                switch h.spherical_transmit_delay_model
                                    % Please see the reference below for a documentation and definition of these three models and the differences.
                                    % Rindal, O. M. H., Rodriguez-Molares, A., & Austeng, A. (2018). A simple , artifact-free , virtual source model. 
                                    % IEEE International Ultrasonics Symposium, IUS, 1–4.        
                                    case spherical_transmit_delay_model.spherical
                                        % Use conventional virtual source model
                                        transmit_delay(:,n_wave)=transmit_delay(:,n_wave)+h.channel_data.sequence(n_wave).source.distance;
                                    case spherical_transmit_delay_model.unified
                                        % Transmit delay model introduced in  Nguyen, N. Q., & Prager, R. W. (2016). High-Resolution Ultrasound 
                                        % Imaging With Unified Pixel-Based Beamforming. IEEE Trans. Med. Imaging, 35(1), 98-108.
                                        if n_wave == 1
                                            % Get an apodization mask to identify the four different regions of the transmit wave field
                                            mask_apod = uff.apodization();
                                            mask_apod.window = uff.window.boxcar;
                                            mask_apod.sequence = h.channel_data.sequence;
                                            mask_apod.minimum_aperture = [0,0];
                                            mask_apod.focus = h.scan;
                                            if isa(h.scan,'uff.sector_scan')
                                                mask_apod.f_number = h.transmit_apodization.f_number; %This should be set according to the actually transmitted f number
                                                mask_all_waves = reshape(mask_apod.data,h.scan.N_depth_axis,h.scan.N_azimuth_axis,numel(h.channel_data.sequence));
                                            elseif isa(h.scan,'uff.linear_scan')
                                                mask_apod.f_number = h.transmit_apodization.f_number; %This should be set according to the actually transmitted f number
                                                mask_all_waves = reshape(mask_apod.data,h.scan.N_z_axis,h.scan.N_x_axis,numel(h.channel_data.sequence));
                                            else
                                                error('Only linear scan and sector scan in 2D is supported for the unified spherical transmit delay model.');
                                            end
                                        end
                                        transmit_delay(:,n_wave) = tools.calculate_unified_delay_model(transmit_delay(:,n_wave),logical(mask_all_waves(:,:,n_wave)),...
                                                                            h.scan,h.channel_data.sequence(n_wave).source);
                                    case spherical_transmit_delay_model.hybrid
                                        % Our hybrid model assume that the transmit wave propagates as a plane-wave in a small region h.pw_margin around the transmit focus
                                        % Calculate the "plane wave" in the transmit direction
                                        if isa(h.scan,'uff.linear_scan')
                                            plane_delay = (-1).^(h.scan.z<h.channel_data.sequence(n_wave).source.z).*sqrt((h.channel_data.sequence(n_wave).source.z-h.scan.z).^2) + h.channel_data.sequence(n_wave).source.distance;
                                        elseif isa(h.scan,'uff.sector_scan')
                                            plane_delay = h.scan.z*cos(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.x*sin(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.y*sin(h.channel_data.sequence(n_wave).source.elevation);
                                        else
                                            error('Only linear scan and sector scan in 2D is supported for the hybrid spherical transmit delay model.');
                                        end
                                        
                                        % Find region in front of and after the focus with margins given by h.margin_in_m
                                        z_mask = logical(h.scan.z < (h.channel_data.sequence(n_wave).source.z + h.pw_margin)) & logical(h.scan.z > (h.channel_data.sequence(n_wave).source.z - h.pw_margin));
                                        
                                        % Replace the region in the virtual source delay with the plane delay in
                                        % the region indicated by the mask
                                        transmit_delay_temp = transmit_delay(:,n_wave)+h.channel_data.sequence(n_wave).source.distance;
                                        transmit_delay_temp(z_mask) = plane_delay(z_mask);
                                        transmit_delay(:,n_wave) = transmit_delay_temp;
                                end
                            end  
                        end
                    % plane wave
                    case uff.wavefront.plane
                        transmit_delay(:,n_wave)=h.scan.z*cos(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.x*sin(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.y*sin(h.channel_data.sequence(n_wave).source.elevation);
                        
                    % photoacoustic wave
                    case uff.wavefront.photoacoustic
                        transmit_delay(:,n_wave)=zeros(N_pixels,1);
                        
                    otherwise
                        error('Unknown wavefront. Check available options at uff.wavefront');
                end
                % convert to seconds and include wave delay
                transmit_delay(:,n_wave) = transmit_delay(:,n_wave)./h.channel_data.sequence(n_wave).sound_speed - h.channel_data.sequence(n_wave).delay;
            end
            % convert to single
            transmit_delay = single(transmit_delay);
            
            % Saving receive tx and rx delay to a public parameter to be able to plot it
            h.transmit_delay = transmit_delay;
            
            ch_data=single(h.channel_data.data);
            if (abs(w0)<eps)
                ch_data = hilbert(ch_data);
            end
            
            % create beamformed data class
            h.beamformed_data=uff.beamformed_data();
            h.beamformed_data.scan=h.scan;
            
            % Check that enough RAM can be allocated to beamformed data
            switch h.dimension
                case dimension.none
                    tools.check_memory(prod([N_pixels, N_channels, N_waves, N_frames, 8]));
                case dimension.receive
                    tools.check_memory(prod([N_pixels, N_waves, N_frames, 8]));
                case dimension.transmit
                    tools.check_memory(prod([N_pixels, N_channels, N_frames, 8]));
                case dimension.both
                    tools.check_memory(prod([N_pixels, N_frames, 8]));
            end
            
            % Run Delay & Sum
            if any(ch_data, 'all') % only process if any data > 0
                switch h.code
                    case code.mex
                        str = "MEX C";
                    case code.mex_gpu
                        str = "MEX GPU";
                    case code.matlab
                        str = "MATLAB";
                    case code.matlab_gpu
                        str = "MATLAB GPU";
                end
                                        
                fprintf(1, "USTB %s beamformer...", str)
                tic()

                switch h.code
                    %% MEX
                    case code.mex
                        bf_data=mex.das_c(ch_data,...
                                          sampling_frequency,...
                                          initial_time,...
                                          tx_apodization,...
                                          rx_apodization,...
                                          transmit_delay,...
                                          receive_delay,...
                                          modulation_frequency,...
                                          int32(h.dimension));
                    %% MEX CUDA                  
                    case code.mex_gpu
                        bf_data=mex.das_cuda(ch_data,...
                                             sampling_frequency,...
                                             initial_time,...
                                             tx_apodization,...
                                             rx_apodization,...
                                             transmit_delay,...
                                             receive_delay,...
                                             modulation_frequency,...
                                             int32(h.dimension), ...
                                             int32(h.gpu_id));
                    
                    %% MATLAB
                    case code.matlab

                        bf_data = tools.matlab_beamformer(ch_data,...
                                                          single(h.channel_data.time),...
                                                          tx_apodization,...
                                                          rx_apodization,...
                                                          transmit_delay,...
                                                          receive_delay,...
                                                          w0,...
                                                          h.dimension);

                    
                    %% MATLAB GPU 
                    case code.matlab_gpu
                        bf_data = tools.matlab_gpu_beamformer(ch_data,...
                                                             single(h.channel_data.time),...
                                                             tx_apodization,...
                                                             rx_apodization,...
                                                             transmit_delay,...
                                                             receive_delay,...
                                                             w0,...
                                                             h.dimension, ...
                                                             h.gpu_id+1);
                        
                    otherwise
                        error('Unknown code implementation requested');
                end % end switch          

                h.elapsed_time = toc();
                fprintf(1, "Completed in %.2f seconds.\n", h.elapsed_time);
            end % end if

            % assign phase according to 2 times the receive propagation distance
            if (abs(w0) > eps)
                bf_data = bf_data .* exp(-1j*2*w0*h.scan.reference_distance/h.channel_data.sound_speed);
            end

            % copy data to object
            h.beamformed_data.data = bf_data;
            
            % pass a reference
            beamformed_data = h.beamformed_data;
            
        end % end go()
    end
end
