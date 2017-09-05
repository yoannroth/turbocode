classdef TurboCode
    %TurboCode Class with main methods or properties of a turbo coded
    %system
    
    %% Parameters which can be controlled by user
    properties (SetAccess = private)     
        backup          = true      % Save data and parameters to a mat file (default true)     
        blkLength       = 1024      % Frame size (default 1024) 
        k               = 4         % Constraint length (default 4) 
        poly            = [13 15]   % Polynomials [Feedback Generator] (default [13 15])    
        puncturing      = false     % Activate puncturing (default false)
        algorithm       = 'MAP'     % Decoding algorithm (MAP or maxlog) (default MAP)
        numIter         = 10        % Number of decoding iterations (default 10)  
    end
    
    %% Internal parameters of the system
    properties (SetAccess = private)    % Properties set by the constructor
        rate						% Rate of the code
        idStr						% Simulation ID
        saveStr						% String for saves
    end
    properties (GetAccess = private, Hidden = true)    % Internal properties
        interleaverIdx				% Indexes for the internal interleaver
        nextStates					% Next states indices for the trellis
        prevStates					% Previous states indices for the trellis
        codewords					% Possible codewords
        nOutputBits					% Number of output bits
    end
    
    %% Methods
    methods (Access = public)
		function obj = TurboCode(varargin)
		  %TURBOCODE Build object TurboCode using user parameters
		  % 
		  %		obj = TurboCode('paramName', paramValue, ...)	  where 
          % paramNam is the parameter to be selected and paramValue the 
          % value to be assigned. Default values can be edited in the 
          % class.

            %% Retrieve User values
            if ~isempty(varargin)
                for idxArg = 1:2:length(varargin)
                    if isprop(obj, varargin{idxArg})
                        % Check if parameter has correct type
                        eval(['classCurrent = class(obj.' varargin{idxArg} ');']);
                        if strcmp(classCurrent, class(varargin{idxArg+1}))
                            % Evaluate depending of parameter form
                            switch class(varargin{idxArg+1})
                                case 'double'
                                    eval(['obj.' varargin{idxArg} ' = [' num2str(varargin{idxArg+1}) '];'])
                                case 'char'
                                    eval(['obj.' varargin{idxArg} ' = ''' varargin{idxArg+1} ''';'])
                                case 'logical'
                                    eval(['obj.' varargin{idxArg} ' = ' num2str(varargin{idxArg+1}) ';'])
                                otherwise
                                    error(['Parameter modification not defined for this type of argument (' ...
                                        varargin{idxArg} ')'])
                            end
                        else  
                            error(['Wrong type for property "' varargin{idxArg} '", ' classCurrent ' expected']);
                        end
                    else
                        warning('CommSystem:badArgument', [' "' varargin{idxArg} '" was specified but is not a valid property'])
                    end
                end
            end
            
            %% Check property values
            if ~(obj.blkLength > 0)
                error('Property blkLength: Select a strictly positive value for blkLength');
            end
            if ~strcmp(obj.algorithm, 'MAP') && ~strcmp(obj.algorithm, 'maxlog')
                error('Property algorithm: Selected algorithm not defined (existing: MAP, maxlog)');
            end
            if length(obj.poly) ~= 2
                error('Property poly: Only 1/2 RSC codes are supported');
            end
            if obj.numIter < 1
                error('Property numIter: There must be at least 1 iteration');
            end
                
            %% Compute internal properties
            % Interleaver indices
            if obj.blkLength>5114 || obj.blkLength<40
                % If longer block length, random interleaving is used
                obj.interleaverIdx = randperm(obj.blkLength);
            else
                obj.interleaverIdx = get_UMTS_interleaver(obj.blkLength)';
            end

            % Get Trellis & codewords
            obj = obj.getRSCTrellis();

            % Get the previous states
            [~, obj.prevStates] = sort(obj.nextStates, 2);

            % Get code rate
            obj.nOutputBits     = (3*obj.blkLength + 4*(obj.k-1));
            obj.rate            = obj.blkLength/obj.nOutputBits;  

            % Puncturing : changes the rate
            if obj.puncturing
                obj.rate    = obj.rate*3/2;
                puncStr     = 'punc';
            else
                puncStr     = [];
            end
            
            % String name
            obj.saveStr = ['tc_' num2str(obj.poly(1)) '_' num2str(obj.poly(2)) ...
                puncStr '_N' num2str(obj.blkLength) '_' obj.algorithm];
            obj.idStr   = ['Turbo Code [' num2str(obj.poly(1)) ' ' num2str(obj.poly(2)) ...
                '], rate ' num2str(obj.rate) ' with ' obj.algorithm ...
                ' and ' num2str(obj.blkLength) ' bits'];
            
            
%             % Compile mex file
%             if ~exist('ForwardBackward.mexa64', 'file')
%                 mex ForwardBackward.c
%             end
            
        end
        
        function [ber, bler] = computeBer(obj, ebno, varargin)
        %computeBer Compute BER and BLER for the given ebno values 
        % and using the parameters of the object.
        %
        % [ber, bler] = computeBer(obj, ebno)  where ebno are the values of
        % ebno to use. For each value of ebno, simulations run until 1000
        % errors (default value) are found. Simulations stop if the ber is
        % under 10^(-6) (default value).
        %
        % [ber, bler] = computeBer(obj, ebno, maxNumErrs)  where ebno are 
        % the values of ebno to use. For each value of ebno, simulations 
        % run until maxNumErrs errors are found. Simulations stop if the 
        % ber is under 10^(-6) (default value).
        %
        % [ber, bler] = computeBer(obj, ebno, maxNumErrs, minBer)  where 
        % ebno are  the values of ebno to use. For each value of ebno, 
        % simulations run until maxNumErrs errors are found. Simulations 
        % stop if the  ber is under minBer.
        
            % Arg check
            switch length(varargin)
                case 0
                    maxNumErrs  = 1000;
                    minBer      = 1e-6;
                case 1
                    maxNumErrs  = varargin{1};
                    minBer      = 1e-6;
                case 2
                    maxNumErrs  = varargin{1};
                    minBer      = varargin{2};
                otherwise
                    error('Too many input argument')
            end
            % Sort ebno values
            ebno = sort(ebno);
            
            % Initialize
            snr         = obj.rate*10.^(ebno/10);   % SNR
            noiseVar    = 1./(snr*2);               % Noise variance
            maxNumBlks  = 1e7;                      % max tries
            
            % Outputs
            ber         = zeros(length(ebno), obj.numIter);
            bler        = zeros(length(ebno), obj.numIter);
            
            % Computations
            for varIdx = 1:length(noiseVar)
                numBlks = 0;
                while (ber(varIdx, end) < maxNumErrs && numBlks < maxNumBlks)
                    % Compute chain
                    errors = obj.computeChain(noiseVar(varIdx), false);
                    
                    % Compute BER/BLER
                    ber(varIdx, :)  = ber(varIdx, :) + errors;
                    bler(varIdx, :) = bler(varIdx, :) + (errors>0)*1.0;
                    numBlks = numBlks + 1;
                end
                
                ber(varIdx, :)  = ber(varIdx, :)/(obj.blkLength*numBlks);
                bler(varIdx, :) = bler(varIdx, :)/(numBlks);

                clc;
                %disp(strDisp);
                disp(['EbN0 = ' num2str(ebno(varIdx))]);
                disp(['BER  = ' num2str(ber(varIdx, end), '%1.3e')]);
                
                if ber(varIdx, end) < minBer
                    disp('Minimum BER reached: stopping simulation')
                end
                if obj.backup
                    % Save data script
                    obj.saveData('ber', ebno, ber, bler)
                end
            end
        end
        
        function [ia, ie] = computeExit(obj, ebno, varargin)
        %computeExit Compute EXIT for the given ebno value 
        % and using the parameters of the object.
        %
        % [ia, ie] = computeExit(obj, ebno)  where ebno is the value of
        % ebno to use. 10 (default value) simulations are run for values of
        % ia equal to 0:0.05:1 (default value)
        %
        % [ia, ie] = computeExit(obj, ebno, nbSim)  where ebno is the value
        % of ebno to use. nbSim simulations are run for values of ia equal 
        % to 0:0.05:1 (default value)
        %
        % [ia, ie] = computeExit(obj, ebno, nbSim, iaVal)  where ebno
        % is  the value of ebno to use. nbSim simulations are run for 
        % values of ia equal to iaVal.
        
            % Arg check
            switch length(varargin)
                case 0
                    nbSim   = 10;
                    ia   = 0:0.05:1;
                case 1
                    nbSim   = varargin{1};
                    ia   = 0:0.1:1;
                case 2
                    nbSim   = varargin{1};
                    ia   = varargin{2};
                otherwise
                    error('Too many input argument')
            end
            if nbSim < 1
                error('There must be at least 1 EXIT simulation');
            end
            if ia < 0
                error('Property iaVal: Ia values must be positive');
            end
            if length(ebno) ~= 1
                error('EbN0 must be a scalar')
            end
            if obj.blkLength<1e4
                warning('This size of blkLength (%d) may lead to imprecise results', obj.blkLength)
            end
            
            % Initialize
            noiseVar    = 1./(obj.rate*10.^(ebno/10)*2);        % Noise variance
            
            % apriori variance values for exit computation
            % reference values for a priori
            refAprioriVar  = [0 linspace(1e-2, 0.9, 6) logspace(0, 2, 200) 110:10:250];  
            % Compute Ia for all those values
            iaJ = zeros(1, length(refAprioriVar));
            for i = 1:length(refAprioriVar)
                iaJ(i) = obj.computeMI('J', sqrt(refAprioriVar(i)));
            end
            % Avoid overflow for last value (infinite variance)
            ia(end) = ia(end)-1e-4;
            % Interpolate the variance for each Ia values set as parameter
            aprioriVar = interp1(iaJ, refAprioriVar, ia, 'linear');
            ie = zeros(length(aprioriVar), 1);
            
            % Computations
            for idxVar = 1:length(aprioriVar)
                disp(['Ia = ' num2str(ia(idxVar))])
                for idxSim = 1:nbSim
                    % Extrinsic mutual information
                    ie(idxVar) = ie(idxVar) + obj.computeExMI(noiseVar, aprioriVar(idxVar))/nbSim;
                end
                if ie(idxVar) < ia(idxVar) && false   % Use to save computation time
                    disp('EXIT cuts the diagonal line')
                    ie((idxVar+1):end) = NaN;
                    break
                end
            end
            if obj.backup
                % Save data script
                obj.saveData('exit', ebno, ia, ie)
            end
        end
        
        function [ia, ie] = computeTraj(obj, ebno, varargin)
            %computeTraj Compute trajectory for the given ebno value 
            % and using the parameters of the object.
            %
            % [ia, ie] = computeTraj(obj, ebno)  where ebno is the value of
            % ebno to use. 1 (default value) trajectory is computed.
            %
            % [ia, ie] = computeTraj(obj, ebno, nbTraj)  where ebno is the 
            % value of ebno to use. nbTraj trajectories are computed.
            %
            % [ia, ie] = computeTraj(obj, ebno, nbTraj, true)  where ebno 
            % is the value of ebno to use. An average over nbTraj 
            % trajectories is computed.
            
            % Arg check
            switch length(varargin)
                case 0
                    nbTraj      = 1;
                    averageTraj = false;
                case 1
                    nbTraj      = varargin{1};
                    averageTraj = false;
                case 2
                    nbTraj      = varargin{1};
                    averageTraj = varargin{2};
                otherwise
                    error('Too many input argument')
            end
            if nbTraj < 1
                error('There must be at least 1 trajectory');
            end
            if ~islogical(averageTraj)
                error('argument "averageTraj" must be logical');
            end
            if length(ebno) ~= 1
                error('EbN0 must be a scalar')
            end
        
            % Initialization variables
            ie      = zeros(2, obj.numIter, nbTraj);
            ia      = zeros(size(ie));
            errors  = zeros(obj.numIter, nbTraj);
            
            % Noise variance
            noiseVar    = 1./(obj.rate*10.^(ebno/10)*2);        

            % Computations
            for idxTraj = 1:nbTraj
                disp(['Trajectory ' num2str(idxTraj) '/' num2str(nbTraj)])
                [errors(:,idxTraj), ia(:,:,idxTraj), ie(:,:,idxTraj)] = obj.computeChain(noiseVar, true);
            end
            if averageTraj
                ie = mean(ie, 3);
                ia = mean(ia, 3);
            end
        end
    
        function h = plotBer(obj, ebno, ber, varargin)
            %plotBer Generate a figure with BER results
            %
            % h = plotBer(obj, ebno, ber) where ebno is the ebno vector for
            % the simulation and ber the BER results. The result for the
            % last iteration (default value) is plotted, versus the Eb/N0
            % (default value).
            %
            % h = plotBer(obj, ebno, ber, 'snr') where ebno is the ebno 
            % vector for the simulation and ber the BER results. The result 
            % for the last iteration (default value) is plotted, versus the
            % SNR.
            %
            % h = plotBer(obj, ebno, ber, 'snr', iter) where ebno is the 
            % ebno vector for the simulation and ber the BER results. The 
            % result for the iterations iter is plotted, versus the
            % SNR.
            
            % Check args
            switch length(varargin)
                case 0
                    plotType    = 'ebno';
                    iter        = obj.numIter;
                case 1
                    plotType    = varargin{1};
                    iter        = obj.numIter;
                case 2 
                    plotType    = varargin{1};
                    iter        = varargin{2};
                otherwise
                    error('Too many input arguments')
            end
            switch plotType
                case 'ebno'
                    dispSnr = 0;
                case 'snr'
                    dispSnr = 1;
                otherwise
                    error('Parameter "plotType" can only be "ebno" or "snr"')
            end
            if iter<1
                error('Iteration(s) number(s) must be a strictly positive value')
            end
            % Remove iteration above numIter
            iter(iter>obj.numIter) = [];
            
            % Legend
            if size(iter,1)>1
                iter = iter';
            end
            legendStr = strcat({'Iteration '}, num2str(iter));
            % Plot
            h = figure('Name', 'BER');
            semilogy(ebno + dispSnr*10*log10(obj.rate), ber(:,iter))
            ylabel('BER')
            if dispSnr
                xlabel('SNR (dB)')
            else
                xlabel('Eb/N0 (dB)')
            end
            grid on
            legend(legendStr)
            title([{'BER'};{obj.idStr}])
            
        end
        
        function h = plotTraj(obj, ebno, ia, ie)
            %plotTraj Generate a figure displaying EXIT chart and trajectories
            
            % Check if corresponding EXIT has been computed (whatever the
            % size)
            exitFile = ['saves/exit_' obj.saveStr '_ebno_' num2str(ebno*100) 'e-2.mat'];
            exitFile = strrep(exitFile, ['N' num2str(obj.blkLength)], 'N*');
            files = dir(exitFile);
            if isempty(files)
                disp(['No EXIT chart found for this configuration'])
                disp(['Computing EXIT Chart (blkLength=100000)'])
                tc = obj;
                tc.backup = true;
                tc.blkLength = 100000;
                tc.computeExit(ebno);
                exitFile = ['saves/exit_' tc.saveStr '_ebno_' num2str(ebno*100) 'e-2.mat'];
                loadExit = load(exitFile);
            else
                % Search
                loadExit = load(files(1).name);
            end
            
            % Plot preparation
            nbTraj = size(ia, 3);
            trajx = zeros(2*size(ia, 2), nbTraj);
            trajy = zeros(size(trajx));

            for idxTraj = 1:nbTraj
                trajxtmp = [ia(1, :, idxTraj); ie(2, :, idxTraj)];
                trajx(:, idxTraj) = trajxtmp(:)';

                trajytmp = [ia(2, :, idxTraj); ie(1, :, idxTraj)];
                trajy(:, idxTraj) = trajytmp(:)';
            end
            
            % Plot
            h = figure('Name', 'Trajectories');
            plot(loadExit.ia, loadExit.ie, 'k')
            hold on
            plot(loadExit.ie, loadExit.ia, 'k')
            plot(trajx, trajy)
            xlabel('I_A')
            ylabel('I_E')
            axis([0 1 0 1])
            title([{'Trajectories'};{obj.idStr}])
            
        end    
        
        function h = plotExit(obj, ebno, ia, ie)
            %plotExit Generate a figure with EXIT chart results
        
            % Plot
            h = figure('Name', 'EXIT Chart');
            plot(ia, ie, 'k')
            xlabel('I_A')
            ylabel('I_E')
            axis([0 1 0 1])
            title([{['Exit Chart for E_b/N_0 = ' num2str(ebno) 'dB']};{obj.idStr}])
        end
        
    end
    
    methods (Access = private)
        function saveData(obj, dataType, varargin)
            % saveData Save special data and object parameters to a file
            %
            % obj.saveData('ber', ebno, ber, per) saves BER and PER data,    
            % along with EbN0 values and parameters of the object.
            %
            % obj.saveData('exit', ebno, ia, ie) saves IA and IE data,    
            % along with the EbN0 value and parameters of the object.
            
            % Get User param in a struct, to be saved
            Param.dataType      = dataType;
            Param.blkLength     = obj.blkLength;
            Param.k             = obj.k;
            Param.poly          = obj.poly;
            Param.puncturing    = obj.puncturing;
            Param.algorithm     = obj.algorithm;
            Param.numIter       = obj.numIter;
            Param.rate          = obj.rate;
            Param.idStr         = obj.idStr;
            
            switch dataType
                case 'ber'
                    if length(varargin) ~= 3
                        error('Not enough input arguments')
                    else
                        Param.ebno  = varargin{1};
                        ber         = varargin{2};
                        per         = varargin{3};
                    end
                    filename    = ['saves/ber_' obj.saveStr];
                    save(filename, 'Param', 'ber', 'per');
                
                case 'exit'
                    if length(varargin) ~= 3
                        error('Not enough input arguments')
                    else
                        Param.ebno  = varargin{1};
                        ia          = varargin{2};
                        ie          = varargin{3};
                    end
                    filename    = ['saves/exit_' obj.saveStr '_ebno_' num2str(Param.ebno*100) 'e-2'];
                    save(filename, 'Param', 'ie', 'ia');
                
                otherwise
                    error('Unknown argument for "dataType"');
            end
        end
        
        function obj = getRSCTrellis(obj)
            %getRSCTrellis Compute elements of the trellis using some parameters
            %   
            %   obj.getRSCTrellis( ) uses obj.k and obj.poly to generate 
            %   obj.codewords and obj.nextStates.
            %   
            %   obj.k is the constraint length (equals to 1 + number of 
            %   memories) 
            % 
            %   obj.poly is a 2 elements vector defined as [Feedback, Generator] 
            %   (where Feedback is the feedback polynomial (in octal) and 
            %   Generator the generator polynomial (in octal) )
            %   
            %   obj.nextStates is a 2x2^(obj.k-1) matrix, where row index is the 
            %   input (0 or 1) and column index is the initial state
            %
            %   obj.codewords gives the outputed codewords. It is a 2x2^k 
            %   matrix, where the 2^(k-1) first columns give the codewords 
            %   if the input is 0 and initial state is the column index, 
            %   while the other 2^(k-1) columns give the codewords if the 
            %   input is 1 and initial state is (column index - 2^(k-1))

            % Errors check
            if length(obj.poly)>2
                error('Feedback and Generator polynomials must be specified; only rate 1/2 are supported')
            end
            for i = 1:length(obj.poly)
                if base2dec(num2str(obj.poly(i)), 8) > (2^obj.k-1)
                    error('Generator polynomial [%d] does not match the constraint length k = %d', obj.poly(i), obj.k)
                end
            end

            % Outputs
            obj.nextStates = zeros(2,2^(obj.k-1));
            obj.codewords = zeros(2,2^obj.k);

            % Convert poly into binary
            Feedback = myde2bi(base2dec(num2str(obj.poly(1)), 8), obj.k, 'left-msb')';
            Generator = myde2bi(base2dec(num2str(obj.poly(2)), 8), obj.k, 'left-msb')';

            % Check if k matches the polynomials
            if Feedback(1) == 0 && Generator(1) == 0
                error('k value (%d) is oversized for polynomials [%d %d]', obj.k, obj.poly(1), obj.poly(2))
            end

            % For each state
            for state = 1:2^(obj.k-1)
                % For each possible inputs
                for input = [0 1]
                    % Set memory to the current memory
                    mem = myde2bi(state-1, obj.k-1, 'left-msb');
                    % Feedback polynomial output
                    outFeedback = mod([input, mem]*Feedback, 2);
                    % Generator polynomial output
                    output = mod([outFeedback, mem]*Generator, 2);
                    % Update memory units
                    mem(2:end) = mem(1:end-1);
                    mem(1) = outFeedback;
                    % Outputed codeword
                    obj.codewords(:, state + input*2^(obj.k-1)) = [input; output];
                    % Next state
                    obj.nextStates(input+1, state) = mybi2de(mem, 'left-msb')+1;
                end

            end
        end
        
        function [outputBits, terminationBits] = rscEncoder(obj, inputBits)
            %RSCencoder Encode input stream with RSC constituent
            %   
            %   [ outputBits, terminationBits ] = obj.rscEncoder( inputBits )
            %   where inputBits are the bits to be encoded

            % Errors check
            if length(obj.poly)>2
                error('Feedback and Generator polynomials must be specified; only rate 1/2 are supported')
            end
            for i = 1:length(obj.poly)
                if base2dec(num2str(obj.poly(i)), 8) > (2^obj.k-1)
                    error('Generator polynomial [%d] does not match the constraint length k=%d', obj.poly(i), obj.k)
                end
            end

            % Output vector
            outputBits = zeros(1,length(inputBits)+(obj.k-1));
            terminationBits = zeros(1,obj.k-1);

            % Initialize memory units
            mem = zeros(1, obj.k-1);

            % Convert Poly into binary
            Feedback = myde2bi(base2dec(num2str(obj.poly(1)), 8), obj.k, 'left-msb')';
            Generator = myde2bi(base2dec(num2str(obj.poly(2)), 8), obj.k, 'left-msb')';
            
            
            % Find output bits
            for i = 1:length(inputBits)
                % Feedback polynomial output
                outFeedback = mod([inputBits(i), mem]*Feedback, 2);
                % Generator polynomial output
                outputBits(i) = mod([outFeedback, mem]*Generator, 2);
                % Update memory units
                mem(2:end) = mem(1:end-1);
                mem(1) = outFeedback;
            end

            % Find termation bits needed to go back to state 0
            for i = 1:(obj.k-1)
                % Find termination bits to force feedback to 0
                terminationBits(i) = mod(mem*Feedback(2:end), 2);
                % Generator polynomial output
                outputBits(length(inputBits)+i) = mod([0, mem]*Generator, 2);
                % Update memory units
                mem(2:end) = mem(1:end-1);
                mem(1) = 0;
            end
        end
        
        function [output] = encoder(obj, data)
            %ENCODER Perform turbo encoding of the input stream
            % 
            %   [ outputBits ] = obj.encoder ( binaryData ) where
            %   binaryData is the data to be encoded
            
            % First Encoder (no interleaving)
            [out1, end1] = obj.rscEncoder(data);

            % Second Encoder
            [out2, end2] = obj.rscEncoder(data(obj.interleaverIdx));

            % Encoded bits
            output = [data end1; data end2; out1; out2];
            % Note : data is sent twice but this a no effect since variance is
            % computed with the actual code rate; however only one should be 
            % used as systematic input for the decoder
        end
        
        function logAPP = bcjr(obj, llr)
            %BCJR Apply BCJR algorithm on the input stream
            %
            %   [ logAPP ] = obj.bcjr(llr) 

            % Variable
            nStates     = 2^(obj.k-1);              % Number of states
            N           = size(llr,2);          % Number of transitions/codewords


            switch(obj.algorithm)
                case 'MAP'
                    % gamma : quantity linked to the transition probabilities
                    gamma = exp(llr'*(1-2*(obj.codewords))/2);

                    % Norm 
                    gamma = gamma./repmat(sum(gamma, 2), 1, nStates*2);

                    % Execute forward & Backward computations with Cmex function
                    [alpha, beta] = ForwardBackward(gamma', int32((obj.nextStates-1)'), int32((obj.prevStates-1)'), nStates, N, 0);

                    % Update transition (ie codewords) probabilities with alpha & beta
                    P = zeros(size(gamma));
                    for trans = 1:nStates
                        P(:,trans) = alpha(trans,1:end-1)'.*gamma(:,trans).*beta(obj.nextStates(1, trans),2:end)';
                        P(:,nStates+trans) = alpha(trans,1:end-1)'.*gamma(:,nStates+trans).*beta(obj.nextStates(2, trans),2:end)';
                    end
                    P = P./repmat(sum(P, 2), 1, nStates*2);

                    % A posteriori llr
                    logAPP = (log(sum(P(:,1:nStates), 2) + 10^(-44)) - log(sum(P(:,(nStates+1):end), 2) + 10^(-44)))';

                case 'maxlog'
                    % gamma : quantity linked to the transition probabilities
                    gamma = llr'*(1-2*(obj.codewords))/2;

                    % Execute forward & Backward computations with Cmex function
                    [alpha, beta] = ForwardBackward(gamma', int32((obj.nextStates-1)'), int32((obj.prevStates-1)'), nStates, N, 1);

                    % Update transition (ie codewords) probabilities with alpha & beta
                    P = zeros(size(gamma));
                    for trans = 1:nStates
                        P(:,trans) = alpha(trans,1:end-1)' + gamma(:,trans) + beta(obj.nextStates(1, trans),2:end)';
                        P(:,nStates+trans) = alpha(trans,1:end-1)' + gamma(:,nStates+trans) + beta(obj.nextStates(2, trans),2:end)';
                    end

                    % Log ratio des bits (using max log rule)
                    logAPP = (max(P(:,1:nStates), [], 2) - max(P(:,(nStates+1):end), [], 2))';

            %         keyboard
                otherwise
                    error('Unknown algorithm for BCJR')
            end

            % Remove termination bits
            logAPP = logAPP(1:end-(obj.k-1));

        end
        
        function [errors, ia, ie] = tcDecoder(obj, llr, data, withTraj)
            %TCDECODER Perform turbo decoding (optionally, with trajectory)
            %
            %   [ errors ] = obj.tcDecoder( llr, data, false ) performs turbo
            %   decoding and compute the number of errors after obj.numIter
            %   iterations
            %   
            %   [ errors, ia, ie ] = obj.tcDecoder( llr, data, true )
            %   performs turbo decoding and compute the number of errors
            %   after obj.numIter iterations along with the exchanges of MI
            %   between the two constituents
            
            % Output vectors : errors for each iteration
            errors      = zeros(1, obj.numIter);
            if withTraj
                ie          = zeros(2, obj.numIter);
                ia          = zeros(size(ie));
                countStuck  = 0;
            else
                ie          = [];
                ia          = [];
            end
                

            % Reshape the LLR for decoding (one "layer" for each decoder)
            llrCh = zeros(2, obj.blkLength + (obj.k-1), 2);
            % Systematic bits
            llrCh(1, 1:obj.blkLength, 1) = llr(1, 1:obj.blkLength);
            llrCh(1, 1:obj.blkLength, 2) = llr(1, obj.interleaverIdx);
            % Termination bits
            llrCh(1, (obj.blkLength+1):end, 1) = llr(1, (obj.blkLength+1):end);
            llrCh(1, (obj.blkLength+1):end, 2) = llr(2, (obj.blkLength+1):end);
            % Parity bits
            llrCh(2, :, 1) = llr(3, :);
            llrCh(2, :, 2) = llr(4, :);

            % Systematic from channel
            sysCh = llr(1, 1:obj.blkLength);

            % Extrinsic vectors
            E1 = zeros(1, obj.blkLength);
            E2 = zeros(size(E1));

            %% Iterative decoding        
            iter = 1;
            while iter<=obj.numIter
                %%%%% First Decoder
                % Build the matrix containing channel + a priori
                A = llrCh(:,:,1);
                A(1,1:obj.blkLength) = sysCh + E2;

                % Compute the APP
                APP = obj.bcjr(A);

                % Compute the extrinsic
                E1 = APP - E2 - sysCh;

                %%%%% Second Decoder
                % Build the vector containing channel + a priori
                A = llrCh(:,:,2);
                A(1, 1:obj.blkLength) = sysCh(obj.interleaverIdx) + E1(obj.interleaverIdx);

                % Compute the APP
                APP(obj.interleaverIdx) = obj.bcjr(A);

                % Compute the extrinsic
                E2 = APP - E1 - sysCh;

                % Compute number of errors
                errors(iter) = sum(not((APP<0)==data));
                
                if withTraj
                    % Compute Ia of the first decoder
                    if iter == 1
                        ia(1, iter) = 0;
                    else
                        ia(1, iter) = ie(2, iter-1);
                    end
                    % Compute Ie of the first decoder
                    ie(1, iter) = obj.computeMI('histogram', E1, 1-2*data);
                    % Compute Ia of the second decoder
                    ia(2, iter) = ie(1, iter);
                    % Compute Ie of the second decoder
                    ie(2, iter) = obj.computeMI('histogram', E2, 1-2*data);
                    % Display error
                    disp(['iter ' num2str(iter) ' : ' num2str(errors(iter)) ' errors'])
                end
                    
                % If no error then stop iterations
                if errors(iter) == 0
                    if withTraj
                        ie(:, (iter+1:end)) = 1;
                        ia(:, (iter+1:end)) = 1;
                    end
                    break
                end

                if withTraj
                    % Search if decoder has converged to a point different of (1, 1)
                    if iter>1
                        if errors(iter) == errors(iter-1)
                            countStuck = countStuck +1;
                        else 
                            countStuck = 0;
                        end
                        % If stucked for 10 iterations, stop simulation
                        if countStuck == 10
                            ie(1, (iter+1:end)) = ie(1, iter);
                            ie(2, (iter+1:end)) = ie(2, iter);
                            ia(1, (iter+1:end)) = ia(1, iter);
                            ia(2, (iter+1:end)) = ia(2, iter);
                            disp('Decoder is stucked: stopping simulation')
                            break
                        end
                    end
                end
                
                % Next iteration
                iter = iter+1;
            end
            if withTraj && iter==(obj.numIter+1)
                disp('Maximum number of iterations has been reached')
            end
        end
        
        function [errors, ia, ie] = computeChain(obj, noiseVar, withTraj)
            %computeChain Performs operations of a transmission chain
            
            % Generate random bits
            data = randi([0 1], 1, obj.blkLength);

            % Encoded bits
            bEnc = obj.encoder(data);

            % BPSK
            xEnc = 1-2*bEnc; 

            % AWGN Channel
            rData       = obj.awgnChannel(xEnc, noiseVar);

            % LLR
            llrData     = real(rData)*2/noiseVar;

            if obj.puncturing
                % Set every other parity LLR to 0
                llrData(3, 1:2:end) = 0;
                llrData(4, 2:2:end) = 0;
            end

            % Turbo decoding
            [errors, ia, ie] = obj.tcDecoder(llrData, data, withTraj);
        end
        
        function ie = computeExMI(obj, noiseVar, aprioriVar)
            %COMPUTEEXMI Compute the extrinsic MI of one constituent for a 
            % given a priori variance
            
            % Information bits
            bits = randi([0 1], 1, obj.blkLength);

            % Encoding
            [dataEnc, dataEncTerm] = obj.rscEncoder(bits);

            % Output encoded bits
            bEnc = [bits dataEncTerm; dataEnc];

            % BPSK
            x = 1-2*bEnc;

            % Apriori
            llra = (aprioriVar/2)*(1-2*bits) + sqrt(aprioriVar)*randn(1, obj.blkLength);

            % AWGN Channel
            y = obj.awgnChannel(x, noiseVar);

            % BPSK LLR
            llrCh = real(y)*2/noiseVar;
            if obj.puncturing
                % Set every other parity LLR to 0
                llrCh(2, 1:2:end) = 0;
            end

            % Systematic bits from channel
            sysCh = llrCh(1, 1:obj.blkLength);

            % Decoding
            % Input of decoder
            llrIn = llrCh;
            llrIn(1, 1:obj.blkLength) = sysCh + llra;

            % Compute the APP
            llr = obj.bcjr(llrIn);

            % Extrinsic LLR
            llrEx = llr - sysCh - llra;
            
            % Compute mutual information
            ie = obj.computeMI('histogram', llrEx, 1-2*bits);
        end
        
    end
    
    methods (Static, Access = private)
        function [output] = awgnChannel(data, noiseVar)
            %AWGNCHANNEL Add AWGN noise to the input, with specified variance
            noise       = sqrt(noiseVar)*randn(size(data)); 
            output      = data + noise;
        end
        
        function [ I ] = computeMI(method, varargin)
            %COMPUTEMI Compute mutual information usign various methods
            %   
            %   I = computeMI('histogram', llr, bits) compute mutual information using
            %   the histogram method. llr are the log likelihood ratios of the bits and
            %   bits the real value (equals to -1 or 1)
            %
            %   I = computeMI('J', sigma) compute mutual information using the J 
            %   function from ten Brink paper. sigma is the standard deviation.
            %
            %   I = computeMI('averaging', llr, bits) compute mutual information using
            %   the Hagenauer method. llr are the log likelihood ratios of the bits and
            %   bits the real value (equals to -1 or 1)

            % Argument check
            if length(varargin) > 2
                error('Too many input arguments')
            elseif isempty(varargin)
                error('Not Enough input arguments')
            end

            switch method
                case 'histogram'
                    % Check arguments
                    if length(varargin) ~= 2
                        error('Not Enough input arguments for histogram computation')
                    else
                        llr     = varargin{1};
                        bits    = varargin{2};
                    end
                    % Histogram method (ten Brink, real entropy computation)
                    minDistrib = min(llr);
                    maxDistrib = max(llr);

                    % Check if both types of llr are represented
                    if isempty(llr(bits==1)) || isempty(llr(bits==-1))
                        I = 0;
                        return
                    end

                    % Prepare to compute histogram
                    bins = linspace(minDistrib, maxDistrib, 200);
                    deltaBins = abs(bins(2)-bins(1));

                    % Compute histograms
                    pdf0 = histc(llr(bits==1), bins);
                    pdf0 = pdf0/(sum(pdf0)*deltaBins);

                    pdf1 = histc(llr(bits==-1), bins);
                    pdf1 = pdf1/(sum(pdf1)*deltaBins);

                    % Find not 0 values
                    val0 = find(pdf0~=0);
                    val1 = find(pdf1~=0);

                    % Compute the terms to be integrated
                    part0 = pdf0(val0).*log2(2*pdf0(val0)./(pdf0(val0)+pdf1(val0)));
                    part1 = pdf1(val1).*log2(2*pdf1(val1)./(pdf0(val1)+pdf1(val1)));

            %         figure; plot(bins, pdf0); hold on; plot(bins, pdf1, 'r'); legend('P(y|X = 0)', 'P(y|X = 1)')
                    % Compute mutual information
                    I = 0.5*sum(part0)*deltaBins + 0.5*sum(part1)*deltaBins;

                case 'J'
                    % Check arguments
                    if length(varargin) ~= 1
                        error('Too many input arguments for J function computation')
                    else
                        aprioriVar = varargin{1}^2;
                    end
                    % J function
                    if aprioriVar > 0
                        gaussFun = @(x) exp(-(x-aprioriVar/2).^2/(2*aprioriVar))/(sqrt(2*pi*aprioriVar)).*log2(1+exp(-x));
                        I = 1 - integral(gaussFun, -200, 200);
                    else
                        I = 0;
                    end

                case 'averaging'
                    % Check arguments
                    if length(varargin) ~= 2
                        error('Not Enough input arguments for averaging')
                    else
                        llr     = varargin{1};
                        bits    = varargin{2};
                    end
                    % Averaging method (Hagenauer method)... not sure it works
                    I = mean(1 - log2(1 + exp(-bits.*llr)));

                otherwise
                    error('Unknown method')
            end
        end
    end
            
end

