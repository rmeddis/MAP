%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef cHMM
    %HMMCLASS Handles all of the HTK related gubbins
    %   Please see the documentation located in a separate file for further
    %   information.     
    
    %%  *********************************************************
    %  properties                      _   _
    %                                 | | (_)
    %  _ __  _ __ ___  _ __   ___ _ __| |_ _  ___  ___
    % | '_ \| '__/ _ \| '_ \ / _ \ '__| __| |/ _ \/ __|
    % | |_) | | | (_) | |_) |  __/ |  | |_| |  __/\__ \
    % | .__/|_|  \___/| .__/ \___|_|   \__|_|\___||___/
    % | |             | |
    % |_|             |_|
    %************************************************************
    
    %% **********************************************************
    % Public properties - can be set by user
    %************************************************************
    properties(Access = public)
        hmmFolder
        paramType               = 'USER_D_A'; %DELTAS and ACCELERATIONS
        numCoeff                = 27; %9*3 THIS IS FOR PROBABILITY MODEL (not high spont+low spont which would be 18*3=54)
                
        HERestDataPath      = fullfile(pwd, 'def', 'HERest_digit');
        binPath             = fullfile(pwd, 'def', 'bin');
        configFile          = fullfile(pwd, 'def', 'config_STANDARD');
        trainWordListFile   = fullfile(pwd, 'def', 'Grammar_digit', 'words3');
        testWordListFile    = fullfile(pwd, 'def', 'Grammar_digit', 'wordsNoSevenZero');
        wordNetFile         = fullfile(pwd, 'def', 'Grammar_digit', 'wdnetNoSP.slf');
        dictFile            = fullfile(pwd, 'def', 'Grammar_digit', 'noSevenZeroDict');
    end
    
    %% **********************************************************
    % Dependent - never set by user. Only calculated when needed
    %************************************************************
    properties(Dependent = true)
        protoFile              % = fullfile(pwd, 'def', 'proto_RobANonly_9'); %probability only
    end
    
    %%  *********************************************************
    % methods        _   _               _
    %               | | | |             | |
    % _ __ ___   ___| |_| |__   ___   __| |___
    %| '_ ` _ \ / _ \ __| '_ \ / _ \ / _` / __|
    %| | | | | |  __/ |_| | | | (_) | (_| \__ \
    %|_| |_| |_|\___|\__|_| |_|\___/ \__,_|___/
    %************************************************************
    
    methods
        %% **********************************************************
        % Constructor
        %************************************************************
        function obj = cHMM(hmmFolder)
            if nargin > 0
                obj.hmmFolder = hmmFolder;
            end
        end % ------ OF CONSTRUCTOR
        
        %% **********************************************************
        % genProto - generate task specific prototype
        %************************************************************
        function genProto(obj)
            % models_1mixsil.exe - takes input (hmmdef) and copies it making a
            % basis of one, two three etc. etc.
            obj.models_1mixsilMat(fullfile(obj.hmmFolder,'hmm0','hmmdef'), fullfile(obj.hmmFolder,'hmm0','models'));            
        end % ------ OF GENPROTO
        
        %% **********************************************************
        % istrained
        %************************************************************
        function boolans = istrained(obj)
            boolans = numel(dir(fullfile(obj.hmmFolder,'hmm36','models')));
        end
        
        %% **********************************************************
        % Train
        %************************************************************
        function train(obj, trainFeatureFolder)
                       
            % Most of the following code block can be replaced by doing a
            % find and replace across the  code below. I didn't
            % want to mess with the working code too much, so I just
            % copied the object properties needed into the variable names
            % expected.
            
            ED_CMDFILE1 = fullfile(obj.HERestDataPath, 'sil1.hed');
            ED_CMDFILE2 = fullfile(obj.HERestDataPath, 'mix2_16.hed');
            ED_CMDFILE3 = fullfile(obj.HERestDataPath, 'mix3_16.hed');
            ED_CMDFILE4 = fullfile(obj.HERestDataPath, 'mix5_16.hed');
            ED_CMDFILE5 = fullfile(obj.HERestDataPath, 'mix7_16.hed');
            
            NUM_COEF        = obj.numCoeff;
            PAR_TYPE        = obj.paramType;
            LIST_FILE       = fullfile(obj.hmmFolder, 'tmp.list');
            
            word_list       = obj.trainWordListFile;%fullfile(obj.grammarPath, 'words');
            word_listSP     = word_list; % for use in hmm4 onwards - UGLY HACK NOW SP ABANDONED
            proto           = obj.protoFile;%Does not exist on disk just yet probably - see a few lines down
            config          = obj.configFile;
            train_list      = fullfile(trainFeatureFolder, 'list.scp');
            labels          = fullfile(trainFeatureFolder, 'labels.mlf');
            labelssp        = labels; % for use in hmm4 onwards - UGLY HACK NOW SP ABANDONED
            hmm_dir         = obj.hmmFolder;
            
            FEAT_ROOT       = trainFeatureFolder;
            
            % Now for the actual HMM training code
            mkdir(hmm_dir)
            for I = 0:36
                h = fullfile(hmm_dir,['hmm' num2str(I)]);
                mkdir(h);
            end
            obj.makeProtoHmm(proto, obj.paramType, obj.numCoeff, 18);
            
            fid = fopen(train_list,'r');
            disp(train_list)
            S = textscan(fid,'%s','Delimiter','\n');
            fclose(fid);
            
            fid = fopen(LIST_FILE,'w');
            
            for I = 1:size(S{1},1)
                str  = fullfile(FEAT_ROOT,S{1}{I});
                fprintf(fid,'%s\n',str);
            end;
            fclose(fid);
            
            % HCompV just gets the vfloor stuff out so we can begin approximating
            cmd = ['"HCompV" -T 2 -D -C "' config '" -o hmmdef -f 0.01 -m -S "' LIST_FILE '" -M "' hmm_dir filesep 'hmm0" "' proto '"'];
            system(cmd);                        
            
            %cmd = [BINDIR filesep 'macro' binExt ' ' num2str(NUM_COEF) ' ' PAR_TYPE ' ' '"' hmm_dir filesep 'hmm0' filesep 'vFloors' '" "' hmm_dir filesep 'hmm0' filesep 'macros' '"'];
            %system(cmd);            
            obj.macroMat(NUM_COEF,PAR_TYPE, fullfile(hmm_dir, 'hmm0', 'vFloors'), fullfile(hmm_dir, 'hmm0', 'macros'));
            
            %MAKE THE INITIAL MODEL PROTOTYPE
            genProto(obj);
                        
            disp('Seed HMM successfully Produced.....');
            
            %Training
            for I = 1:3
                disp(I)
                j = I-1;
                cmd = ['HERest -D -C ' config ' -I ' labels ' -t 250.0 150.0 1000.0 -S ' LIST_FILE  ' -H ' hmm_dir filesep 'hmm' num2str(j)  filesep 'macros -H ' hmm_dir filesep 'hmm' num2str(j) filesep 'models -M ' hmm_dir filesep 'hmm' num2str(I) ' ' word_list];
                %                 disp(cmd)
                system(cmd);
            end
            
            disp('3 iterations complete');
            
            rmdir ([hmm_dir filesep 'hmm4'],'s')
            copyfile ([hmm_dir filesep 'hmm3'], [hmm_dir filesep 'hmm4'])                                                
            
            % The following command takes state 3 from the silence model
            % and appends it to the end of the model as state 2 of the
            % short pause model.
            % Original:
            % cmd = [BINDIR filesep 'spmodel_gen' binExt ' ' hmm_dir filesep 'hmm3' filesep 'models ' hmm_dir filesep 'hmm4' filesep 'models'];
            % system(cmd);
            % New:
            obj.spmodel_genMat(fullfile(hmm_dir,'hmm3','models'), fullfile(hmm_dir,'hmm4','models'));                                     
            
            cmd = ['HHEd  -T 2 -H ' hmm_dir filesep 'hmm4' filesep 'macros -H ' hmm_dir filesep 'hmm4' filesep 'models -M ' hmm_dir filesep 'hmm5 ' ED_CMDFILE1 ' ' word_listSP ];
            system(cmd);
            disp ('SP model fixed')
            
            % after the spmodel_gen command - the word_list is changed to
            % word_listSP. The sp model is just ignored currently
            
            for I = 6:8
                disp(I)
                j = I-1;
                cmd = ['HERest -C ' config ' -I ' labelssp ' -S ' LIST_FILE  ' -H ' hmm_dir filesep 'hmm' num2str(j)  filesep 'macros -H ' hmm_dir filesep 'hmm' num2str(j) filesep 'models -M ' hmm_dir filesep 'hmm' num2str(I) ' ' word_listSP];
                system(cmd);
            end
            
            disp('6 iterations complete');
            
            cmd = ['HHEd  -T 2 -H ' hmm_dir filesep 'hmm8' filesep 'macros -H ' hmm_dir filesep 'hmm8' filesep 'models -M ' hmm_dir filesep 'hmm9 ' ED_CMDFILE2 ' ' word_listSP ];
            system(cmd);
            disp ('2 gaussians per mixture')
            
            for I = 10:12
                disp(I)
                j = I-1;
                cmd = ['HERest -C ' config ' -I ' labelssp ' -S ' LIST_FILE  ' -H ' hmm_dir filesep 'hmm' num2str(j)  filesep 'macros -H ' hmm_dir filesep 'hmm' num2str(j) filesep 'models -M ' hmm_dir filesep 'hmm' num2str(I) ' ' word_listSP];
                system(cmd);
            end
            
            disp ('9 iterations completed')
            
            cmd = ['HHEd  -T 2 -H ' hmm_dir filesep 'hmm12'  filesep 'macros -H ' hmm_dir filesep 'hmm12' filesep 'models -M ' hmm_dir filesep 'hmm13 ' ED_CMDFILE3 ' ' word_listSP ];
            system(cmd);
            disp ('3 gaussians per mixture')
            
            for I = 14:20
                disp(I)
                j = I-1;
                cmd = ['HERest -C ' config ' -I ' labelssp ' -S ' LIST_FILE  ' -H ' hmm_dir filesep 'hmm' num2str(j)  filesep 'macros -H ' hmm_dir filesep 'hmm' num2str(j) filesep 'models -M ' hmm_dir filesep 'hmm' num2str(I) ' ' word_listSP];
                system(cmd);
            end
            
            disp ('16 iterations completed')
            
            cmd = ['HHEd  -T 2 -H ' hmm_dir filesep 'hmm20'  filesep 'macros -H ' hmm_dir filesep 'hmm20' filesep 'models -M ' hmm_dir filesep 'hmm21 ' ED_CMDFILE4 ' ' word_listSP ];
            system(cmd);
            disp ('5 gaussians per mixture')
            
            for I = 22:28
                disp(I)
                j = I-1;
                cmd = ['HERest -C ' config ' -I ' labelssp ' -S ' LIST_FILE  ' -H ' hmm_dir filesep 'hmm' num2str(j)  filesep 'macros -H ' hmm_dir filesep 'hmm' num2str(j) filesep 'models -M ' hmm_dir filesep 'hmm' num2str(I) ' ' word_listSP];
                system(cmd);
            end
            
            disp ('23 iterations completed')
            
            cmd = ['HHEd  -T 2 -H ' hmm_dir filesep 'hmm28'  filesep 'macros -H ' hmm_dir filesep 'hmm28' filesep 'models -M ' hmm_dir filesep 'hmm29 ' ED_CMDFILE5 ' ' word_listSP ];
            system(cmd);
            disp ('7 gaussians per mixture')
            
            for I = 30:36
                disp(I)
                j = I-1;
                cmd = ['HERest -C ' config ' -I ' labelssp ' -S ' LIST_FILE  ' -H ' hmm_dir filesep 'hmm' num2str(j)  filesep 'macros -H ' hmm_dir filesep 'hmm' num2str(j) filesep 'models -M ' hmm_dir filesep 'hmm' num2str(I) ' ' word_listSP];
                system(cmd);
            end
            
            disp ('30 iterations completed')
            
        end % ------ OF TRAIN
        
        %% **********************************************************
        % Test
        %************************************************************
        function test(obj, testFeatureFolder)
            flags = '-p 0.0 -s 0.0';
            test_word_list = obj.testWordListFile;
            
            net = obj.wordNetFile;
            dict = obj.dictFile;
            
            LIST_FILE       = fullfile(testFeatureFolder, 'tmp.list');
            config          = obj.configFile;
            
            disp ('Now testing with 7 mixture HMMs')
            mod_file = fullfile(obj.hmmFolder, 'hmm36', 'models');
            mac_file = fullfile(obj.hmmFolder, 'hmm36', 'macros');
            
            
            RESULTS_DIR = testFeatureFolder;
            TEST_FEAT_ROOT  = testFeatureFolder;
            
            N1list =  fullfile(testFeatureFolder, 'list.scp');
            
            ftest = fopen(N1list,'r');
            S = textscan(ftest,'%s','Delimiter','\n');
            fclose(ftest);
            
            flist = fopen(LIST_FILE,'w');
            for I = 1:size(S{1},1)
                str  = fullfile(TEST_FEAT_ROOT,S{1}{I});
                fprintf(flist,'%s\n',str);
            end;
            fclose(flist);
            
            cmd = ['HVite -D -H ' mac_file ' -H ' mod_file ' -S ' LIST_FILE ' -C ' config ' -w ' net ' -l ''*'' -i ' RESULTS_DIR filesep 'result.mlf ' flags ' ' dict ' ' test_word_list];
            system(cmd);
        end % ------ OF TEST
        
        %% **********************************************************
        % Get methods determining feature vector related gubbins
        %************************************************************
        function value = get.protoFile(obj)
            value = fullfile(obj.hmmFolder, 'proto_AutoGen');
        end
        
    end % ------ OF METHODS
    
    %%  *********************************************************
    %      _        _   _                       _   _               _
    %     | |      | | (_)                     | | | |             | |
    %  ___| |_ __ _| |_ _  ___   _ __ ___   ___| |_| |__   ___   __| |___
    % / __| __/ _` | __| |/ __| | '_ ` _ \ / _ \ __| '_ \ / _ \ / _` / __|
    % \__ \ || (_| | |_| | (__  | | | | | |  __/ |_| | | | (_) | (_| \__ \
    % |___/\__\__,_|\__|_|\___| |_| |_| |_|\___|\__|_| |_|\___/ \__,_|___/
    %************************************************************
    
    methods(Static)        
        %% **********************************************************
        % macroMat (matlab port of macro C code)
        % ported by NC - Nov 2011
        %************************************************************
        function macroMat(VECSIZE, PARAMETER_TYPE, infile, outfile)
          % This function takes the vFloors file that is created after
          % invoking HCompV and generates a Macro file required for further  
          % HMM training 
            
            ofp = fopen(outfile,'w');
                        
            fprintf(ofp, '~o\n<STREAMINFO> 1 %d\n', VECSIZE);
            fprintf(ofp, '<VECSIZE> %d\n', VECSIZE);
            fprintf(ofp, '<NULLD>\n<%s>\n', PARAMETER_TYPE);
            
            fp = fopen(infile);            
            tline = fgets(fp);
            while ischar(tline)
                fprintf(ofp,tline);
                tline = fgets(fp);
            end
            fclose(fp);
            fclose(ofp);
        end %---- of MACROMAT
        
        %% **********************************************************
        % spmodel_genMat (matlab port of spmodel_gen C code)
        % ported by NC - Nov 2011
        %************************************************************
        function spmodel_genMat(infile, outfile)
          % This function copies the middle state (3) from the silence model
          % and makes a sp model out of it by copying it to state 2 of the sp model.
          
          % From: http://www.voxforge.org/home/dev/acousticmodels/linux/create/htkjulius/tutorial/monophones/step-7
            %  In the last step you created HMM models that did not include an "sp"
            %  (short pause) silence model - which refers to the types of short pauses
            %  that occur between words in normal speech.  However, you did create a
            %  "sil" silence model - sil silence models are typically of longer
            %  duration, and refer to the pauses occur at the end of a sentence.
            %
            % The HTK book says that the sp model needs to have its "emitting state
            % tied to the centre state of the silence model".  What this means is that
            % you need to create a new sp model in your hmmdefs, that it will use the
            % centre state of sil, and then they both need to be 'tied' together.  For
            % a bit of background on HMMs and states, see this example.
            %
            % This can be done by copying the centre state from the sil model in your
            % hmmdefs file and adding it to the sp model, and then running a special
            % tool called HHED to 'tie' the sp model to the sil model so that they
            % share the same centre state.  The HTK book provides some background on
            % what this means, but you need an understanding of the basics of Hidden
            % Markov Modelling before tackling the HTK Book explanations 
            
            ofp = fopen(outfile,'a+'); % we append this time
            
            fprintf(ofp,'~h "sp"\n');
            fprintf(ofp,'<BEGINHMM>\n<NUMSTATES> 3\n<STATE> 2\n');
            
            %-- This block gets the hmmdef file to the MODEL
            lNow = [];
            fp = fopen(infile);
            while ~(strcmpi(lNow, '~h "sil"'))
                lNow = fgetl(fp);
            end
            %------------------------------
            
            %-- This block gets the hmmdef file to the STATE
            lNow = [];
            fp = fopen(infile);
            while ~(strcmpi(lNow, ['<STATE> ' num2str(3)]))
                lNow = fgetl(fp);
            end
            %------------------------------
            
            %%% This block puts a copy of the hmmdef file in from the
            %%% correct line as found above
            tline = fgetl(fp);
            while ~(strcmpi(tline, ['<STATE> ' num2str(4)]))
                fprintf(ofp,'%s\n', tline);
                tline = fgetl(fp);
            end
            fclose(fp);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf(ofp,'<TRANSP> 3\n0.000000e+00 1.000000e+00 0.000000e+00\n');
            fprintf(ofp,'0.000000e+00 5.000000e-01 5.000000e-01\n');
            fprintf(ofp,'0.000000e+00 0.000000e+00 0.000000e+00\n<ENDHMM>\n');                        

            fclose(ofp);
        end %---- of spmodel_genMat
        
        %% **********************************************************
        % models_1mixsilMat (matlab port of models_1mixsil C code)
        % ported by NC - Nov 2011
        %************************************************************
        function models_1mixsilMat(infile, outfile)
            % This function takes the hmmdef file and      
            % generates a HMM Model file       
            
            ofp = fopen(outfile,'w');
            for ii = 1:11
                
                %%% This block gets the hmmdef file to the correct line
                lNow = [];                                
                fp = fopen(infile);
                while ~(strcmpi(lNow, '~h "hmmdef"'))
                    lNow = fgetl(fp);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%% This block puts a numbered header
                if ii==1;  fprintf(ofp,'~h "one"\n'); end
                if ii==2;  fprintf(ofp,'~h "two"\n'); end
                if ii==3;  fprintf(ofp,'~h "three"\n'); end
                if ii==4;  fprintf(ofp,'~h "four"\n'); end
                if ii==5;  fprintf(ofp,'~h "five"\n'); end
                if ii==6;  fprintf(ofp,'~h "six"\n'); end
                if ii==7;  fprintf(ofp,'~h "seven"\n'); end
                if ii==8;  fprintf(ofp,'~h "eight"\n'); end
                if ii==9;  fprintf(ofp,'~h "nine"\n'); end
                if ii==10; fprintf(ofp,'~h "oh"\n'); end
                if ii==11; fprintf(ofp,'~h "zero"\n'); end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%% This block puts a copy of the hmmdef file in from the
                %%% correct line as found above
                tline = fgets(fp);
                while ischar(tline)
                    fprintf(ofp,tline);
                    tline = fgets(fp);
                end                
                fclose(fp); %close it as we reached EOF
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
            %%% THIS IS THE SILENCE MODEL @ THE END
            fprintf(ofp,'~h "sil"\n');
            fprintf(ofp,'<BEGINHMM>\n<NUMSTATES> 5\n');
            for kk = 2:4
                fprintf(ofp,'<STATE> %d\n<NUMMIXES> 1\n', kk);
                %-- This block gets the hmmdef file to the correct line
                lNow = [];                                
                fp = fopen(infile);
                while ~(strcmpi(lNow, ['<STATE> ' num2str(kk)]))
                    lNow = fgetl(fp);
                end
                %------------------------------
                
                %%% This block puts a copy of the hmmdef file in from the
                %%% correct line as found above
                tline = fgetl(fp);
                while ~(strcmpi(tline, ['<STATE> ' num2str(kk+1)]))
                    fprintf(ofp,'%s\n', tline);
                    tline = fgetl(fp);
                end                 
                fclose(fp);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
            fprintf(ofp,'<TRANSP> 5\n0.000000e+00 1.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00\n');
            fprintf(ofp,'0.000000e+00 6.000000e-01 4.000000e-01 0.000000e+00 0.000000e+00\n');
            fprintf(ofp,'0.000000e+00 0.000000e+00 6.000000e-01 4.000000e-01 0.000000e+00\n');
            fprintf(ofp,'0.000000e+00 0.000000e+00 0.000000e+00 7.000000e-01 3.000000e-01\n');
            fprintf(ofp,'0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00\n<ENDHMM>\n');
            
            fclose(ofp);
        end %---- of models_1mixsilMat
        
        %% **********************************************************
        % makeProtoHmm (Make a prototype HMM)
        %************************************************************
        function makeProtoHmm(filename,featureType,numFeatures,numStates)
            % filename is obvious
            % featureType is usually USER_D_A
            % numFeatures is number of features (including differences if used)
            % numStates is usually 18
            
            ofp = fopen(filename,'w');
            
            fprintf(ofp,'<BeginHMM>\n');
            fprintf(ofp,' <NumStates> %d <VecSize> %d <%s> <nullD> <diagC>\n',numStates,numFeatures,featureType');
            fprintf(ofp,' <StreamInfo> 1 %d\n',numFeatures);
            for state=2:numStates-1,
                fprintf(ofp,' <State> %d <NumMixes> 1\n',state);
                fprintf(ofp,'  <Stream> 1\n');
                fprintf(ofp,'  <Mixture> 1 1.0\n');
                fprintf(ofp,'    <Mean> %d\n',numFeatures);
                fprintf(ofp,'      ');
                fprintf(ofp,'%1.1f ',zeros(1,numFeatures));
                fprintf(ofp,'\n');
                fprintf(ofp,'   <Variance> %d\n',numFeatures);
                fprintf(ofp,'      ');
                fprintf(ofp,'%1.1f ',ones(1,numFeatures));
                fprintf(ofp,'\n');
            end
            fprintf(ofp,'<TransP> %d\n',numStates);
            transp = zeros(numFeatures);
            transp(1,2)=1;
            for state=2:numStates-2,
                transp(state,state)=0.6;
                transp(state,state+1)=0.4;
            end
            transp(numStates-1,numStates-1)=0.9;
            transp(numStates-1,numStates)=0.1;
            for state=1:numStates,
                fprintf(ofp,'%1.3e  ',transp(state,1:numStates));
                fprintf(ofp,'\n');
            end
            fprintf(ofp,'<EndHMM>\n');                      
            fclose(ofp);
        end %---- of MAKEPROTOHMM
        
        %% **********************************************************
        % createMLF - master label file - belongs with hmm class
        %************************************************************
        function createMLF(mapFileFolder)
            mlfFileName = 'labels';
            d = dir(fullfile(mapFileFolder, '*.map'));
            fid = fopen([fullfile(mapFileFolder,mlfFileName) '.mlf'],'w');
            fprintf(fid,'#!MLF!#\n');
            
            for I = 1:size(d,1)
                fprintf(fid,['"*/' d(I).name(1:end-3) 'lab"\n']);
                fprintf(fid,'sil\n');
                labels = d(I).name(5:end-5);
                for J = 1:length(labels);
                    switch labels(J)
                        case 'O'
                            S = 'oh';
                        case '1'
                            S = 'one';
                        case '2'
                            S = 'two';
                        case '3'
                            S = 'three';
                        case '4'
                            S = 'four';
                        case '5'
                            S = 'five';
                        case '6'
                            S = 'six';
                        case '7'
                            S = 'seven';
                        case '8'
                            S = 'eight';
                        case '9'
                            S = 'nine';
                        case 'Z'
                            S = 'zero';
                    end
                    fprintf(fid,S);
                    fprintf(fid,'\n');
                end;
                fprintf(fid,'sil\n.\n');
            end;
            fclose(fid);
        end % ------ OF CREATEMLF
        
        %% **********************************************************
        % createSCP - Just a file list
        %************************************************************
        function createSCP(mapFileFolder)
            scpFileName = 'list';
            d = dir(fullfile(mapFileFolder, '*.map'));
            fid = fopen([fullfile(mapFileFolder,scpFileName) '.scp'],'w');
            
            for I = 1:size(d,1)
                fprintf(fid,[d(I).name(1:end-3) 'map\n']);
            end;
            
            fclose(fid);
        end% ------ OF createSCP
        
        %% **********************************************************
        % score - Simple non-dynamic scoring
        %************************************************************
        function score(testMLFpath)            
            % First get the test result files and store the file locations
            testIdx = 1;
            temp = dir(testMLFpath);
            for nn = 1:numel(temp)
                if length(temp(nn).name) > 3
                    if strcmp(temp(nn).name(end-3:end), '.mlf')
                        testMLFfiles{testIdx} = temp(nn).name;
                        testIdx = testIdx+1;
                    end
                end
            end
            
            for fileIdx = 1:numel(testMLFfiles)
                % Extract file name and regonised words as strings -> store in recResults.
                % This code block is Matt Robertson's
                fmlf = fopen(fullfile(testMLFpath , testMLFfiles{fileIdx}),'r');
                
                MLF = textscan(fmlf,'%s','delimiter','\n');
                replaceDownPattern = '[0-9\-\.\s]';
                
                for I  = 1:size(MLF{1},1);
                    if strfind(MLF{1}{I},'.rec') >0
                        M{I} = MLF{1}{I};
                    else
                        if strfind(MLF{1}{I},'MLF') >0;
                            M{I} = MLF{1}{I};
                        else
                            M{I} = regexprep(MLF{1}{I},replaceDownPattern,'');
                        end
                    end;
                end;
                
                fclose(fmlf);
                
                A = 1;
                for I = 1:size(M,2)
                    if not(strcmpi(M{I},'sp')) && not(strcmpi(M{I},'sil')) && not(strcmpi(M{I},''));
                        M2{A} = M{I};
                        A = A + 1;
                    end;
                end;
                
                recResults = [];
                X=0;
                for I = 1:length(M2)
                    if strfind(M2{I},'MLF')
                    else
                        if strfind(M2{I},'.rec')
                            X = X+1;
                            
                            % found FileName
                            recResults(X).fileName = M2{I};
                            recResults(X).utterance = {};
                        else
                            recResults(X).utterance = [recResults(X).utterance M2{I}];
                        end;
                    end;
                end;
                
                % Make separate 2D arrays of input and output digits
                inputDigits  = zeros(numel(recResults),3);
                outputDigits = zeros(numel(recResults),3);
                for nn = 1:numel(recResults)
                    ipStr = recResults(nn).fileName(end-8:end-6);
                    
                    for kk = 1:3
                        inputDigits(nn,kk)  = cHMM.htk_str2num(ipStr(kk)); %see local function @ bottom of script
                        opStr = recResults(nn).utterance(kk);
                        outputDigits(nn,kk) = cHMM.htk_str2num(opStr{:}); %bit of a hack to do cell2str as it were
                    end
                end
                
                % Now do the scoring (simple whan data in the right format)
                % :D
                scoreArray = (inputDigits == outputDigits);
                pcSent(fileIdx) = 100*sum(all(scoreArray,2)) / numel(recResults);
                pcWord(fileIdx) = 100*sum(scoreArray(:)) / numel(scoreArray);
            end
            
            %will output table to console if used on *nix
            xlsdataFull = [{'-- File Name --', '%Sent', '%Word'}; testMLFfiles' num2cell(pcSent') num2cell(pcWord')];
            disp(xlsdataFull)
            fid = fopen(fullfile(testMLFpath,['score__' num2str(pcWord(1), '%0.1f') '__.txt']),'w');            
            fclose(fid);
        end% ------ OF SCORE
        
        %% **********************************************************
        % scoreWhole folder - make my life easier @ command line
        %************************************************************
        function scoreWholeFolder(folderToScore, searchString)
            if nargin < 2
                searchString = '*featR*';
            end            
            dirInfo = dir(fullfile(folderToScore, searchString));
            numFolders = numel(dirInfo);
            for nn = 1:numFolders;
                currentScoring = fullfile(folderToScore, dirInfo(nn).name);
                disp('');
                disp(currentScoring);
                cHMM.score(currentScoring);
            end
        end                        
        
        %**************************************************************************
        % htk_str2num - Convert strings to integers
        %**************************************************************************
        function opNum = htk_str2num(ipString)            
            switch ipString
                case {'oh' , 'O'}
                    opNum = 0;
                case {'one' , '1'}
                    opNum = 1;
                case {'two' , '2'}
                    opNum = 2;
                case {'three' , '3'}
                    opNum = 3;
                case {'four' , '4'}
                    opNum = 4;
                case {'five' , '5'}
                    opNum = 5;
                case {'six' , '6'}
                    opNum = 6;
                case {'eight' , '8'}
                    opNum = 8;
                case {'nine' , '9'}
                    opNum = 9;
                otherwise
                    assert(0) % throw error
            end
        end % ------ OF htk_str2num                      
        
    end % ------ OF STATIC METHODS
    
end % ------ OF CLASS
