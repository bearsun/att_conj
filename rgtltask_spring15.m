function rgtltask_spring15(debug)

% Upadted 3/19/15 by Liwei, this version mixed high-load and low-load in
% the same block. So for each block there will be 12 trials of high-load
% trials and 12 trials of low-load trials, and before each trial there will
% be a cue for what kind of stimuli they will see.
% Updated 3/19/15 afternoon by Liwei: also make the rsvp begins ahead of
% the main task & exclude "O" from all letters
% Upadted 3/28/15 change some features back, now rsvp and conjunction task
% start together, for pretest of RSVP, use 5 back to check accuracy,
% h_threshold changed to .4, l_threshold changed to .8

global abbreviatedFilename; % for the EYETRACKER
global monitorh
global distance
global rect
monitorh=30; %12;% in cm
distance=55; %25;% in cm
mainscreen=1;
rng('shuffle');
nblocks=6;
ttrials=150;
ntrialsperb=ttrials/nblocks;
ntletterange=0:4;%the number of options should be mod(ttrials,*)==0
h_thres=2/5;
l_thres=4/5;
nback=5;
trialduration=2.5;% in s
%pretrialduration=1;% in s, for rsvp
framerate=Screen('FrameRate',mainscreen);
reward=.05;
penalty=0;
money=0;

% Keyboard setting
kspace = KbName('space'); kesc = KbName('Escape');
kleft = KbName('Left'); kright = KbName('Right');% response:left/yes,right/no
kdown = KbName('Down');kup=KbName('Up');
kreturn=KbName('Return');
kback = KbName('BackSpace');
% Num Keyboard
kn1=KbName('KP_End');
kn2=KbName('KP_Down');
kn3=KbName('KP_Next');
kn4=KbName('KP_Left');
kn5=KbName('KP_Begin');
kn6=KbName('KP_Right');
kn7=KbName('KP_Home');
kn8=KbName('KP_Up');
kn9=KbName('KP_Prior');
kn0=KbName('KP_Insert');

knseries=[kn0 kn1 kn2 kn3 kn4 kn5 kn6 kn7 kn8 kn9];
possiblekn=knseries;% all keys are possible
possiblekn2=[kleft kright kdown];
letters=['A':'N','P','R':'Z'];%no Q,O
target=numel(letters); %temp targets

% specify how many targets should show in each trial
ntletters=BalanceTrials((nblocks+2)*ntrialsperb,1,ntletterange);
ntletters=ntletters(1:((nblocks+2)*ntrialsperb));
ntletters=reshape(ntletters,[(nblocks+2),ntrialsperb]);

% construct the letter array without targets, template for all sequences
%pretrialframes=pretrialduration*framerate;
trialframes=trialduration*framerate;%in frames
frameperletter=1:12;%in frames, adjust to threshold
ndurations=numel(frameperletter);
SeqLetters=cell(nblocks+2,ntrialsperb,ndurations);

for block=1:(nblocks+2)
    disp(block);
    for trial=1:ntrialsperb
        ntargets=ntletters(block,trial);
        for kduration=1:ndurations
            duration=frameperletter(kduration);
            nletters=ceil(trialframes/duration);
            % line up non-target letters first
            c=1;
            letterspt=Randi(numel(letters)-1,[1,nletters]);%1:23 non-targets, 24 for target
            while c
                c=0;
                letterspt=Shuffle(letterspt);
                for letter=2:nletters  %checking for repeat one by one
                    if letterspt(letter)==letterspt(letter-1)%if repeat, reshuffle
                        c=1;break;
                    end
                end
            end
            % put targets in
            targetspt=zeros(1,nletters-2);%no targets at the beginning
            targetspt(1:ntargets)=1;
            % calculate the attentional blink period (150ms-450ms)*changed
            % to (0-450ms) according to Peter
            blink0=0;%:duration:150;
            blinkend=0:duration:30;
            nblink0=numel(blink0);
            nblinkend=numel(blinkend)-1;
            r=1;
            while r
                r=0;
                targetspt=Shuffle(targetspt);
                for letter=1:(nletters-3)
                    if targetspt(letter)==1
                        if targetspt(letter+1)==1 %no repeat
                            r=1;break;
                        end
                        if letter+nblink0>(nletters-2) % totally out of range
                            continue;
                        elseif letter+nblinkend>(nletters-2) % the end of the blink out of range
                            rnblinkend=nletters-2-letter;
                        else
                            rnblinkend=nblinkend;
                        end
                        
                        for blink=nblink0:rnblinkend
                            if targetspt(letter+blink)==1
                                r=1;h=1;break;
                            end
                        end
                        if exist('h','var')&&h==1
                            h=0;break;
                        end
                    end
                end
            end
            targetspt=[0 0 targetspt];
            letterspt(targetspt==1)=target;
            SeqLetters(block,trial,kduration)={letterspt};
        end
    end
end

disp('pass_rng_letters');
% Colors
gray = [127 127 127 255]; white = [255 255 255 255]; black = [0 0 0 255];
% red=[175,87,87,255];
% green=[69,138,69,255];
blue=[95,95,189,255];
yellow=[121,121,61,255];
truegreen=[0 255 0 255]; truered=[255 0 0 255];

fixsi = 8;
textcolor=black;
nrings=2;
stimPerRing=8;
nballs=nrings*stimPerRing;

targetoptions=[1:nballs,zeros(1,nballs/2+1)]; %catch trials 36%,inner ring 32%, outer ring 32%
repeat=1;
while repeat
    [targetindex,loadseq]=BalanceTrials(ttrials,1,targetoptions,['h','l']);
    targetindex=(reshape(targetindex,ntrialsperb,nblocks))' ;
    repeat=0;
    for i=2:ntrialsperb
        w=targetindex(:,i-1)==targetindex(:,i);
        if any(w)&&any(targetindex(w,i-1)~=0)
            repeat=1;
            break;
        end
    end
end

corkeyarr=(targetindex<9)*kleft;%inncer ring
corkeyarr(targetindex>8)=kright;%outer ring
corkeyarr(targetindex==0)=kdown;%catch

loadseq=(reshape([loadseq;BalanceTrials(2*ntrialsperb,1,['h','l'])],ntrialsperb,nblocks+2))';

disp('pass_rng_ball_target_load');


%for target letter for each trial
ltargets=Shuffle(1:numel(letters));
for i=2:ceil(ntrialsperb*(nblocks+2)/numel(letters))
    while 1
        k=Shuffle(1:numel(letters));
        if k(1)~=ltargets(end);
            break;
        end
    end
    ltargets=[ltargets,k]; %#ok<AGROW>
end

ltargets=ltargets(1:(ttrials+2*ntrialsperb));
ltargets=reshape(ltargets',ntrialsperb,nblocks+2);
ltargets=ltargets';

disp('pass_rng_letter_target');

% Login screen
prompt = {'Outputfile', 'Subject', 'age', 'gender','group','session'};
defaults = {'Training', '99', '18','M','1','1'};
answer = inputdlg(prompt, 'Training', 2, defaults);
[output, subnum, subage, gender,group,session] = deal(answer{:});

ksession=str2double(session);

if ksession==1&&~exist([pwd,'/subinfo/',subnum,'_rg.mat'],'file')
    rg=IsoRGBY;
    save([pwd,'/subinfo/',subnum,'_rg.mat'],'rg');
    red=[rg(1,:),255];
    green=[rg(2,:),255];
else
    rg=load([pwd,'/subinfo/',subnum,'_rg.mat']);
    
    red=[rg.rg(1,:),255];
    green=[rg.rg(2,:),255];
end


outputname = [output '-' subnum gender subage '-g' group '-s' session];
if exist(outputname,'file')==2&&(str2double(subnum)~=99)
    fileproblem = input('That file already exists! Append a .x (1), overwrite (2), or break (3/default)?');
    if isempty(fileproblem) || fileproblem==3
        return;
    elseif fileproblem==1
        outputname = [outputname '.x'];
    end
end
outfile = fopen(outputname,'w');
fprintf(outfile,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t \n','subnum', 'subage', 'gender', 'group', 'session', 'block', 'trial', 'load','answer', 'keypressed', 'cor', 'corate','duration','rgans','rgkp','rgcor','money');
abbreviatedFilename=[subnum,'_',datestr(today,'mmdd')];
[mainwin,rect]=Screen('OpenWindow',mainscreen,gray);
Screen('BlendFunction', mainwin, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

pbuffer=NaN(trialframes,1);
for ibuffer=1:trialframes
    [pbuffer(ibuffer),~]=Screen('OpenOffscreenWindow',mainwin,gray);
    Screen('BlendFunction', pbuffer(ibuffer), GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
end

%GENERATE T AND L STIMULI

[tMat,~,tALPHA]=imread('newT.png','PNG');
[lMat,~,lALPHA]=imread('newL.png','PNG');
tImg=cat(3,tMat,tALPHA);
lImg=cat(3,lMat,lALPHA);

T=Screen('MakeTexture',mainwin,tImg);
L=Screen('MakeTexture',mainwin,lImg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calibration
if ~debug
    Eyelink('Shutdown');
    Eyelink('Initialize');
    HideCursor;
    Screen('Fillrect', mainwin, gray)
    Screen('Flip',mainwin);
    
    Eyelink('StartSetup')
    pause(2)
    
    
    whichKey=0;
    
    keysWanted=[kspace kreturn kback];
    flushevents('keydown');
    while 1
        pressed = 0;
        while pressed == 0
            [pressed, ~, kbData] = kbcheck;
        end;
        
        for keysToCheck = 1:length(keysWanted)
            if kbData(keysWanted(keysToCheck)) == 1
                
                keyPressed = keysWanted(keysToCheck);
                if keyPressed == kback
                    whichKey=9;
                    flushevents('keydown');
                    waitsecs(.1)
                elseif keyPressed == kspace
                    whichKey=1;
                    flushevents('keydown');
                    waitsecs(.1)
                elseif keyPressed == kreturn
                    whichKey=5;
                    flushevents('keydown');
                    waitsecs(.1)
                else
                end
                flushevents('keydown');
                
            end;
        end;
        
        if whichKey == 1
            whichKey=0;
            [~, tx, ty] = Eyelink('TargetCheck');
            Screen('FillRect', mainwin ,black, [tx-20 ty-5 tx+20 ty+5]);
            Screen('FillRect', mainwin ,black, [tx-5 ty-20 tx+5 ty+20]);
            Screen('Flip', mainwin);
        elseif whichKey == 5
            whichKey=0;
            Eyelink('AcceptTrigger');
        elseif whichKey == 9
            break;
        end
    end;
    status = Eyelink('OpenFile',abbreviatedFilename);
    if status
        error(['openfile error, status: ',status]);
    end
    Eyelink('StartRecording');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assigning group/session
g=str2double(group)-1;
leftcolor=red;
rightcolor=green;
% if s==0
%     bsequence=repmat([1,0],1,nblocks/2);%0 for R/G, 1 for T/L
% elseif s==1
%     bsequence=repmat([0,1],1,nblocks/2);%0 for R/G, 1 for T/L
% else
%     error('wrong s...');
% end

if g==0 %R/G with low load, T/L with high load
    btl='h';%0 for High load, 1 for low load
    brg='l';
elseif g==1 %R/G with high load, T/L with low load
    btl='l';
    brg='h';
else
    error('wrong g...');
end

% basic screen
center = [rect(3)/2 rect(4)/2];
fixRect = CenterRect([0 0 fixsi fixsi], rect);

% construct stimuli
corticalStimSize=5;%in mm
proximalStimDist=.5;% in degrees!!
stimSeparation=.6; %in degrees, to be scaled
jitterDistance=.4;%in degrees; to be scaled

%empty variable
stimLocation = NaN(nblocks,ntrialsperb,nrings,stimPerRing,4);
stimSize=NaN(2,4);
eccentricity=NaN(nrings);

separationAngle=360/stimPerRing;
compass = separationAngle:separationAngle:360; %rectangle
jpix=ang2pix(jitterDistance);

%enlarge for placeholder
jph=[-jpix,-jpix,jpix,jpix]';

%ecc in degree to size
for ring=3:4
    ecc=proximalStimDist+stimSeparation*ring^2;
    gratingSize=CorticalScaleFactor(corticalStimSize,ecc);
    stimSize(ring-2,:)=[0 0 1 1] * ang2pix(gratingSize);
    eccentricity(ring-2)=ang2pix(ecc);
end

%regenerate positions each block with jitters
for block=1:nblocks
    for trial=1:ntrialsperb
        for ring=1:nrings
            for stimIndex=1:stimPerRing
                stimX=eccentricity(ring)*cosd(compass(stimIndex));
                stimY=eccentricity(ring)*sind(compass(stimIndex));
                stimLocation(block,trial,ring,stimIndex,:)=CenterRectOnPoint(stimSize(ring,:),center(1)+stimX,center(2)+stimY);
            end
        end
    end
end

%for the fillarc function
righthalf = 180;
lefthalf = -180;
startAngle = 0;
fullAngle=360;

%Welcome screen
HideCursor;
Screen('TextSize', mainwin, 20);
Screen('DrawText', mainwin, 'Welcome to the Tse Lab.', center(1)-350, center(2)-300, textcolor);
Screen('DrawText', mainwin, 'Please focus on the center of the screen and count the presence of the letter specified before every trial.', center(1)-350, center(2)-100, textcolor);
Screen('DrawText', mainwin, 'When the question mark comes up, please press the number key to report the count.', center(1)-350, center(2)-50, textcolor);
Screen('DrawText', mainwin, 'Press to start!', center(1)-350, center(2)+320, textcolor);
Screen('Flip', mainwin);

while 1 %wait to start
    [keyIsDown, ~, keyCode] = KbCheck;
    FlushEvents('keyDown');
    if keyIsDown
        nKeys = sum(keyCode);
        if nKeys == 1
            if keyCode(kesc)
                session_end;return
            else keyCode(kspace)
                break;
            end
        end
    end
end
if ~debug
    Eyelink('Message','duratoin_test');
end

nlostall=0;
nlost2all=0;
nmissall=0;
% duartion test
h_acor=NaN(nback,1);
l_acor=NaN(nback,1);
h_kduration=ceil(88*framerate/1000); %6 frames 66ms
l_kduration=ceil(132*framerate/1000);%8 frames 132ms
l_slot=1;
h_slot=1;
if ~(debug==2)
    for block=7:8
        if ~debug
            Eyelink('Message','block_start');
        end
        
        for trial=1:ntrialsperb
            nlost=0;
            %for letters
            if loadseq(block,trial)=='h'
                kduration=h_kduration;
                thres=h_thres;
                acor=h_acor;
                slot=h_slot;
            elseif loadseq(block,trial)=='l'
                kduration=l_kduration;
                thres=l_thres;
                acor=l_acor;
                slot=l_slot;
            else
                error('check loadseq');
            end
            sequence=SeqLetters{block,trial,kduration};
            ntargets=ntletters(block,trial);
            ansn=knseries(ntargets+1);
            ltarget=ltargets(block,trial);
            if ~(ltarget==target) %if it is not Z, switch Z and the target letter
                mtar=(sequence==target);
                sequence(sequence==ltarget)=target;
                sequence(mtar)=ltarget;
            end
            
            if ~debug
                Eyelink('Message','trial_cue');
            end
            
            Screen('DrawText',mainwin,letters(ltarget),center(1)-5,center(2)-5,black);
            Screen('Flip',mainwin);
            
            WaitSecs(.5);
            while 1 %wait to start
                [keyIsDown, ~, keyCode] = KbCheck;
                FlushEvents('keyDown');
                if keyIsDown
                    nKeys = sum(keyCode);
                    if nKeys == 1
                        if keyCode(kesc)
                            session_end;return
                        elseif keyCode(kup)
                            break;
                        end
                    end
                end
            end
            disp(frameperletter(kduration));
            if ~debug
                Eyelink('Message','trial_start');
            end
            tStart=GetSecs;
            for i=1:trialframes
                seq=ceil(i/(frameperletter(kduration)));
%                 disp('duration');
%                 disp(frameperletter(kduration));
%                 disp('seq');
%                 disp(seq);
                Screen('DrawText',mainwin,letters(sequence(seq)),center(1)-5,center(2)-5,black);
                Screen('Flip',mainwin);
                tEnd=GetSecs;
                td=tEnd-tStart;
                if td>.02
                    nlost=nlost+1;
%                     disp(['lost frame No.: ' num2str(nlost)]);
                end
                tStart=GetSecs;
            end
            
%             disp(['lost frame this trial: ' num2str(nlost)]);
            nlostall=nlostall+nlost;
            if ~debug
                Eyelink('Message','test');
            end
            
            Screen('DrawText',mainwin,'?',center(1)-5,center(2)-5);
            Screen('Flip',mainwin);
            while 1
                [keyIsDown, ~, keyCode] = KbCheck;
                FlushEvents('keyDown');
                if keyIsDown
                    nKeys = sum(keyCode);
                    if nKeys == 1
                        if keyCode(kesc)
                            session_end;return
                        elseif any(keyCode(possiblekn))
                            keypressed=find(keyCode);
                            break;
                        end
                    end
                end
            end
            if keypressed==ansn
                money=max(money+reward,0);
                DrawFormattedText(mainwin,['Yes!\n+',num2str(reward),'\nBalance: ',num2str(money)],'center','center',truegreen);
                cor=1;
                
                if ~debug
                    Eyelink('Message','correct');
                end
            else
                money=max(money+penalty,0);
                DrawFormattedText(mainwin,['No!\n',num2str(penalty),'\nBalance: ',num2str(money)],'center','center',truered);
                cor=0;
                if ~debug
                    Eyelink('Message','incorrect');
                end
            end
            Screen('Flip',mainwin);
            %output
            acor(slot)=cor;
%             disp('cor');
%             disp(cor);
%             disp('acor');
%             disp(acor);
            corate=nanmean(acor);
            disp('corate');
            disp(corate);
            
            key=find(knseries==keypressed);
            keyansn=find(knseries==ansn);
            
            corkey=NaN;
            kp2=NaN;
            cor2=NaN;
            
            WaitSecs(.5);
            fprintf(outfile, '%s\t %s\t %s\t %s\t %s\t %d\t %d\t %s\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t \n', ...
                subnum, subage, gender, group, session, block, trial, loadseq(block,trial), keyansn, key, cor,corate,frameperletter(kduration),corkey,kp2,cor2,money);
            WaitSecs(.5);
            
            if ~rem(slot,nback)%&&~(trial>(ntrialsperb-5)&&block==8)
                if corate>thres&&kduration>1
                    kduration=kduration-1;
                elseif corate<thres&&kduration<ndurations
                    kduration=kduration+1;
                end
            end
            
            slot=slot+1;
            if slot>nback
                slot=slot-nback;
            end
            if loadseq(block,trial)=='h'
                h_kduration=kduration;
                h_thres=thres;
                h_acor=acor;
                h_slot=slot;
            elseif loadseq(block,trial)=='l'
                l_kduration=kduration;
                l_thres=thres;
                l_acor=acor;
                l_slot=slot;
            else
                error('check loadseq');
            end
            
            
            Screen('FillRect',mainwin,black,fixRect);
            Screen('Flip',mainwin);
            
            while 1 %wait to start
                [keyIsDown, ~, keyCode] = KbCheck;
                FlushEvents('keyDown');
                if keyIsDown
                    nKeys = sum(keyCode);
                    if nKeys == 1
                        if keyCode(kesc)
                            session_end;return
                        elseif keyCode(kup)
                            break;
                        end
                    end
                end
            end
            if ~debug
                Eyelink('Message','trial_end');
            end
            WaitSecs(.5);
        end
        
        if ~debug
            Eyelink('Message','block_end');
        end
        
        if block==8
            break;
        end
        
        Screen('DrawText', mainwin, 'End of the Test 1 block. Press spacebar to start the next.', center(1)-350, center(2)+320, textcolor);
        Screen('Flip', mainwin);
        
        
        
        while 1 %wait to start
            [keyIsDown, ~, keyCode] = KbCheck;
            FlushEvents('keyDown');
            if keyIsDown
                nKeys = sum(keyCode);
                if nKeys == 1
                    if keyCode(kesc)
                        session_end;return
                    else keyCode(kup)
                        break;
                    end
                end
            end
        end
        WaitSecs(.5);
    end
    
    Screen('DrawText', mainwin, 'End of the Test 2 block. Press to start the dual-task session.', center(1)-350, center(2)+320, textcolor);
    Screen('Flip', mainwin);
    while 1 %wait to start
        [keyIsDown, ~, keyCode] = KbCheck;
        FlushEvents('keyDown');
        if keyIsDown
            nKeys = sum(keyCode);
            if nKeys == 1
                if keyCode(kesc)
                    session_end;return
                else keyCode(kup)
                    break;
                end
            end
        end
    end
    WaitSecs(.5);
    
end
%%%%%%%%%%%%%%%session start
if ~debug
    Eyelink('Message','session_start');
end

WaitSecs(1);

for block=1:nblocks
    
    if ~debug
        Eyelink('Message','block_start');
    end
    
    for trial=1:ntrialsperb
        
        if loadseq(block,trial)=='h'
            kduration=h_kduration;
            if strcmp(btl,'h')
                sti='tl';
                cuecolor=black;
            elseif strcmp(brg,'h')
                sti='rg';
                cuecolor=truered;
            end
        elseif loadseq(block,trial)=='l'
            kduration=l_kduration;
            if strcmp(btl,'l')
                sti='tl';
                cuecolor=black;
            elseif strcmp(brg,'l')
                sti='rg';
                cuecolor=truered;
            end
        else
            error('check loadseq');
        end
        nlost=0;
        nlost2=0;
        nmiss=0;
        %for balls
        ti=targetindex(block,trial);
        tring=ceil(ti/stimPerRing);
        tpos=ti-(tring-1)*stimPerRing;
        %for letters
        sequence=SeqLetters{block,trial,kduration};
        ntargets=ntletters(block,trial);
        ansn=knseries(ntargets+1);
        ltarget=ltargets(block,trial);
        if ~(ltarget==target) %if it is not Z, switch Z and the target letter
            mtar=(sequence==target);
            sequence(sequence==ltarget)=target;
            sequence(mtar)=ltarget;
        end
        
        if ~debug
            Eyelink('Message','trial_cue');
        end
        disp(frameperletter(kduration));
        jit=zeros(nrings,stimPerRing,4);
        randangle=rand(nrings,stimPerRing)*360;
        % draw everything on buffer
        for ibuffer=1:trialframes
            Screen('FillRect',pbuffer(ibuffer),gray,rect);
            seq=ceil(ibuffer/(frameperletter(kduration)));
%             disp(frameperletter(kduration));
%             disp(seq);
            for ring=1:nrings
                for stimIndex=1:stimPerRing
                    if ibuffer==1
                        xjit=jpix*(rand*2-1);
                        yjit=(-1)^randi(2)*sqrt(jpix^2-xjit^2);
                        jit(ring,stimIndex,:)=[xjit;yjit;xjit;yjit];
                    end
                    if strcmp(sti,'rg')
                        if stimIndex==tpos&&ring==tring
                            Screen('FillArc',pbuffer(ibuffer),leftcolor,squeeze(stimLocation(block,trial,ring,stimIndex,:))+squeeze(jit(ring,stimIndex,:)),startAngle,lefthalf);
                            Screen('FillArc',pbuffer(ibuffer),rightcolor,squeeze(stimLocation(block,trial,ring,stimIndex,:))+squeeze(jit(ring,stimIndex,:)),startAngle,righthalf);
                        else
                            Screen('FillArc',pbuffer(ibuffer),rightcolor,squeeze(stimLocation(block,trial,ring,stimIndex,:))+squeeze(jit(ring,stimIndex,:)),startAngle,lefthalf);
                            Screen('FillArc',pbuffer(ibuffer),leftcolor,squeeze(stimLocation(block,trial,ring,stimIndex,:))+squeeze(jit(ring,stimIndex,:)),startAngle,righthalf);
                        end
                    elseif strcmp(sti,'tl')
                        if stimIndex==tpos&&ring==tring
                            Screen('DrawTexture',pbuffer(ibuffer),T,[],squeeze(stimLocation(block,trial,ring,stimIndex,:))+squeeze(jit(ring,stimIndex,:)),randangle(ring,stimIndex));
                        else
                            Screen('DrawTexture',pbuffer(ibuffer),L,[],squeeze(stimLocation(block,trial,ring,stimIndex,:))+squeeze(jit(ring,stimIndex,:)),randangle(ring,stimIndex));
                        end
                    else
                        error('check_sti');
                    end
                    Screen('DrawArc',pbuffer(ibuffer),black,squeeze(stimLocation(block,trial,ring,stimIndex,:))+jph,startAngle,fullAngle);
                end
            end
            
            Screen('DrawText',pbuffer(ibuffer),letters(sequence(seq)),center(1)-5,center(2)-5,black);
%             if ibuffer==1
%                 image1=Screen('GetImage',pbuffer(ibuffer));
%                 imwrite(image1,'demo.jpg');
%             end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        Screen('DrawText',mainwin,letters(ltarget),center(1)-5,center(2)-5,cuecolor);
        Screen('Flip',mainwin);
        
        WaitSecs(.5);
        while 1 %wait to start
            [keyIsDown, ~, keyCode] = KbCheck;
            FlushEvents('keyDown');
            if keyIsDown
                nKeys = sum(keyCode);
                if nKeys == 1
                    if keyCode(kesc)
                        session_end;return
                    elseif keyCode(kup)
                        break;
                    end
                end
            end
        end
        
        if ~debug
            Eyelink('Message','trial_start');
        end
        tStart=GetSecs;
        tlast=tStart;
        for i=1:trialframes
            
            Screen('DrawTexture',mainwin,pbuffer(i));
            [~,tonset,~,miss,~]=Screen('Flip',mainwin);
            tEnd=GetSecs;
            td=tEnd-tStart;
            if td>.02
                nlost=nlost+1;
%                 disp(['lost frame No.: ' num2str(nlost)]);
            end
            td2=tonset-tlast;
            if td2>.02
                nlost2=nlost2+1;
%                 disp(['lost2 frame No.: ' num2str(nlost2)]);
            end
            tlast=tonset;
            tStart=GetSecs;
            if miss>0
                nmiss=nmiss+miss;
%                 disp(['nmiss: ',num2str(nmiss)]);
            end
        end
        
%         disp(['lost frame this trial: ' num2str(nlost)]);
        nlostall=nlostall+nlost;
        nlost2all=nlost2all+nlost2;
        nmissall=nmissall+nmiss;
        if ~debug
            Eyelink('Message','test');
        end
        
        Screen('DrawText',mainwin,'?',center(1)-5,center(2)-5,black);
        Screen('Flip',mainwin);
        while 1
            [keyIsDown, ~, keyCode] = KbCheck;
            FlushEvents('keyDown');
            if keyIsDown
                nKeys = sum(keyCode);
                if nKeys == 1
                    if keyCode(kesc)
                        session_end;return
                    elseif any(keyCode(possiblekn))
                        keypressed=find(keyCode);
                        break;
                    end
                end
            end
        end
        if keypressed==ansn
            money=max(money+reward,0);
            DrawFormattedText(mainwin,['Yes!\n+',num2str(reward),'\nBalance: ',num2str(money)],'center','center',truegreen);
            cor=1;
            if ~debug
                Eyelink('Message','correct');
            end
        else
            money=max(money+penalty,0);
            DrawFormattedText(mainwin,['No!\n',num2str(penalty),'\nBalance: ',num2str(money)],'center','center',truered);
            cor=0;
            if ~debug
                Eyelink('Message','incorrect');
            end
        end
        Screen('Flip',mainwin);
        %output
        slot=rem(trial,nback)+1;
        acor(slot)=cor;
%         disp('cor');
%         disp(cor);
%         disp('acor');
%         disp(acor);
        corate=nanmean(acor);
        disp('corate');
        disp(corate);
        
        key=find(knseries==keypressed);
        keyansn=find(knseries==ansn);
        
        WaitSecs(1);
        
        %for visual search
        Screen('FillRect',mainwin,black,fixRect);
        Screen('Flip',mainwin);
        while 1
            [keyIsDown, ~, keyCode] = KbCheck;
            FlushEvents('keyDown');
            if keyIsDown
                nKeys = sum(keyCode);
                if nKeys == 1
                    if keyCode(kesc)
                        session_end;return
                    elseif any(keyCode(possiblekn2))
                        kp2=find(keyCode);
                        break;
                    end
                end
            end
        end
        
        corkey=corkeyarr(block,trial);
        if kp2==corkey
            %Screen('FillRect',mainwin,truegreen,fixRect);
            if cor==1
                money=max(money+reward,0);
                DrawFormattedText(mainwin,['Yes!\n+',num2str(reward),'\nBalance: ',num2str(money)],'center','center',truegreen);
            else
                DrawFormattedText(mainwin,['Yes!\n+0\nBalance: ',num2str(money)],'center','center',truegreen);
            end
            Screen('Flip',mainwin);
            cor2=1;
            if ~debug
                Eyelink('Message','rgcorrect');
            end
        else
            %Screen('FillRect',mainwin,truered,fixRect);
            DrawFormattedText(mainwin,['No!\n0\nBalance: ',num2str(money)],'center','center',truered);
            Screen('Flip',mainwin);
            cor2=0;
            if ~debug
                Eyelink('Message','rgincorrect');
            end
        end
        
        fprintf(outfile, '%s\t %s\t %s\t %s\t %s\t %d\t %d\t %s\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t \n', ...
            subnum, subage, gender, group, session, block, trial,loadseq(block,trial), keyansn, key, cor,corate,frameperletter(kduration),corkey,kp2,cor2,money);
        
        WaitSecs(1);
        
        Screen('FillRect',mainwin,black,fixRect);
        Screen('Flip',mainwin);
        
        while 1 %wait to start
            [keyIsDown, ~, keyCode] = KbCheck;
            FlushEvents('keyDown');
            if keyIsDown
                nKeys = sum(keyCode);
                if nKeys == 1
                    if keyCode(kesc)
                        session_end;return
                    elseif keyCode(kup)
                        break;
                    end
                end
            end
        end
        if ~debug
            Eyelink('Message','trial_end');
        end
        WaitSecs(.5);
    end
    if ~debug
        Eyelink('Message','block_end');
    end
    if block==6
        break;
    end
    
    Screen('DrawText', mainwin, ['End of the ' num2str(block) ' block. Press spacebar to start the next.'], center(1)-350, center(2)+320, textcolor);
    Screen('Flip', mainwin);
    while 1 %wait to start
        [keyIsDown, ~, keyCode] = KbCheck;
        FlushEvents('keyDown');
        if keyIsDown
            nKeys = sum(keyCode);
            if nKeys == 1
                if keyCode(kesc)
                    session_end;return
                else keyCode(kup)
                    break;
                end
            end
        end
    end
    WaitSecs(.5);
    
end


Screen('TextSize', mainwin, 20);
Screen('DrawText', mainwin, 'That''s the end of this session. Thank you.', center(1)-350, center(2)+320, textcolor);
DrawFormattedText(mainwin,['Your Balance: ',num2str(money)],'center','center',textcolor);
Screen('Flip', mainwin);

while 1 %wait to quit
    [keyIsDown, ~, keyCode] = KbCheck;
    FlushEvents('keyDown');
    if keyIsDown
        nKeys = sum(keyCode);
        if nKeys == 1
            session_end;return
        end
    end
end

    function session_end
        if ~debug
            Eyelink('Message','session_end');
            Eyelink('Stoprecording');
            Eyelink('CloseFile');
            Eyelink('ReceiveFile');
        end
        fclose(outfile);
        ShowCursor;
        Screen('CloseAll');
        disp(['lost frame No.: ' num2str(nlostall)]);
        disp(['lost2 frame No.: ' num2str(nlost2all)]);
        disp(['nmiss: ',num2str(nmissall)]);
        return
    end

    function pixels=ang2pix(ang)
        pixpercm=rect(4)/monitorh;
        pixels=tand(ang/2)*distance*2*pixpercm;
    end

    function stimSize=CorticalScaleFactor(corticalSize,eccentricity)
        M=.065*eccentricity+.054;
        stimSize=M*corticalSize;
    end
end