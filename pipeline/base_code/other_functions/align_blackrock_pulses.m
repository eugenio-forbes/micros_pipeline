function [behavioral_offset_ms,blackrock_offset] = align_blackrock_pulses(behavioral_pulses,blackrock_pulses,windSize)
%PULSEALIGN - Pick matching behavioral and eeg pulses.
%
% This method picks matching behavioral and eeg pulses from the
% beginning and end of the behavioral period for use with the
% logalign function to align behavioral and eeg data.
%
% This is josh's pimped version of this code that should work a lot better.... 12/9/09
%
% FUNCTION:
%   function [beh_ms,eeg_offset] = pulsealign2(beh_ms,pulses,pulseIsMS)
%
% INPUT ARGS:
%   beh_ms = beh_ms;   % A vector of ms times extracted from the
%                      %  log file
%   pulses = pulses;   % Vector of eeg pulses extracted from the eeg
%
% OUTPUT ARGS:
%   beh_ms- The truncated beh_ms values that match the eeg_offset
%   eeg_offset- The trucated pulses that match the beh_ms
%
% 1/15/07 MvV: added functionality to align neuralynx pulses. These
% pulses are already in ms. .


%JJ on how this algorithm works:
%  Step through the recorded sync pulses in chunks of  windsize.  Use corr to find the chunks of behavioral pulse times where the inter-pulse intervals are correlated.  When the maximum correlation is greater than corrThresh, then it indicates that the pairs match.

%note that sampling rate never comes in here. this is how alignment should work---it should be entirely sampling-rate independent....

%!%%
%these are parameters that one could potentially tweak....
%windSize=30;
corrThresh=.99;        
%%%Eugenio: Adjusting window size helps reduce the
%maximum deviation down to subms, so changed code to change window size
%until it gets it right (max 10 tries?

%%%

eegBlockStart=1:windSize:length(blackrock_pulses)-windSize;

% $$$ 
% $$$ for b=1:length(eegBlockStart)
% $$$   eeg_pulseTimes=pulses(eegBlockStart(b)+[0:windSize-1]);
% $$$ 
% $$$   for i=1:length(beh_ms)-length(eeg_pulseTimes)+1
% $$$     beh_pulseTimes=beh_ms(i+[0:windSize-1]);    
% $$$ %    [r(i),p(i)]=corr(diff(eeg_pulseTimes),diff(beh_pulseTimes));
% $$$     [r(i),p(i)]=corr(diff(eeg_pulseTimes),diff(beh_pulseTimes),'type','kendall');
% $$$   end
% $$$   [blockR(b),blockBehMatch(b)]=max(r);
% $$$ end



beh_d=diff(behavioral_pulses);
%  beh_d(beh_d>10*1000)=0; %if any interpulse differences  are greater than five seconds, throw them out!
pulse_d=diff(blackrock_pulses);
for b=1:length(eegBlockStart)
    eeg_d=pulse_d(eegBlockStart(b)+[0:windSize-1]);
    r=zeros(1,length(beh_d)-length(eeg_d));p=r;
    for i=1:(length(beh_d)-length(eeg_d))
        %     [r(i),p(i)]=corr(eeg_d,beh_d(i+[0:windSize-1]));
        %     r(i)=corr(eeg_d,beh_d(i+[0:windSize-1]));
        r(i)=fastCorr(eeg_d,beh_d(i+[0:windSize-1]));
    end
    [blockR(b),blockBehMatch(b)]=max(r);
end
%now, for each block, check if it had a good correlation.  if so, then add the set of matching pulses into the output

blackrock_offset=[];
behavioral_offset_ms=[];

for b=find(blockR>corrThresh)
  x=blackrock_pulses(eegBlockStart(b)+[0:windSize-1]);
  blackrock_offset=[blackrock_offset;x];
  y=behavioral_pulses(blockBehMatch(b)+[0:windSize-1]);
  behavioral_offset_ms=[behavioral_offset_ms; y];
end
end

function r=fastCorr(x,y)
%josh's faster version of corr
c=cov(x,y);
r=c(1,2)./(std(x)*std(y));
end

