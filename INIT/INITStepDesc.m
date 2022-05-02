classdef INITStepDesc
   properties
      PosRevOff
      AllowMetSecr
      HowToUsePrevResults %'ignore', 'exclude', 'essential'
      RxnsToIgnoreMask %[exch;import;simple transport;advanced transport;spontaneous;extracellular;custom;all]
                       %These only refers to reactions without GPRs (except custom)
      MetsToIgnore
      MILPParams %Cell array of MILPparams - dictates how many iterations that will be run in this step.
                 %Typically, MIPGap and TimeLimit is specified
      AbsMIPGaps %If the objective is close to zero, a percentage of that is very small. 
                 %Therefore, also set an absolut value for this (typically 10 or 20). 
                 %For practical reasons, the first number is not used
   end
   methods
      function obj = INITStepDesc(posRevOff_, AllowMetSecr_, howToUsePrevResults_, rxnsToIgnoreMask_, metsToIgnore_, MILPParams_, absMIPGaps_)
         if nargin > 0
            obj.PosRevOff = posRevOff_;
         else
            obj.PosRevOff = false;
         end
         if nargin > 1
            obj.AllowMetSecr = AllowMetSecr_;
         else
            obj.AllowMetSecr = false;
         end
         if nargin > 2
            obj.HowToUsePrevResults = howToUsePrevResults_;
         else
            obj.HowToUsePrevResults = 'essential';
         end
         if nargin > 3
            obj.RxnsToIgnoreMask = rxnsToIgnoreMask_;
         else
            obj.RxnsToIgnoreMask = [1;0;0;0;0;0;0;0];
         end
         if nargin > 4
            obj.MetsToIgnore = metsToIgnore_;
         else
            obj.MetsToIgnore = [1;0;0;0;0;0;0;0];
         end
         if nargin > 5
            obj.MILPParams = MILPParams_;
         else
            params = struct();
            params.TimeLimit = 5000;
            params.MIPGap = 0.0004;            
            obj.MILPParams = {params};
         end
         
         if nargin > 6
            obj.AbsMIPGaps = absMIPGaps_;
         else
            obj.AbsMIPGaps = 10;            
         end
         
         
         
      end
   end
end