classdef INITStepDesc
% Describes a step in the ftINIT algorithm. A cell array of objects of this
% class is used as input to ftINIT to specify how the algorithm should be run.
   properties
      PosRevOff
      AllowMetSecr
      HowToUsePrevResults %'ignore', 'exclude', 'essential'
      RxnsToIgnoreMask %Specifies reactions to leave outside the problem in
                       %the MILP.
                       % [b1,b2,b3,b4,b5,b6,b7,b8], bx is either 0 or 1, 
                       % where 1 means that the group is excluded.
                       % b1 - Exchange rxns
                       % b2 - Import rxns without GPRs (from s into the cell)
                       % b3 - Simple transport reactions without GPRs (moves one metabolite between compartments)
                       % b4 - Advanced transport reactions without GPRs (moves metabolites between compartments, more complex function such as antiporter)
                       % b5 - Spontaneous reactions
                       % b6 - Reactions in the s compartment without GPRs
                       % b7 - Customly specified rxns (sent in when generating prepData)
                       % b8 - All rxns without GPRs

      MetsToIgnore % Structure describing mets that can be removed from the model
                   % before running ftINIT, such as water etc.
                   % .simpleMets
                   %      .mets        Names of metabolites to remove
                   %      .compsToKeep Compartments for which metabolites should be kept.
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