function structtest = verifyEqualOct(testCase,actual,expected)
    % TODO: doc
if ~isOctave && ~verLessThan('matlab','8.1') % function introduced R2013a    
    verifyEqual(testCase,actual,expected);
else
    structtest=-1;
    structtest=structTest(expected,actual,structtest);
    if lt(structtest,1)
        error('Output is not as expected, more extensive debugging required')
    end
end
end

function structtest = structTest(expected,actual,structtest)
if structtest==0
    return
end
try
    
switch class(expected)
    case {'double', 'single', 'logical', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'}
        if all(((actual./expected)-1)<0.000001) % allow 0.0001% deviation,
            % as Octave imports with less precision?
            structtest=true;
        else
            structtest=false;
        end
    case {'char', 'string'}
        if all(strcmp(expected,actual))
            structtest=true;
        else
            structtest=false;
        end
    case 'struct'
        fn = fieldnames(expected);
        for i=1:numel(fn)
            structtest=structTest({expected.(fn{i})},{actual.(fn{i})},structtest);
        end
    case 'cell'
        if iscellstr(expected)
            if all(strcmp(expected,actual))
                structtest=true;
            else
                structtest=false;
            end
        else
            [a,b]=size(expected);
            for i=1:a
                for j=1:b
                    structtest=structTest(expected{i,j},actual{i,j},structtest);            
                end
            end
        end
    end
catch
    error('Output is not as expected, more extensive debugging required')
end
end
