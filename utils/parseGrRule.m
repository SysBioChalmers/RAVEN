function tree=parseGrRule(rule)
% parseGrRule  Parse a grRule into a syntax tree.
%
% Parses a gene-protein-reaction rule into a tree of AND/OR/gene nodes, so
% that callers can reason about its structure instead of searching the raw
% string. Gene IDs are matched as whole tokens, so IDs that contain "and" or
% "or" as a substring (RAND1), or that are prefixes of one another (10 and
% 100), are never confused with operators or with each other.
%
% Parameters
% ----------
% rule : char or string
%     a grRule, e.g. "(G1 and G2) or G3". Both the word forms ("and", "or",
%     any case) and the symbol forms ("&", "&&", "|", "||") are accepted.
%     An empty rule yields an empty tree.
%
% Returns
% -------
% tree : struct
%     a scalar struct with fields:
%
%     - type     : char — 'gene', 'and' or 'or'.
%     - id       : char — the gene ID; '' for 'and'/'or' nodes.
%     - children : cell — child nodes; {} for 'gene' nodes.
%
%     An empty rule returns a 0x0 struct array, for which isempty is true.
%
% Notes
% -----
% The grammar is the usual precedence, with AND binding tighter than OR:
%
%     expr   := term ( OR term )*
%     term   := factor ( AND factor )*
%     factor := gene | '(' expr ')'
%
% Malformed rules raise RAVEN:badGrRule rather than being silently
% misparsed. Callers that lint user input (findPotentialErrors) catch this
% and report it instead of failing.
%
% Examples
% --------
%     tree = parseGrRule('(G1 and G2) or G3');
%     tree.type            % 'or'
%     numel(tree.children) % 2
%
% See also: isDnfGrRule, grRuleToDNF, grRuleToString

if nargin<1 || isempty(rule)
    tree=emptyTree();
    return
end
rule=char(rule);
if all(isspace(rule))
    tree=emptyTree();
    return
end

toks=tokenizeGrRule(rule);
if isempty(toks)
    tree=emptyTree();
    return
end

[tree,pos]=parseExpr(toks,1,rule);
if pos<=numel(toks)
    error('RAVEN:badGrRule',['Unexpected "' tokenText(toks{pos}) '" in grRule: ' rule]);
end
end


function tree=emptyTree()
tree=struct('type',{},'id',{},'children',{});
end


function toks=tokenizeGrRule(rule)
% Split into '(' , ')' , 'AND' , 'OR' and gene tokens. Operators are matched
% as whole words, which is what keeps gene IDs containing "and"/"or" intact.
toks={};
i=1;
n=numel(rule);
while i<=n
    c=rule(i);
    if isspace(c)
        i=i+1;
        continue
    end
    if c=='(' || c==')'
        toks{end+1}=c; %#ok<AGROW>
        i=i+1;
        continue
    end
    % Read a word: everything up to whitespace or a bracket.
    j=i;
    while j<=n && ~isspace(rule(j)) && rule(j)~='(' && rule(j)~=')'
        j=j+1;
    end
    word=rule(i:j-1);
    switch lower(word)
        case {'and','&','&&'}
            toks{end+1}='AND'; %#ok<AGROW>
        case {'or','|','||'}
            toks{end+1}='OR'; %#ok<AGROW>
        otherwise
            toks{end+1}=struct('gene',word); %#ok<AGROW>
    end
    i=j;
end
end


function [node,pos]=parseExpr(toks,pos,rule)
% expr := term ( OR term )*
[node,pos]=parseTerm(toks,pos,rule);
if pos>numel(toks) || ~isOp(toks{pos},'OR')
    return
end
children={node};
while pos<=numel(toks) && isOp(toks{pos},'OR')
    [child,pos]=parseTerm(toks,pos+1,rule);
    children{end+1}=child; %#ok<AGROW>
end
node=struct('type','or','id','','children',{children});
end


function [node,pos]=parseTerm(toks,pos,rule)
% term := factor ( AND factor )*
[node,pos]=parseFactor(toks,pos,rule);
if pos>numel(toks) || ~isOp(toks{pos},'AND')
    return
end
children={node};
while pos<=numel(toks) && isOp(toks{pos},'AND')
    [child,pos]=parseFactor(toks,pos+1,rule);
    children{end+1}=child; %#ok<AGROW>
end
node=struct('type','and','id','','children',{children});
end


function [node,pos]=parseFactor(toks,pos,rule)
% factor := gene | '(' expr ')'
if pos>numel(toks)
    error('RAVEN:badGrRule',['grRule ends after an operator: ' rule]);
end
tok=toks{pos};
if isOp(tok,'(')
    [node,pos]=parseExpr(toks,pos+1,rule);
    if pos>numel(toks) || ~isOp(toks{pos},')')
        error('RAVEN:badGrRule',['Unbalanced parentheses in grRule: ' rule]);
    end
    pos=pos+1;
    return
end
if isstruct(tok)
    node=struct('type','gene','id',tok.gene,'children',{{}});
    pos=pos+1;
    return
end
error('RAVEN:badGrRule',['Expected a gene but found "' tokenText(tok) '" in grRule: ' rule]);
end


function tf=isOp(tok,op)
tf=ischar(tok) && strcmp(tok,op);
end


function txt=tokenText(tok)
if isstruct(tok)
    txt=tok.gene;
else
    txt=tok;
end
end
