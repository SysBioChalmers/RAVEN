function cModel = ravenToCobraModel(rModel)
	%requiredEquiv={'rxns','mets','S','c','genes','lb','ub','rev','grRules','rxnGeneMat','metFormulas','description'};
	%required={'metCharge','rules','subSystems'};
	%optionalEquiv={'rxnNames','metNames'};
	%optionalOtherName=containers.Map({'inchis',''eccodes'},{'metInchiString','rxnECNumbers'});
	
	[ST, I]=dbstack('-completenames');
	thisPath=fileparts(ST(I).file);
	load(fullfile(thisPath,'convFieldsCobra'));
	
	% equivalent required
	fieldOther = rmfield(rModel, intersect(fieldnames(rModel), f.requiredEquiv));
	cModel = rmfield(rModel,fieldnames(fieldOther));

	% different name required
	cModel=setfield(cModel,'mets',convertMets(rModel.mets,rModel.comps(rModel.metComps)));

	% new required
	cModel=setfield(cModel,'metCharge',zeros(numel(rModel.mets),1));
	cModel=setfield(cModel,'rules',convertRules(rModel.rxnGeneMat));
	if (isfield(rModel,'subSystems')) 
		cModel=setfield(cModel,'subSystems',rModel.subSystems);
	else
		cModel=setfield(cModel,'subSystems',repmat({''},size(rModel.rxns,1)));
	end

	% equivalent optional
	%fieldOther = rmfield(rModel, setdiff(fieldnames(rModel), f.optionalEquiv));
	fields = rmfield(rModel,fieldnames(setdiff(fieldnames(rModel), f.optionalEquiv)));
	cModel=structUpdate(cModel,fields);

	% different name optional
	if (isfield(rModel,'inchis')) cModel=setfield(cModel,f.optionalOtherName('inchis'),rModel.inchis); end
	if (isfield(rModel,'eccodes')) cModel=setfield(cModel,f.optionalOtherName('eccodes'),rModel.eccodes); end

end

function mets=convertMets(mets,metComps)
	mets=arrayfun(@(x) [mets{x} metComps{x}],[1:numel(mets)],'UniformOutput',false);
end

function rules=convertRules(rxnGeneMat)
	rules=repmat({''},size(rxnGeneMat,1),1);
	for i=1:size(rxnGeneMat,1)
		tmp=arrayfun(@(x) num2str(x),find(rxnGeneMat(i,:)>0),'UniformOutput',false);
		if(~isempty(tmp))
			tmp=cellfun(@(x) ['x(' x ')'],tmp,'UniformOutput',false);
			rules{i}=['(' strjoin(tmp,'| ') ')'];
		end
	end
end

function s_merged=structUpdate(s_old,s_new)
 	%// Remove overlapping fields from first struct%// Obtain all unique names of remaining fields,%// Merge both structs
 	s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
 	names = [fieldnames(s_merged); fieldnames(s_new)];
 	s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end