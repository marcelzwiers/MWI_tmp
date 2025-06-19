function ConcatStructure = ConcactStructure(Struct1,dim)
% ConcatStructure = ContactStructure(Struct1, dim)
% function transsforms an array of structures and automatically concatenates
% the fields of the structures in a given dimension.
if ~isstruct(Struct1)
    error('Input must be a structure array.');
end
if ~isscalar(dim) || ~isnumeric(dim) || dim < 1
    error('Dimension must be a positive integer.');
end


fieldsStruct1 = fieldnames(Struct1);
ConcatStructure = struct();
for i = 1:length(fieldsStruct1)
    ConcatStructure.(fieldsStruct1{i}) = [];
end
for k= 1:length(Struct1)

    for i = 1:length(fieldsStruct1)

        ConcatStructure.(fieldsStruct1{i}) = ...
        cat(dim, ConcatStructure.(fieldsStruct1{i}) , Struct1(k).(fieldsStruct1{i}));
   
    end

end

