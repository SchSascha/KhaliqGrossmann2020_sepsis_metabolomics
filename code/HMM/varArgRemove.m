function [value, remainder] = varArgRemove(key, defaultValue, keyValueArray)

for k = 1 : length(keyValueArray) - 1,
    if isstr(keyValueArray{k}),
        if strcmpi(key, keyValueArray{k}),
            value = keyValueArray{k+1};
            remainder = keyValueArray([1:k-1, k+2:length(keyValueArray)]);
            return;
        end
    end
end
value = defaultValue;
remainder = keyValueArray;
