function featSelData = selectFeature(data, selGenes)

featSelData = cell(1, 2);
for m = 1 : 2
    featSelData{m} = cell(length(data{m}),1);
    for q = 1 : length(data{m})
        featSelData{m}{q} = data{m}{q}(selGenes, :);
    end
end
