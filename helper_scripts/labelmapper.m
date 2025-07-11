function target_img = labelmapper(label_img, source_value, target_value)

target_img = zeros(size(label_img));

for ii = 1:length(source_value)
    target_img(label_img==source_value(ii)) = target_value(ii);
end