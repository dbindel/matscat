function [elt] = extract_eltV1(elt)

for k = 1:length(elt)
  elt(k).V     = elt(k).V1;
  elt(k).Vtype = elt(k).V1type;
end
