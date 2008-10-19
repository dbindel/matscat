function [elt] = extract_eltV2(elt)

for k = 1:length(elt)
  elt(k).V     = elt(k).V2;
  elt(k).Vtype = elt(k).V2type;
end
