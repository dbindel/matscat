function [elt] = extract_eltR(elt)

for k = 1:length(elt)
  elt(k).V     = elt(k).R;
  elt(k).Vtype = elt(k).Rtype;
end
