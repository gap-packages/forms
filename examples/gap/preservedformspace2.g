#example of PreservedFormspace with given scalars
G := Sp(4,5);
group := Group(G.1 * Z(5), G.2 * Z(5)^3);
PreservedFormspace(group, [Z(5)^0, Z(5)^0], false);
PreservedFormspace(group, [Z(5)^2, Z(5)^2], false);
quit;
