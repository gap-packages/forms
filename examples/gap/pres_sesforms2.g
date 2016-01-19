#Preserved Sesquilinear forms: example 2
a := [ [ -1, 0, 0, -1, 0, 1 ], [ 0, -1, -1, 0, 0, 1 ], 
       [ -1, 0, 0, 1, 0, 0 ],  [ 0, -1, 1, 0, 0, -1 ], 
       [ 0, 0, 0, 0, 0, -1 ], [ 0, -1, -1, 1, 1, 1 ] ] * One(GF(3));;
b := [ [ 1, -1, 1, -1, 1, -1 ], [ 1, 1, -1, 1, 1, 0 ], 
       [ -1, 0, 1, 0, 0, 0 ], [ 0, -1, 0, 0, 0, 1 ], 
       [ 1, 1, 1, 1, 1, 1 ], [ -1, 1, 1, 1, -1, 0 ] ] * One(GF(3));;
g := Group( a, b );
forms := PreservedSesquilinearForms( g );
Display( forms[1] );
m := GModuleByMats( [a,b], GF(3) );;
usemeataxe := MTX.InvariantBilinearForm(m);
quit;
