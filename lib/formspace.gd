# #! @Arguments matrix group, scalars, unitary
# #! @Returns a basis of the formspace preserved by <A>group</A> modulo <A>scalars</A> consisting of bilinear forms if <A>unitary = false</A> and unitary forms if <A>unitary = true</>
DeclareOperation("PreservedFormspace", [IsMatrixGroup, IsVector and IsFFECollection, IsBool]);

DeclareOperation("PreservedFormspace", [IsMatrixGroup]);