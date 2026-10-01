This file describes changes in the Forms package.

## 1.3.0 (2026-05-17)

- Add `ConformalSymplecticGroup` and some generalizations (#80, #86)
- Fix scalar normalization in forms recognition, which affected
  `PreservedSesquilinearForms` (#84)

## 1.2.14 (2026-01-26)

- Determine the field of definition of invariant forms of classical groups the
  same way as the GAP library, preparing for GAP storing it with
  `InvariantBilinearForm` etc. (#82)

## 1.2.13 (2025-05-05)

- Rename the undocumented helpers `IsOrthogonalMatrix`, `IsSymplecticMatrix` and
  `IsHermitianMatrix` with a `FORMS_` prefix, as their names were misleading
- Fix an error in `BaseChangeToCanonical` for forms in dimension 1 (#74)
- Use `PrimitiveRoot` instead of `PrimitiveElement` in `DiscriminantOfForm` and
  forms recognition, as only the former is guaranteed to generate the
  multiplicative group (#75)

## 1.2.12 (2024-08-30)

- Remove the undocumented `FormByPolynomial` (#47)
- Fix a too restrictive argument check in classical group constructors with a
  prescribed form, which rejected forms equal up to a scalar (#68)
- Update and correct part of the forms recognition code; keep the old
  recognition functions for compatibility with the recog package

## 1.2.11 (2024-04-05)

- Fix `BaseChangeOrthogonalBilinear`, broken in 1.2.10 (#62)
- Fix `BaseChangeSymplectic` for degenerate forms
- Speed up `BaseChangeHermitian` and `BaseChangeOrthogonalQuadratic` for large
  matrices (#55, #56)
- Make `Forms_RESET` compatible with FinInG (#61)

## 1.2.10 (2024-03-21)

- Speed up `BaseChangeSymplectic` substantially, e.g. 55 times for a 200x200
  matrix over GF(17) (#42)
- Speed up other base change computations (#46)

## 1.2.9 (2022-10-14)

- Fix an error in `BaseChangeOrthogonalBilinear` over fields with more than 256
  elements

## 1.2.8 (2022-07-09)

- Speed up `BaseChangeOrthogonalBilinear` by orders of magnitude for matrices
  with more than 100 rows, with similar changes to `BaseChangeHermitian` and
  `BaseChangeOrthogonalQuadratic`
- Validate the input of `BaseChangeOrthogonalBilinear`

## 1.2.7 (2022-03-02)

- Let classical group constructors with a prescribed form accept both plain
  matrices and matrix objects

## 1.2.6 (2021-07-29)

- Extend the GAP library functions `GO`, `GU`, etc. to accept a prescribed
  invariant form; this needs GAP 4.12 (#4)
- Add `Omega` methods taking a field `GF(q)` instead of `q`
- Start to support matrix objects; require GAP >= 4.9 (#10)
- Keep `BaseField` working for forms in GAP >= 4.11
- Drop the dependency on GAPDoc (#3)
- Fix errors in the manual (#18)

## 1.2.5 (2018-09-27)

## 1.2.4 (2017-08-26)

## 1.2.3.4 (2016-01-19)

## 1.2.3 (2015-10-26)
