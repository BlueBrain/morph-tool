Changelog
=========

Version 2.5.1
-------------
- Update to use NeuroM v2 (#62)
- Fix apical_point_section_segment and point_to_section_segment with NeuroM v2 (#72)
- Fix compatibility with NeuroM>=2.1 (#65)
- Introduce 'sanitize' flag for conversion (#66)
- Introduce axon point module (#54)
- Allow to pass numpy rangom generators to graft_axon (#69)
- Reduce the dependencies by putting 'bluepyopt' dependency behind 'nrn' extras (#63)

Version 2.4.7
-------------
Use this version if you want to work with NeuroM v1 and MorphIO v2

- Fix upper constrains of MorphIO and NeuroM to keep old release stable

Version 2.4.1
-------------
- Fix NestedPool for py38 and CI (#60)

Version 2.4.0
-------------
- Make MorphDB.df hashable (#53)
- Improve feature (#55)
- Fix resampling issue from numerical imprecision (#58)
- automatic morphology alignement (#51)
- Introduce installation extras to lighten the package as dependency by default (#57)
