Changelog
=========

Version 2.9.0
-------------
- use NeuroM dependency of version >= 3.0
- use pytest instead of nosetests

Version 2.8.0
-------------
- Use numpy.isclose in spatial ``point_to_section_segment``. Add new keyword arguments ``rtol``,``atol`` to ``point_to_section_segment``.  (#84)
- Add functionality to plot morphologies with synapses. Move all plot functionality under
  ``plot`` package. (#82)

Version 2.7.0
-------------
- Fix compatibility with numpy>=1.21 (#80)
- Add missing entries of CHANGELOG (#79)
- Revert "Introduce a function to get section ids of sections with cut leaves" (#78)

Version 2.6.0
-------------
- Add sanitize to folder convert (#75)
- move MorphToolException back into package init.py so as not to break api (#75)
- Introduce a function to get section ids of sections with cut leaves (#76)

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
- automatic morphology alignment (#51)
- Introduce installation extras to lighten the package as dependency by default (#57)
