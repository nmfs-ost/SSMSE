# SSMSE

# SSMSE 0.3.0

## Major changes

- [Change column names per r4ss updates](https://github.com/nmfs-ost/SSMSE/issues/216). Back-compatible with existing user input, although output names are as in [r4ss version 1.50.0](https://github.com/r4ss/r4ss/releases/tag/v1.50.0). SSMSE is now dependent on r4ss 1.50.0 or greater.
- [Add Mac OS arm64 binary](https://github.com/nmfs-ost/SSMSE/issues/208).


## Minor improvements and bug fixes

- [Refactor code](https://github.com/nmfs-ost/SSMSE/issues/192)
- Fix small issues with package warnings and notes. [PR #191](https://github.com/nmfs-ost/SSMSE/pull/191), [PR #195](https://github.com/nmfs-ost/SSMSE/pull/195). Thanks to @kellijohnson-NOAA for working on this!
- [Fix an issue with incorrect parameter devs default value](https://github.com/nmfs-ost/SSMSE/pull/187). Thanks @charliehinchliffe for identifying this issue and solution!
- [Enable default reporting of derived model quanities F, B_ratio, and SPR_ratio](https://github.com/nmfs-ost/SSMSE/pull/184)
- [Fix a bug with multiarea models](https://github.com/nmfs-ost/SSMSE/pull/181). Thanks @skylersagarese-NOAA for identifying this issue and solution!
- [Fix a issue with reading in mean body size](https://github.com/nmfs-ost/SSMSE/pull/179). Thanks @skylersagarese-NOAA for identifying this issue and solution!
- [Fix issue with reading in recdev2](https://github.com/nmfs-ost/SSMSE/pull/174). Thanks @CassidyPeterson-NOAA for identifying this issue and solution!
- [Fix issue withfleets of type=3 in the list of fleet names](https://github.com/nmfs-ost/SSMSE/pull/174/commits/52f4d04db41a454ee6acf4a7bd4618eaf58eb463). Thanks @CassidyPeterson-NOAA for identifying this issue and solution!
- [Fix issue with mean size at age data](https://github.com/nmfs-ost/SSMSE/pull/173)


# SSMSE 0.2.8

Version submitted to JOSS and Zenodo.