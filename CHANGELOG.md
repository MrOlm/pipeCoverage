# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project (attempts to) adhere to [Semantic Versioning](http://semver.org/).

## [1.0] - 2017-09-21
### changed
- Combined pipe_coverage.awk version 0.3 and pcov_tools.py version 0.5 to make combined version 1.0

### added
- Version file

## [0.3] - 2016-11-01
### fixed
- Fixed issue in pipe_coverage.awk. When very few reads map, it was still determining the read length based on the RC called in the beginning. Now, I update RC at the end to become j (which will be RC if that value is reached). This effects the start of END

## [0.2] - 2016-10-15
### added
- Added this changelog

### fixed
- Fixed issue in pipe_coverage.awk. When dealing with buckets, I'm now subtracting 1 from the value before dividing. This is because the value 3000 belongs in bucket 0, not bucket 1. I'm essentially just making ita zero-based value, which is fine because we'll never get a length 0 scaffold or position. This effects lines 11 and 20. 
