# Contributing

This document describes how to contribute to this repository. Pull
requests containing bug fixes, updates, and extensions to the existing
tools and tool suites in this repository will be considered for
inclusion.

## How to Contribute

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Make sure you have `git` [installed](https://help.github.com/articles/set-up-git)
* Fork the repository on [GitHub](https://github.com/galaxyproject/tools-iuc/fork)
* Make the desired modifications - consider using a [feature branch](https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches).
* Make sure you have added the necessary tests for your changes and they pass.
* Make sure submitted tools meet IUC [Best Practices](https://galaxy-iuc-standards.readthedocs.io/en/latest/)
* Open a [pull request](https://help.github.com/articles/using-pull-requests)
  with these changes.

## What to contribute

* Wrappers for new [OSI licensed](https://opensource.org/licenses/alphabetical) tools
* Visualization Plugins
* Updates for tools
* Enhancements for tools (e.g. supporting new parameters)
* Bug fixes
* Documentation improvements
* New test cases

### Abandoned Tools

* For tools of general interest, the Galaxy metabolomics community is usually willing to adopt tools that
  you (the developer) are abandoning.
* If there are tools that you find useful but seem to be abandoned and not
  updated, you're welcome to create an issue recommending that the IUC adopt
  that tool.

## What not to contribute

* Things already wrapped and currently maintained by other users
* Wrappers without tests
* New datatypes
    * When possible, new datatypes should be added directly to the Galaxy
      codebase.

## Tests

Contributed tools should include test cases for all tools. They need not
necessarily cover all uses of the program, but should ensure that it is
generally working. The Galaxy Wiki has a
[page](https://wiki.galaxyproject.org/Admin/Tools/WritingTests) on writing
tests.

It is recommended to test it with [planemo](https://github.com/galaxyproject/planemo/), which provides a simple command line utility for testing functionality.

```console
$ planemo test --install_galaxy my_tool.xml
```

## Requirements for Contributions

Before a PR will be accepted, there are [some requirements](https://wiki.galaxyproject.org/Tools/BestPractices) on the
submitted code (which we will be happy to help you achieve if you need the
assistance).

* Tools must contain tests (and test-data)
* Continuous integration tests must pass: 
    * The tests must pass (`planemo test --install_galaxy my_tool.xml`)
    * The tools must pass linting by planemo (`planemo lint my_tool.xml`)
    * Any Python or R script must pass linting with `flake8` and `lintr`, respectively
* If there's a relevant paper for the tool, it should be cited in a [citation](https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-citations) block
* The tool must be licensed to allow use by anyone. The OSI maintains a [list of appropriate](https://opensource.org/licenses/alphabetical) licenses
* At least one approving review by an member of the repository

