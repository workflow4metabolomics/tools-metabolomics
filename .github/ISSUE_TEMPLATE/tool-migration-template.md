---
name: Tool Migration Template
about: Provide a form to collect all the steps and information needed to migrate a
  Galaxy tool wrapper from outside this repository into this repository.
title: "[MIGRATION]: "
labels: ''
assignees: ''

---

# Summary

- Tool:
- Link to external repository containing the wrapper XML: 
- Link to the toolshed repository (if there is no shed.yml file)

# Steps
These are the default steps to be completed - if something is not applicable, please remove it or mark it as completed.

## Setup
- [ ] Copy the code from the other repository into a sub-folder under tools/
- [ ] Check `planemo lint .` and `planemo test .`

## Make `planemo lint .` and `planemo test .` pass
- [ ] Check and/or update dependencies
- [ ] Update the Galaxy profile version
- [ ] Update the tool wrapper to recent standards
- [ ] Add citations -> default should be latest W4M publication
- [ ] Update tool wrapper version

## Improvements

This is an [example](https://github.com/workflow4metabolomics/tools-metabolomics/blob/master/tools/xcms/xcms_plot_raw.xml) of a recent tool wrapper that complies with recent best practices.

- [ ] Make the tool work with Pulsar (see [including required files](https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-required-files-include)
- [ ] Add [EDAM](https://edamontology.github.io/edam-browser/#topic_0091) topics and operations (see [here](https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-edam-topics))
- [ ] Add creator tag (see [here](https://docs.galaxyproject.org/en/latest/dev/schema.html#tool-creator))
```
<organization
  url="https://workflow4metabolomics.org/"
  email="workflow4metabolomics@proton.me"
  name="Workflow4Metabolomics" />
```
- [ ] Update help text and parameter documentation

## Cleanup
- [ ] Archive the other repository and update the README.md file there
