# Galaxy Tool Documentation: Batch Cohort Correction

## Overview
This Galaxy tool is designed to correct batch and cohort effects in intensity measurements from scientific studies. Using a mixed-model approach, it adjusts intensity values while accounting for batch and injection order effects.

---

## Table of Contents
1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Inputs](#inputs)
5. [Outputs](#outputs)
6. [Usage Example](#usage-example)
7. [Commands Executed by Galaxy](#commands-executed-by-galaxy)
8. [Important Notes](#important-notes)
9. [Contributing](#contributing)
10. [License](#license)
11. [About](#about)

---

## Prerequisites
- **Galaxy Platform**: Ensure access to a functional Galaxy instance.  
- **R version 4.2.2**: The tool relies on R for computations.  
- Required R packages: `r-optparse`, `r-dplyr`, `r-lme4`.

---

## Installation
Download the tool from the Galaxy repository or install it directly on your Galaxy instance:

```bash
git clone https://github.com/your_name/your_project.git
```

---

## Inputs
The input file should be in CSV format and include the following columns:
- **Batch**: Batch identifier (optional for batch correction).
- **SampleID**: Sample identifier.
- **Injection_Order**: Injection order (mandatory for correction).
- **Ion1, Ion2, ...**: Intensity columns to be corrected.

**Sample Input File:**
```csv
SampleID,Batch,Injection_Order,Ion1,Ion2
1,1,5,500,300
2,1,15,520,310
3,2,25,490,290
4,2,35,505,295
```

---

## Outputs
The output will also be in CSV format, with corrected intensity values.

**Sample Output File:**
```csv
SampleID,Batch,Injection_Order,Ion1,Ion2
1,1,5,-0.2464,-0.2464
2,1,15,1.3362,1.3362
3,2,25,-0.5720,-0.5719
4,2,35,0.3269,0.3268
```

---

## Usage Example
1. Upload your CSV file to Galaxy.  
2. Select the **Batch Cohort Correction** tool in your workflow.  
3. Specify the input file and set a name for the output file.  
4. Run the job and retrieve the corrected output file.

---

## Commands Executed by Galaxy
The process will run the following command:

```bash
Rscript $__tool_directory__/executable_func.R --input $input --output $output
```

---

## Important Notes
- **Injection_Order**: Mandatory for accurate corrections.  
- **CSV Format**: Ensure the file is properly formatted with columns separated by commas.  
- Malformed or improperly formatted files will result in explicit errors.

---

## Contributing
1. Fork the repository.  
2. Create a branch for your updates.  
3. Submit a pull request.  
4. Report bugs or suggest improvements in the Issues section.

---

## License
_here we will write about the_ `LICENSE`

---

## About
### Authors:  
- **Elfried Salanon**  
📅 **Date:** 2025  
- **Marie Lefebvre**  
📅 **Date:** 2025  