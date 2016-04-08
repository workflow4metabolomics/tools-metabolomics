Instructions to integrate the "NMR bucketing" tool into a local instance of Galaxy
Version February 2015 M Tremblay-Franco


## --- R bin and Packages : --- ##
R version 3.0.2 (2013-09-25) -- "Frisbee Sailing"
Platform: x86_64-redhat-linux-gnu (64-bit)

Install the "batch" library, necessary for parseCommandArgs function and the "pracma" library, nessecary for cumtrapz function:
 - Download package source (*.tar.gz file) from your favorite CRAN (http://www.r-project.org/)
For example: http://cran.univ-lyon1.fr/

 - Install package in your R session
install.packages("path/package_name.tar.gz",lib="path",repos=NULL)
For Example: install.packages("/usr/lib64/R/library/pracma_1.8.3.tar",lib="/usr/lib64/R/library",repos=NULL)

 - Finally, load the package into your R session
library(batch)
library(pracma)


## --- Config : --- ##
 - Edit the file "/galaxy/dist/galaxy-dist/tool_conf.xml" and add 
<section id="id_name" name="Name">
  <tool file="path/NmrBucketing_xml.xml" />
</section>
to create a new section containing the NMR_Bucketing tool
or add
  <tool file="path/NmrBucketing_xml.xml" />
in an existing section

 - Put the three files NmrBucketing_xml.xml, NmrBucketing_wrapper.R and NmrBucketing_script.R in a same directory
For example, path=/galaxy/dist/galaxy-dist/tools/stats

 - Edit the NmrBucketing_xml.xml file and change the path in the following lines
    # R script
    R --vanilla --slave --no-site-file --file=path/NmrBucketing_wrapper.R --args
    
    ## Library name for raw files storage
    library path/$library

## --- XML help part --- ##
one image: 
Copy the 'Mth_Architecture_Repertoire_Bruker.png' file within the directory to your galaxy-dist/static/images/


 - Activate the "user_library_import_dir" in your /galaxy/dist/galaxy-dist/universe_wsgi.ini and create the users directories in this path, for example:

	#In universe_wsgi.ini
 		user_library_import_dir = /projet/sbr/galaxy/import/user

	#Create the user "myaccount" in this path

		User path: /projet/sbr/galaxy/import/user/myaccount@sb-roscoff.fr




Finally, restart Galaxy