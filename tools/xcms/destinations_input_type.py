import logging


'''
This file must be placed in lib/galaxy/jobs/rules/
'''


def input_type(job):
    '''
    This function checks the input file format/extension and decide which
    destination in the job_conf.xml using
     - If it's a zip file, we will launch the job in multi-thread mode
        (-pe thread 8)
     - If it's an individual file (mzxml, mzml, mzdata or netcdf),
        the job will use for instance (-pe thread 1)
    '''
    log = logging.getLogger(__name__)
    indata = dict([(da.name, da.dataset) for da in job.input_datasets])
    indata.update([(da.name, da.dataset) for da in job.input_library_datasets])
    input_extension = indata["input"].extension
    if 'input' in indata:
        input_extension = indata["input"].extension
        log.debug("The input extension is %s" % input_extension)
        if input_extension in ["mzxml", "mzml", "mzdata", "netcdf"]:
            return 'thread1-mem_free8'
    # zip file
    return 'thread8-mem_free16'
