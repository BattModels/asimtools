from asimtools.job import Job

def do_nothing_asimmodule(calc_input, **kwargs):
    ''' Does nothing but creates input and output files '''
    job = Job(calc_input, {})
    job.start()
    job.complete()
