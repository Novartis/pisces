#!/usr/bin/env python
from signal import signal, SIGINT, SIGTERM
from subprocess import Popen, PIPE, call
from tqdm import tqdm
import logging
import sys
import os
import stat
import pickle
import pandas as pd
import time
from pkg_resources import get_distribution

__version__ = get_distribution("novartis_pisces").version

def _submit_drmaa(args, unknown_args):
    """ Submit multiple 'pisces run' jobs to the cluster using libdrmaa """

    if args.local:
        submit_local(args.metadata, args.config, args.max_memory, args.runtime, unknown_args, args.dry_run, args.debug, args.workdir)
    else:
        submit_drmaa(args.metadata, args.config, args.max_memory, args.runtime, unknown_args, args.dry_run, args.debug, args.workdir)

def create_job_scripts(sample_metadata, config, max_memory, max_runtime, metadata_file, unknown_args, working_dir, summarize=True):
    import pisces.cli
    
    logging.info('Creating job templates.')
    if working_dir:
        lock_file = os.path.join(working_dir, '.pisces')
    else:
        lock_file = os.path.join(os.getcwd(), '.pisces')
    for index, row in sample_metadata.iterrows():
        sample_id = row['SampleID']
        fq1 = None
        fq2 = None
        sra = None
        try:
            fq1 = row['Fastq1']
        except:
            sra = row['SRA']
        try:
            fq2 = row['Fastq2']
        except:
            fq2 = None
        if fq2 and not fq1:
            raise RuntimeError("Fastq1 column is mandatory when Fastq2 is populated at row {row} in {csv}.".format(
                row=str(row + 1), csv=metadata.filename))
        elif sra and (fq1 or fq2):
            raise RuntimeError("SRA column cannot be populated when Fastq1 or Fastq2 are used at row {row} in {csv}.".format(
                row=str(row + 1), csv=metadata.filename))
        output_dir = row['Directory']
        cmd = ['--config', config, 'run']
        if sra:
            cmd = cmd + ['-sra']
            cmd.extend(sra.split(';'))
        if fq1:
            cmd = cmd + ['-fq1']
            cmd.extend(fq1.split(';'))
        if fq2:
            cmd = cmd + ['-fq2']
            cmd.extend(fq2.split(';'))
        cmd = cmd + ['-o', output_dir]
        cmd = cmd + ['-n', sample_id] + unknown_args
        logging.info('command: %s', cmd)
        logging.info("Created job template for %s", sample_id)

        run_parser = pisces.cli.create_parser()
        run_args, run_unknown_args = run_parser.parse_known_args(cmd)

        submit_script_path = os.path.join(
            lock_file, "pisces_" + row['SampleID'] + '.sh')
        with open(submit_script_path, 'w') as submit_script:
            submit_script.write("#!/bin/sh\n")
            submit_script.write("cd %s\n" % os.getcwd())
            submit_script.write("source %s/env.sh\n" % lock_file)
            submit_script.write("/usr/bin/env python " + pisces.cli.__file__ + ' ' + ' '.join(cmd))
        st = os.stat(submit_script_path)
        os.chmod(submit_script_path, st.st_mode | stat.S_IEXEC)
        logging.info("Created script for %s", sample_id)
        yield (submit_script_path, run_args, run_unknown_args)
    if summarize:
        submit_script_path = os.path.join(lock_file, "pisces_summarize.sh")
        with open(submit_script_path, 'w') as submit_script:
            submit_script.write("#!/bin/sh\n")
            submit_script.write("cd %s\n" % os.getcwd())
            submit_script.write("source %s/env.sh\n" % lock_file)
            cmd = ['--config', config, 'summarize-expression', '-m', metadata_file]
            submit_script.write("/usr/bin/env python " + pisces.cli.__file__ + ' ' + ' '.join(cmd) + '\n')
            cmd = ['--config', config, 'summarize-qc', '-m', metadata_file]
            submit_script.write("/usr/bin/env python " + pisces.cli.__file__ + ' ' + ' '.join(cmd) + '\n')
        st = os.stat(submit_script_path)
        os.chmod(submit_script_path, st.st_mode | stat.S_IEXEC)
        logging.info("Created script for experiment summary")
    
def submit_local(metadata, config, max_memory, max_runtime, unknown_args, dry_run=False, working_dir=None, debug=False):
    sample_metadata = pd.read_csv(metadata)
    scripts = [s[0] for s in create_job_scripts(sample_metadata, config, max_memory, max_runtime, metadata, unknown_args, working_dir)]
    with tqdm(total=len(scripts), leave=False, unit='jobs run', unit_scale=True, position=0) as jobs_run:
        for script in scripts:
            p = Popen(['/bin/sh', script], shell=True)
            p.wait()
            jobs_run.update(1)
            
def submit_drmaa(metadata, config, max_memory, max_runtime, unknown_args, dry_run=False, debug=False, working_dir=None, blocking=True, summarize=True):
    """ args:
                metadata - open filehandle or bytes buffer
                config - dictionary of config values from pisces.cli
                unknown_args - "pisces run arguments as a String"
    """
    import drmaa
    job_tracker = {}
    
    
    if not working_dir:
        lock_file = os.path.join(os.getcwd(), '.pisces')
    else:
        lock_file = os.path.join(working_dir, '.pisces')

    def delete_jobs(signal=None, frame=None):
        delete = ''
        while delete.upper() not in ('Y', 'N'):
            delete = input("Delete jobs? (y/n)")
        if delete.upper() == 'Y':
            logging.info("Deleting jobs...")
            try:
                for jobid in job_tracker.keys():
                    s.control(jobid, drmaa.JobControlAction.TERMINATE)
                logging.info("Deleted %s jobs.", len(jobids))
            except NameError:
                logging.info("Jobs could not be deleted.")
        elif delete.upper() == 'N':
            logging.info(
                "You may check on the status of your jobs by running 'pisces submit' in this directory."
            )
        with open(os.path.join(lock_file, 'jobs'), 'wb') as jobs:
            pickle.dump(job_tracker, jobs)
        sys.exit()

    def submit_jobs(s, sample_metadata, summarize):
        for submit_script_path, run_args, run_unknown_args in create_job_scripts(sample_metadata, config, max_memory, max_runtime, metadata.name, unknown_args, working_dir):
            if dry_run:
                cmd = "qsub -o {out} -cwd -j y -l h_rt={runtime},m_mem_free={memory}G -pe smp {threads} -binding linear:{threads} {script}".format(
                    out=os.path.join(lock_file, "logs"),
                    runtime=str(max_runtime),
                    memory=str(max_memory // run_args.threads),
                    threads=str(run_args.threads),
                    script=submit_script_path)
                logging.info(cmd)
            else:
                jt = s.createJobTemplate()
                remote_cmd = os.path.join(submit_script_path)
                jt.remoteCommand = remote_cmd
                native_spec = ' -l h_rt={runtime},m_mem_free={memory}G -pe smp {threads} -binding linear:{threads}'
                native_spec = native_spec.format(
                    runtime=str(max_runtime),
                    memory=str(max_memory // run_args.threads),
                    threads=str(run_args.threads))
                jt.nativeSpecification = native_spec.format()
                jt.workingDirectory = os.getcwd()
                jt.joinFiles = True
                jt.outputPath = os.path.join(":" + lock_file, "logs", run_args.name + ".log")
               

                jobid = s.runJob(jt)
                s.deleteJobTemplate(jt)
                job_tracker[jobid] = drmaa.JobState.UNDETERMINED
        if summarize:
            submit_script_path = os.path.join(lock_file, "pisces_summarize.sh")
            jt = s.createJobTemplate()
            remote_cmd = os.path.join(submit_script_path)
            jt.remoteCommand = remote_cmd
            native_spec = ' -l h_rt={runtime},m_mem_free={memory}G -pe smp {threads} -binding linear:{threads} -hold_jid {hold_jobs}'
            native_spec = native_spec.format(
                runtime=str(max_runtime),
                memory=str(max_memory // run_args.threads),
                threads=str(run_args.threads),
                hold_jobs=','.join(job_tracker.keys()))
            jt.nativeSpecification = native_spec.format()
            jt.workingDirectory = working_dir
            jt.joinFiles = True
            jt.outputPath = os.path.join(":" + lock_file, "logs", "summarize.log")
           

            jobid = s.runJob(jt)
            s.deleteJobTemplate(jt)
            job_tracker[jobid] = drmaa.JobState.UNDETERMINED
        if dry_run:
            sys.exit()

    def track_job_progress(job_tracker, s):
        with tqdm(
                total=len(job_tracker),
                leave=False,
                unit='jobs run',
                unit_scale=True,
                position=0) as jobs_run:
            with tqdm(
                    total=len(job_tracker),
                    leave=False,
                    unit='jobs finished',
                    unit_scale=True,
                    position=1) as jobs_done:
                run_status = set((drmaa.JobState.DONE, drmaa.JobState.RUNNING,
                                  drmaa.JobState.FAILED))
                done_status = set((drmaa.JobState.DONE, drmaa.JobState.FAILED))

                while True:
                    jobs_run.update(
                        sum([
                            status in run_status
                            for status in job_tracker.values()
                        ]) - jobs_run.n)
                    jobs_done.update(
                        sum([
                            status in done_status
                            for status in job_tracker.values()
                        ]) - jobs_done.n)
                    for jobid, status in job_tracker.items():
                        if status in done_status:
                            continue
                        elif s.jobStatus(jobid) == drmaa.JobState.DONE:
                            job_tracker[jobid] = drmaa.JobState.DONE
                        elif s.jobStatus(
                                jobid) == drmaa.JobState.QUEUED_ACTIVE:
                            job_tracker[jobid] = drmaa.JobState.QUEUED_ACTIVE
                        elif s.jobStatus(jobid) == drmaa.JobState.RUNNING:
                            job_tracker[jobid] = drmaa.JobState.RUNNING
                        elif s.jobStatus(jobid) == drmaa.JobState.FAILED:
                            job_tracker[jobid] = drmaa.JobState.FAILED
                            logging.error("Job with job id %s failed",
                                          str(jobid))
                        else:
                            job_tracker[jobid] = 'undetermined'
                    if all(status in done_status
                           for status in job_tracker.values()):
                        logging.info("All jobs have finished.")
                        logging.info("Log files are at %s" % os.path.join(
                            lock_file, "logs"))
                        logging.info("To resubmit more jobs you must delete %s"
                                     % lock_file)
                        if all(status == drmaa.JobState.DONE
                               for status in job_tracker.values()):
                            sys.exit(0)
                        elif any(status == 'undetermined'
                                 for status in job_tracker.values()):
                            sys.exit(
                                "Some jobs may have failed. Please check output files."
                            )
                        else:
                            sys.exit("Some jobs have failed.")

                    time.sleep(10)

    if debug:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s %(name)-8s %(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M')

    logging.info("PISCES version %s", __version__)

    if blocking:
        signal(SIGINT, delete_jobs)
        
    with drmaa.Session() as s:
        if os.path.exists(lock_file):
            logging.info("loading existing jobs file from %s" % lock_file)
            with open(os.path.join(lock_file, 'jobs'), 'rb') as jobs:
                job_tracker = pickle.load(jobs)
            # remove non-existant jobs
            jobids = tuple(job_tracker.keys())
            for jobid in jobids:
                try:
                    s.jobStatus(jobid)
                except:
                    job_tracker.pop(jobid)

            #raise RuntimeError(".pisces lock file exists. If you are not running jobs you may safely delete the file and try again.")
        else:
            logging.info("creating .pisces lock file at %s" % lock_file)
            os.makedirs(lock_file)
            call('export -p > %s/env.sh' % lock_file, shell=True)

            sample_metadata = pd.read_csv(metadata)
            if not all(
                    col in sample_metadata.columns
                    for col in ('SampleID', 'Directory')):
                raise IOError(
                    "('SampleID', 'Directory') columns are all required in metadata.csv"
                )
            elif not any(
                    col in sample_metadata.columns
                    for col in ('SRA', 'Fastq1', 'Fastq2')):
                raise IOError(
                    "One of ('SRA', 'Fastq1, 'Fastq2') columns are required in metadata.csv"
                )
            submit_jobs(s, sample_metadata, summarize=summarize)
            with open(os.path.join(lock_file, 'jobs'), 'wb') as jobs:
                pickle.dump(job_tracker, jobs)
        if blocking:
            track_job_progress(job_tracker, s)
