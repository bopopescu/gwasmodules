"""
GWA project using full sequence data.

Quan, Bjarni, Dazhe, et al.
"""
import sys


def run_parallel(pid, call_method_id, run_id='gwas'):
        """
        If no mapping_method, then analysis run is set up.
        """
        job_id = '%s_%d_%d_%s' % (run_id, x_start_i, x_stop_i, temperature)
        file_prefix = env['results_dir'] + job_id

        #Cluster specific parameters        
        shstr = '#!/bin/bash\n'
        shstr += '#$ -S /bin/bash\n'
        shstr += '#$ -N %s\n' % job_id
        #shstr += '#$ -o %s_job_$JOB_ID.out\n' % file_prefix
        #shstr += '#$ -e %s_job_$JOB_ID.err\n' % file_prefix
        shstr += '#$ -o %s_job.out\n' % file_prefix
        shstr += '#$ -e %s_job.err\n' % file_prefix
        shstr += 'source /etc/modules-env.sh\n'
        shstr += 'module load scipy/GotoBLAS2/0.9.0\n'
        shstr += 'module load matplotlib/1.0.0\n'
        shstr += 'module load mysqldb/1.2.3\n'
        shstr += 'export GOTO_NUM_THREADS=1\n'


        shstr += "python %gwas_project.py %d %d %s %d %s" % \
                        (env['script_dir'], x_start_i, x_stop_i, temperature, call_method_id, run_id)

        #shstr += "> " + file_prefix + "_job.out) >& " + file_prefix + "_job.err\n"
        print '\n', shstr, '\n'
        script_file_name = run_id + ".sh"
        f = open(script_file_name, 'w')
        f.write(shstr)
        f.close()

        #Execute qsub script
        os.system("qsub " + script_file_name)



def run_gwas(pid, call_method_id):
        #Set up GWAS
        pass

def run():
        call_method_id = int(sys.argv[1])
        if len(sys.argv) < 3:
                pid = int(sys.argv[2])



if __name__ == '__main__':
        run()
