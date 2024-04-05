from flow import FlowProject
import math
import time
import numpy as np
from environment import MyUniversityCluster


class MyProject(FlowProject):
    pass

gmx_mpi = "gmx_mpi"
gro_file_start = "plasticizer.gro"
gro_file_min = "confout_min.gro"
gro_file_npt = "confoutNPT.gro"
index_file = "index.ndx"
topology = "topol.top"
min_mdp_file = "minim.mdp"
npt_mdp_file = "md_NPT.mdp"
nvt_mdp_file = "md.mdp"



c5_sim = MyProject.make_group(name="c5_sim")

def condition_c5(job): # defining my condition for C_C5 simulation
    if job.statepoint().get("sim_type") == "C5": 
        return True 





def grompp_min(): 
    ''' helper function for minim grompp '''
    cmd = (
         "{gmx} grompp -c plasticizer.gro -f {min_mdp} -n {index} -p {topol} -o em.tpr"
        .format(
            gmx=gmx_mpi,
            min_mdp = min_mdp_file,
            index = index_file,
            topol = topology,
        )
        
    )
    return cmd
    


def grompp_nvt_c5(): 
    ''' helper function for nvt grompp for c5 sims with fewer cores '''
    cmd = (
         "{gmx} grompp -c confout_min.gro -f {nvt_mdp} -n {index} -p {topol} -o topol.tpr"
        .format(
            gmx=gmx_mpi,
            nvt_mdp = nvt_mdp_file,
            index = index_file,
            topol = topology,
        )
        
    )
    return cmd


def mdrun_min_c5():
    ''' helper function for c5 mdrun with fewer cores'''
    cmd = (
        "mpirun {gmx} mdrun -s em.tpr -tableb table_a1.xvg"
        .format(
            gmx=gmx_mpi)
        )
    return cmd


def mdrun_c5():
    ''' helper function for c5 mdrun with fewer cores'''
    cmd = (
        "mpirun {gmx} mdrun -s topol.tpr -tableb table_a1.xvg"
        .format(
            gmx=gmx_mpi)
        )
    return cmd


def polystat():
    ''' helper function to perform the gmx polystat analysis for radius of gyration tensor values '''
    cmd = (
        """i=3
        echo "$i" | gmx polystat -f traj_comp.xtc -s topol.tpr -n index.ndx"""
        .format(
            gmx=gmx_mpi)
        )
    return cmd 


# flow procedure 


        
    
@c5_sim
@MyProject.pre(condition_c5)
@MyProject.post.isfile("em.tpr")
@MyProject.operation(with_job=True, cmd=True, directives = {"np" : 4})
def perform_grompp_min(job):
    with job:
        return grompp_min()


@c5_sim        
@MyProject.pre(condition_c5)
@MyProject.pre.after(perform_grompp_min)
@MyProject.pre.isfile("em.tpr")
@MyProject.post.isfile("confout.gro") # note that this post condition is potentially an issue
@MyProject.operation(with_job=True, cmd=True, directives = {"np": 4})
def perform_mdrun_c5(job):
    with job:
        return mdrun_min_c5()
    
@c5_sim
@MyProject.pre.after(perform_grompp_min)
@MyProject.pre(condition_c5)
@MyProject.pre.after(perform_mdrun_c5)
@MyProject.post.isfile("confout_min.gro")
@MyProject.operation(with_job=True, cmd=True, directives = {"np": 4})
def move_confout_c5(job):
    with job:
        cmd = "mv confout.gro confout_min.gro"
        return cmd


@c5_sim    
@MyProject.pre(condition_c5)
@MyProject.pre.after(move_confout_c5)
@MyProject.pre.isfile("confout_min.gro")
@MyProject.post.isfile("topol.tpr")
@MyProject.operation(with_job=True, cmd=True, directives = {"np": 4})
def perform_grompp_nvt(job):
    with job:
        return grompp_nvt_c5()

@c5_sim
@MyProject.pre(condition_c5)
@MyProject.pre.after(perform_grompp_nvt)
@MyProject.pre.isfile("topol.tpr")
@MyProject.post.isfile("confout.gro")
@MyProject.operation(with_job=True, cmd=True, directives = {"np": 4})
def perform_mdrun_nvt(job):
    with job:
        return mdrun_c5()
        
        
@c5_sim 
@MyProject.pre.isfile("traj.trr")
@MyProject.pre.isfile("\#*")
@MyProject.operation(with_job=True, cmd=True, directives = {"np": 4})
def clean(job):
    with job:
        """ removing the repeat files and also the .trr file which you don't need for any analysis """
        cmd = (
            """rm \#* 
            rm traj.trr """
            )
        return cmd



@c5_sim
@MyProject.pre(condition_c5)
@MyProject.pre.after(perform_mdrun_nvt)
@MyProject.operation(with_job=True, cmd=True, directives = {"np": 4})
def perform_polystat(job):
    with job:
        return polystat()



@c5_sim 
@MyProject.pre(condition_c5)
@MyProject.pre.after(perform_polystat)
@MyProject.operation(with_job=True, directives = {"np": 4})
def perform_analysis(job):
    with job:
        def rsquare(x,y,z): 
            ''' Calculates the square radius of gyration '''
            
            R = x**2 + y**2 + z**2
            
            return R 
        
        
        def acyli(x,y,z):
            '''Calculates the acylindricity'''
            A = y**2 - z**2
            
            return A
            
        # reading the .xvg file produced by polystat 
        
        no = '#'
        no1 = '@'
        no2 = '&'
        
        myfile = open('format1.txt', 'w')
        with open("polystat.xvg", "r") as test_file:
            for line in test_file:
                if not line.startswith(no):
                    
                    if not line.startswith(no1):
                        if not line.startswith(no2):
                            myfile.write(line)
        myfile.close()
        
        with open("format1.txt", 'r') as f:
            lines = f.readlines()
            time = [float(line.split()[0]) for line in lines]
            re = [float(line.split()[1]) for line in lines]
            s1 = [float(line.split()[3]) for line in lines]
            s2 = [float(line.split()[4]) for line in lines]
            s3 = [float(line.split()[5]) for line in lines]

        rgs = []
        cs = []
        
        for i, j, k in zip(s1, s2, s3):
            if np.isnan(i) == True or np.isnan(j) == True or np.isnan(k) == True:
                continue 
            
            rgs.append(rsquare(i, j, k))
            cs.append(acyli(i, j, k))
            
        av_rgs = np.average(rgs)
        av_cs = np.average(cs)
        
        err_rgs = np.std(av_rgs)/len(rgs)
        err_cs = np.std(av_cs)/len(cs)
        
        job.doc["rg"] = av_rgs
        job.doc["c"] = av_cs
        
        # now we'll calculate the configurational entropy of the PL 
        
        start = 1000000 # time to start at
        space = 200000 # time between tests
        stop = job.statepoint().get("sim_time")
        length = np.arange(start, stop, space)
        length2 = np.arange(start-space, stop-space, space)

        
        def entropy(k_b, prob):
            ''' function to calculate entropy '''
            x = []
            for i in prob:
                if i != 0:
                    x.append(i*np.log(i))
                else:
                    continue
            S = -sum(x)
            
            return S
            
        r = open('endtoend.txt', 'w')
        r.close()
        
        
        for i, p in zip(length, length2):

            end_slice = re[int(round(p/20.4)):int(round(i/20.4))]

            bins = np.histogram_bin_edges(end_slice, bins=600, range = (0, 5))
        
            values, x = np.histogram(end_slice, bins=bins, density=True)
        
            k_b = 1.380649E-23
        
            values /= values.sum()
            entropies.append(entropy(k_b, values))   
            
            wr = entropy(k_b, values)
            
            with open("endtoend.txt", "a") as f:
                f.write(str(i)+ "   " + str(entropy(k_b, values)) + "\n")
                
        job.doc["s"] = np.average(entropies)
        
        return None
        

       


if __name__ == "__main__":
    MyProject().main()

    
