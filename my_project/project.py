from flow import FlowProject
import math
import time
import os
import numpy as np
from scipy import integrate
import re


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
PI_grofile = "PI.gro"
c5_grofile = "plasticizer.gro"

c_sim = MyProject.make_group(name="c5_c5_sim") 

def condition_c_c5(job): # defining my condition for C_C5 simulation
    if job.statepoint().get("sim_type") == "C5_C5": 
        return True 
        

def make_box_c5_c5(c5_num):
    ''' uses gmx insert molecules to make c5_c5 box '''
    cmd = (
        "{gmx} insert-molecules -ci {gro_file} -box 12 12 12 -nmol {c5_num} -try 10000 -o out.gro"
        .format(
            gmx=gmx_mpi,
            gro_file=c5_grofile,
            c5_num=c5_num,
        )
    )
    return cmd


def mdrun_min_c5_c5():
    ''' helper function to perform the long minimisation '''
    cmd = ("""sed -i 's/0.6204/1.1  /' ff-PI5.itp
    sed -i 's/0.6204/1.1  /' ff-PI.itp
    gmx_mpi grompp -f minim.mdp -c out.gro -n index.ndx -p topol.top -o em.tpr
    gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/1.1 /1.0 /' ff-PI5.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/1.0 /0.9 /' ff-PI5.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/0.9 /0.8 /' ff-PI5.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/0.8 /0.7 /' ff-PI5.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/0.7 /0.6204/' ff-PI5.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    mv confout.gro confout_min.gro""")
    
    return cmd 

def grompp_npt(): 
    ''' helper function for npt grompp '''
    cmd = (
        "{gmx} grompp -c {gro_file} -f {npt_mdp} -n {index} -p {topol} -o topol.tpr"
        .format(
            gmx=gmx_mpi,
            gro_file=gro_file_min,
            npt_mdp = npt_mdp_file,
            index = index_file,
            topol = topology,
        )
    )
    return cmd

def grompp_nvt(): 
    ''' helper function for nvt grompp '''
    cmd = (
         "{gmx} grompp -c {gro_file} -f {nvt_mdp} -n {index} -p {topol} -o topol.tpr"
        .format(
            gmx=gmx_mpi,
            gro_file=gro_file_npt,
            nvt_mdp = nvt_mdp_file,
            index = index_file,
            topol = topology,
        )
        
    )
    return cmd


def mdrun():
    ''' helper function for generic mdrun '''
    cmd = (
        "{gmx} mdrun -s topol.tpr -tableb table_a1.xvg"
        .format(
            gmx=gmx_mpi)
        )
    return cmd



@c5_c5_sim
@MyProject.pre(condition_c5_c5)
@MyProject.post.isfile("out.gro")
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def box_c5_c5(job):
    with job:
        flex = job.statepoint().get("flexibility")
        if flex == "r":
            PL_conc = job.statepoint().get("PL_conc")
            back_len = job.statepoint().get("backbone_length")
            side_len = job.statepoint().get("side_chain_length")
            freq = job.statepoint().get("side_chain_frequency")
            if side_len == 0:
                PL_size = back_len
            else:
                PL_size = back_len + math.ceil(back_len/freq)*side_len

        if flex == "f":
            PL_conc = job.statepoint().get("PL_conc")
            back_len = job.statepoint().get("backbone_length") - 2
            side_len = job.statepoint().get("side_chain_length")
            freq = job.statepoint().get("side_chain_frequency")
            if side_len == 0:
                PL_size = back_len + 2
            else:
                PL_size = back_len + math.ceil(back_len/freq)*side_len  + 2

        c5_num = round(((PL_conc/100) * 72 * 300)/ PL_size)

        return make_box_c5_c5(c5_num)
   

@c5_c5_sim
@MyProject.pre(condition_c5_c5)
@MyProject.pre.after(box_c5_c5)
@MyProject.post.isfile("confout_min.gro") 
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def do_min(job):
    with job:
        return mdrun_min()

    
@c5_c5_sim       
@MyProject.pre(condition_c5_c5)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def perform_grompp_npt(job):
    with job:
        return grompp_npt()


@c5_c5_sim       
@MyProject.pre(condition_c5_c5)
@MyProject.pre.after(perform_grompp_npt)
@MyProject.post.isfile("confoutNPT.gro")
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def perform_mdrun_npt(job):
    with job:
        return mdrun()
    
@c5_c5_sim
@MyProject.pre(condition_c5_c5)
@MyProject.pre.after(perform_mdrun_npt)
@MyProject.post.isfile("confoutNPT.gro")
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def move_confout(job):
    with job:
        cmd = "mv confout.gro confoutNPT.gro"
        return cmd

@c5_c5_sim  
@MyProject.pre(condition_c5_c5)
@MyProject.pre.after(move_confout)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def perform_grompp_nvt(job):
    with job:
        return grompp_nvt()

@c5_c5_sim
@MyProject.pre(condition_c5_c5)
@MyProject.pre.after(perform_grompp_nvt)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def perform_mdrun_nvt(job):
    with job:
        return mdrun()


@c5_c5_sim
@MyProject.pre.isfile("traj.trr")
@MyProject.pre.isfile("\#*")
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def clean(job):
    with job:
        """ removing the repeat files and also the .trr file which you don't need for any analysis """
        cmd = (
            """rm \#*
            rm traj.trr """
            )
        return cmd

@c5_c5_sim
@MyProject.pre.after(clean)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def get_ent_files(job):
    """ doing gmx rdf to calculate it at each time """
    with job:
    
        flex = job.statepoint().get("flexibility")
        if flex == "r":
            PL_conc = job.statepoint().get("PL_conc")
            back_len = job.statepoint().get("backbone_length")
            side_len = job.statepoint().get("side_chain_length")
            freq = job.statepoint().get("side_chain_frequency")
            if side_len == 0:
                PL_size = back_len
            else:
                PL_size = back_len + math.ceil(back_len/freq)*side_len
    
        if flex == "f":
            PL_conc = job.statepoint().get("PL_conc")
            back_len = job.statepoint().get("backbone_length") - 2
            side_len = job.statepoint().get("side_chain_length")
            freq = job.statepoint().get("side_chain_frequency")
            if side_len == 0:
                PL_size = back_len + 2
            else:
                PL_size = back_len + math.ceil(back_len/freq)*side_len  + 2
    
        c5_num = round(((PL_conc/100) * 72 * 300)/ PL_size)

        cmd = (
        """for i in `seq 1 {num}`
    do
    
    echo "$i" | gmx_mpi polystat -f traj_comp.xtc -n index_ent.ndx -s topol.tpr -o poly"$i".xvg -b 500000
    
    
    done
     """.format(
                num = C5_num,
            )
        )         
        return cmd
    

@c5_c5_sim
@MyProject.pre.after(get_ent_files)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def perform_analysis(job):
    with job:
        ext = ('.xvg')
        start1 = 'poly'
        space = '\n'
        gap = "     "
        no = '#'
        no1 = '@'
        no2 = '&'
        start = 1000000 # time to start at
        space = 600000 # time between tests
        stop = job.statepoint().get("sim_time")
        length = np.arange(start, stop, space)
        length2 = np.arange(start-space, stop-space, space)
        
        dirname = os.getcwd()
        
        def get_numbers_from_filename(filename):
            return re.search(r'\d+',filename).group(0)   
        
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
        
        entropies = []
        for i, p in zip(length, length2):
            all_end_dists = []
            for index, files in enumerate(os.listdir(dirname)):
            
                if files.endswith(ext) and files.startswith(start1):
                    count = float(get_numbers_from_filename(files))
                
                    myfile = open('hold.txt', 'w')
                    with open(files, "r") as test_file:
        
                        for line in test_file:
                            
                            if not line.startswith(no):
                                if not line.startswith(no1):
                                    if not line.startswith(no2):
                                        myfile.write(line)
                    myfile.close()
                    
                    with open("hold.txt", 'r') as f:
                        lines = f.readlines()
                        time = [float(line.split()[0]) for line in lines]
                        t1 = [float(line.split()[1]) for line in lines]
        
                    end_slice_t1 = t1[int(round(p/20.4)):int(round(i/20.4))]
        
                    for j in end_slice_t1:
                        all_end_dists.append(j)  
        
            bins = np.histogram_bin_edges(all_end_dists, bins=600, range = (0, 5))
        
            values, x = np.histogram(all_end_dists, bins=bins, density=True)
        
            k_b = 1.380649E-23
        
            values /= values.sum()
            entropies.append(entropy(k_b, values))   
            
            wr = entropy(k_b, values)
            
            with open("endtoend.txt", "a") as f:
                f.write(str(i)+ "   " + str(entropy(k_b, values)) + "\n")  
        
        job.doc["entropy_PLPL"] = np.average(entropies)  
            


    
    
@c5_c5_sim
@MyProject.pre.isfile("poly*")
@MyProject.pre.isfile("\#*")
@MyProject.pre.after(perform_analysis)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def clean2(job):
    """ removing the poly files because there are too many """
    with job:
        cmd = (
            """rm poly*
            rm \#* """
            )
        return cmd

if __name__ == "__main__":
    MyProject().main()

    
