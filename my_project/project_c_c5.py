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

c_sim = MyProject.make_group(name="c_sim") 

def condition_c_c5(job): # defining my condition for C_C5 simulation
    if job.statepoint().get("sim_type") == "C_C5": 
        return True 



def make_box_c_c5(c5_num):
    ''' uses gmx insert molecules to insert PL molecules  '''
    cmd = (
        "{gmx} insert-molecules -f {gro_file} -ci {PL_gro} -nmol {c5_num} -try 10000 -o out.gro "
        .format(
            gmx=gmx_mpi,
            gro_file=PI_grofile,
            c5_num=c5_num,
            PL_gro = c5_grofile,
        )
    )
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
        "mpirun {gmx} mdrun -s topol.tpr -tableb table_a1.xvg"
        .format(
            gmx=gmx_mpi)
        )
    return cmd
    
   

def mdrun_min():
    ''' helper function to perform the long minimisation '''
    cmd = ("""sed -i 's/0.6204/1.1  /' ff-PI5.itp
    sed -i 's/0.6204/1.1  /' ff-PI.itp
    gmx_mpi grompp -f minim.mdp -c out.gro -n index.ndx -p topol.top -o em.tpr
    mpirun gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/1.1 /1.0 /' ff-PI5.itp
    sed -i 's/1.1 /1.0 /' ff-PI.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    mpirun gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/1.0 /0.9 /' ff-PI5.itp
    sed -i 's/1.0 /0.9 /' ff-PI.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    mpirun gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/0.9 /0.8 /' ff-PI5.itp
    sed -i 's/0.9 /0.8 /' ff-PI.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    mpirun gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/0.8 /0.7 /' ff-PI5.itp
    sed -i 's/0.8 /0.7 /' ff-PI.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    mpirun gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    sed -i 's/0.7 /0.6204/' ff-PI5.itp
    sed -i 's/0.7 /0.6204/' ff-PI.itp
    gmx_mpi grompp -f minim.mdp -c confout.gro -n index.ndx -p topol.top -o em.tpr
    mpirun gmx_mpi mdrun -s em.tpr -tableb table_a1.xvg
    mv confout.gro confout_min.gro""")
    
    return cmd 



@c_sim
@MyProject.pre(condition_c_c5)
@MyProject.post.isfile("out.gro")
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def box_c_c5_1(job):
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

        return make_box_c_c5(c5_num)
   

@c_sim
@MyProject.pre(condition_c_c5)
@MyProject.pre.after(box_c_c5_1)
@MyProject.post.isfile("confout_min.gro") 
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def do_min(job):
    with job:
        return mdrun_min()

    
@c_sim        
@MyProject.pre(condition_c_c5)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def perform_grompp_npt(job):
    with job:
        return grompp_npt()


@c_sim        
@MyProject.pre(condition_c_c5)
@MyProject.pre.after(perform_grompp_npt)
@MyProject.post.isfile("confoutNPT.gro")
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def perform_mdrun_npt(job):
    with job:
        return mdrun()
    
@c_sim
@MyProject.pre(condition_c_c5)
@MyProject.pre.after(perform_mdrun_npt)
@MyProject.post.isfile("confoutNPT.gro")
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def move_confout(job):
    with job:
        cmd = "mv confout.gro confoutNPT.gro"
        return cmd



@c_sim    
@MyProject.pre(condition_c_c5)
@MyProject.pre.after(move_confout)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def perform_grompp_nvt(job):
    with job:
        return grompp_nvt()

@c_sim
@MyProject.pre(condition_c_c5)
@MyProject.pre.after(perform_grompp_nvt)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def perform_mdrun_nvt(job):
    with job:
        return mdrun()


@c_sim
@MyProject.pre.isfile("traj.trr")
@MyProject.pre.isfile("\#*")
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def clean(job):
    """ removing the repeat files and also the .trr file which you don't need for any analysis """
    with job:
        cmd = (
            """rm \#*
            rm traj.trr """
            )
        return cmd

@c_sim
@MyProject.pre.after(clean)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def get_rdf_files(job):
    """ doing gmx rdf to calculate it at each time """
    with job:
        # first retreive the simulation time 
        time = job.statepoint().get("sim_time")
    
        cmd = (
        """b=100
        for i in `seq 0 100 {sim_time}`
        do
        
        c=$( echo "$i + $b" | bc )
        gmx_mpi rdf -f traj_comp.xtc -n index.ndx -s topol.tpr -o rdf_"$i".xvg -selrpos whole_mol_com -seltype whole_mol_com -rmpbc yes -ref 3 -sel 3 -b "$i" -e "$c" -cn rdf_cn_"$i".xvg -bin 0.012
        done """
            .format(
                sim_time = time,
            )
        )         
        return cmd
    
@c_sim 
@MyProject.pre.after(get_rdf_files)
@MyProject.operation(with_job=True, directives={"np": 40})
def perform_analysis(job):
    """ takes files make in get_rdf_files and calculates the average miscibility parameter """
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

        dirname = os.getcwd()
        r = 3
        pi = math.pi
        Vs = (4/3)*pi*r**3
        
        def get_numbers_from_filename(filename):
            return re.search(r'\d+',filename).group(0)
        
        ext = '.xvg'
        start = 'rdf_cn_'
        
        val_file = open('int_values.txt', 'w')
        
        for filename in os.listdir(dirname):
            if filename.endswith(ext) and filename.startswith(start):
                count = float(get_numbers_from_filename(filename))
                with open(filename, "r") as file:
                    skip_header = True
                    for line in file:
                        if skip_header:
                            if line.startswith(('#', '@', '&')):
                                continue
                            else:
                                skip_header = False
                        # Assuming data is tab-separated
                        try:
                            d = line.split()[0]
                            rdf = line.split()[1]
                            
                            Ns = rdf[249]
                            
                            # now extracting the box dimensions from confoutNPT.gro
                            
                            filename = 'confoutNPT.gro'
    
                            with open(filename, 'r') as file:
                                lines = file.readlines()
                                last_line = lines[-1].strip()  # Read the last line and remove leading/trailing whitespaces
    
                                dimensions = last_line.split()[-3:]  # Extract the last 3 values (x, y, z)
                                # Convert dimensions to floats
                                x = float(dimensions[0])
                                y = float(dimensions[1])
                                z = float(dimensions[2])
    
                            Vt = x * y * z
                            
                            d_out = np.divide(c5_num-Ns-1, Vt-Vs)
                            d_in = np.divide(Ns+1, Vs)
                            
                            ratio = d_in/d_out
                            
     
                            
                            val_file.write(str(count))
                            val_file.write("   ")
                            val_file.write(str(ratio))
                            val_file.write(space)
    
                        except:
                            pass
                        
        
        val_file.close()
        
        
        with open("int_values.txt", 'r') as f:
            lines = f.readlines()
            time = [float(line.split()[0]) for line in lines]
            rdf = [float(line.split()[1]) for line in lines]
        
        new_rdf = []
        for i, j in zip(time,rdf):
            if i > 900000:
                new_rdf.append(j)
        
        if np.average(new_rdf) < 5:
            mis = 0
        else:
            mis = 1
            
        job.doc["miscibility"] = mis
        
        return None 
    
    
@c_sim
@MyProject.pre.isfile("rdf*")
@MyProject.pre.isfile("\#*")
@MyProject.pre.after(perform_analysis)
@MyProject.operation(with_job=True, cmd=True, directives={"np": 40})
def clean2(job):
    """ removing the rdf files because there are too many """
    with job:
        cmd = (
            """rm rdf*
            rm \#* """
            )
        return cmd








if __name__ == "__main__":
    MyProject().main()

    
