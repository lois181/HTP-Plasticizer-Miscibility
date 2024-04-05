import signac
import itertools
from plasticizer_class import Plasticizer  # Import your Plasticizer class here
import os

##################################### ATTENTION ##########################################

                       ### THIS IS A DRAFT FILE ###

###########################################################################################

       
import signac

# Initialize a signac project
project = signac.init_project('my_project')

# Read sublists from file
with open('parameters.txt', 'r') as file:
    for line in file:
        # Split each line into individual values
        sublist = line.split()
        
        # Convert values to appropriate types
        sublist = [int(sublist[0]), int(sublist[1]), int(sublist[2]), sublist[3], int(sublist[4]), sublist[5], int(sublist[6]), int(sublist[7])]
        
        # Create statepoint and job
        statepoint = {
            'side_chain_length': sublist[1],
            'backbone_length': sublist[0],
            'side_chain_frequency': sublist[2],
            'flexibility': sublist[3],
            'sim_time': sublist[4],
            'sim_type': sublist[5],
            'PL_num': sublist[6],
            'PL_conc': sublist[7]
        }
        
        job = project.open_job(statepoint=statepoint)
        job.init()


        


#project = signac.get_project("plasticizer_simulation")

# Iterate through all jobs in the project
for job in project:
    print(f"Job ID: {job._id}")
    
    # Access and print the state variables for the job
    state_variables = job.statepoint()
    #print(state_variables.get("side_chain_length"))
    print("State Variables:")
    for key, value in state_variables.items():
        print(f"{key}: {value}")
 
    
    # Access and print any additional metadata associated with the job
    metadata = job.doc
    print("Metadata:")
    for key, value in metadata.items():
        print(f"{key}: {value}")
    
    print("=" * 40)  # Separate each job's information
    
    # move into each job and use the plasticizer class to produce the required files 
    root = os.getcwd()
    path = job.path
    os.chdir(path)
    
    state_variables = job.statepoint()
    print("Full statepoint:", state_variables)

    
    plasticizer = Plasticizer(state_variables.get("side_chain_length"), 
                              state_variables.get("backbone_length"), 
                              state_variables.get("side_chain_frequency"),
                              state_variables.get("flexibility"),
                              state_variables.get("sim_time"),
                              state_variables.get("sim_type"),
                              state_variables.get("PL_num"),
                              state_variables.get("PL_conc"))

    print(state_variables.get("PL_conc"))
    print(state_variables.get("PL_num"))
    # Generate .gro and .itp files based on state variables and flexibility
    gro_coordinates = plasticizer.get_coordinates()[0]
    print("This is gro coords")
    print(gro_coordinates)
    gro_filename = plasticizer.get_coordinates()[1]
    itp_filename = plasticizer.generate_itp_file()
    ffitp_filename = plasticizer.generate_ffitp_file()
    mdp_filename = plasticizer.generate_mdp_file()
    top_filename = plasticizer.generate_top_file()
    ndx_filename = plasticizer.generate_ndx_file()
    tables_filename = plasticizer.generate_tables()
    
    
    # move back to the root directory for the next iteration 
    os.chdir(root)
