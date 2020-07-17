import os
from tabulate import tabulate
import matplotlib.pyplot as plt

def plot_traj(dir, beta, dist, time):
    """Plot and save beta vs distance and beta vs time graphs."""
    img_file = os.path.join(dir, r'plots.png')
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1); ax2 = fig.add_subplot(2,1,2)
    ax1.plot(dist, beta, 'r')
    ax1.set_xscale('log')
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylabel('v/c')
    ax2.plot(time, beta, 'b')
    ax2.set_xscale('log')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('v/c')
    fig.tight_layout()
    fig.savefig(img_file)

def make_dir(folder_name):
    """Make a directory in current working directory."""
    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, folder_name)
    i = 1
    while os.path.exists(final_directory):
        final_directory = os.path.join(current_directory, folder_name + str(i))
        i += 1
    os.makedirs(final_directory)
    return final_directory

def make_trajfile(dir, beta, dist, time):
    """Make a txt file containing speed, distance and time results."""
    traj_file = os.path.join(dir, r'trajectory.txt')
    with open(traj_file, 'w') as f:
        table = tabulate({"Time (s)": time,"Beta (c)": beta,
            "Distance (m)": dist}, headers="keys", showindex = "always")
        f.write(table)

def make_varfile(dir, sail):
    """Make a txt file containing sail variables."""
    var_file = os.path.join(dir, r'variables.txt')
    with open(var_file, 'w') as f:
        variables = [[key, value] for key, value in sail.__dict__.items()]
        table = tabulate(variables)
        f.write(table)

def write_results(sail, beta, dist, time):
    """Create directory in current working directory and save motion and variable files."""
    #Create directory
    dir = make_dir(sail.name)
    make_trajfile(dir, beta, dist, time)
    make_varfile(dir, sail)
    plot_traj(dir, beta, dist, time)
