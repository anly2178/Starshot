def find_accel_dist_time(states, times, target_beta):
    """
    From an array of states (containing velocity and position information)
    and an array of corresponding times, return the acceleration distance
    and acceleration time to reach a target beta.
    """
    i = 0
    while i < len(times):
        time = times[i]
        beta = states[0,i]
        dist = states[1,i]
        if beta >= target_beta:
            return dist, time
        else:
            i += 1
    #If target beta not reached
    print("Target velocity was not reached.")
