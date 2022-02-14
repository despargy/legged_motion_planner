#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt

def store_multicol_data(data, step, traj_size):

    CoM_x = []
    CoM_y = []
    CoM_z = []

    lin_x = []
    lin_y = []
    lin_z = []

    S0_x = []
    S0_y = []
    S0_z = []

    S1_x = []
    S1_y = []
    S1_z = []

    for pp in range(traj_size):

        CoM_x.append(data.iloc[pp, 0])
        CoM_y.append(data.iloc[pp, 1])
        CoM_z.append(data.iloc[pp, 2])

        lin_x.append(data.iloc[pp, 3])
        lin_y.append(data.iloc[pp, 4])
        lin_z.append(data.iloc[pp, 5])

        S0_x.append(data.iloc[pp, 6])
        S0_y.append(data.iloc[pp, 7])
        S0_z.append(data.iloc[pp, 8])

        S1_x.append(data.iloc[pp, 9])
        S1_y.append(data.iloc[pp, 10])
        S1_z.append(data.iloc[pp, 11])

    list_com = [CoM_x, CoM_y, CoM_z]
    list_lin = [lin_x, lin_y, lin_z]
    list_foot_pos = [ [S0_x, S0_y, S0_z], [S1_x, S1_y, S1_z]]

    print(list_lin)
    return list_com, list_lin, list_foot_pos

def store_onecol_data(data, step, traj_size):

    CoM_x = []
    CoM_y = []
    CoM_z = []

    lin_x = []
    lin_y = []
    lin_z = []

    ang_x = []
    ang_y = []
    ang_z = []

    F0_x = []
    F0_y = []
    F0_z = []

    F1_x = []
    F1_y = []
    F1_z = []

    S0_x = []
    S0_y = []
    S0_z = []

    S1_x = []
    S1_y = []
    S1_z = []

    for pp in range(traj_size):
        CoM_x.append(data.iloc[pp*step + 0, 0])
        CoM_y.append(data.iloc[pp*step + 1, 0])
        CoM_z.append(data.iloc[pp*step + 2, 0])

        lin_x.append(data.iloc[pp*step + 3, 0])
        lin_y.append(data.iloc[pp*step + 4, 0])
        lin_z.append(data.iloc[pp*step + 5, 0])

        ang_x.append(data.iloc[pp*step + 6, 0])
        ang_y.append(data.iloc[pp*step + 7, 0])
        ang_z.append(data.iloc[pp*step + 8, 0])

        F0_x.append(data.iloc[pp*step + 9, 0])
        F0_y.append(data.iloc[pp*step + 10, 0])
        F0_z.append(data.iloc[pp*step + 11, 0])

        F1_x.append(data.iloc[pp*step + 12, 0])
        F1_y.append(data.iloc[pp*step + 13, 0])
        F1_z.append(data.iloc[pp*step + 14, 0])

        S0_x.append(data.iloc[pp*step + 15, 0])
        S0_y.append(data.iloc[pp*step + 16, 0])
        S0_z.append(data.iloc[pp*step + 17, 0])

        S1_x.append(data.iloc[pp*step + 18, 0])
        S1_y.append(data.iloc[pp*step + 19, 0])
        S1_z.append(data.iloc[pp*step + 20, 0])

    list_com_solution = [CoM_x, CoM_y, CoM_z]
    list_lin_solution = [lin_x, lin_y, lin_z]
    list_ang_solution = [ang_x, ang_y, ang_z]
    list_forces_solution = [[F0_x, F0_y, F0_z], [F1_x, F1_y, F1_z]]
    list_foot_pos_solution = [[S0_x, S0_y, S0_z], [S1_x, S1_y, S1_z]]

    return list_com_solution, list_lin_solution, list_ang_solution, list_forces_solution, list_foot_pos_solution

if __name__ == '__main__':

    step = 21

    # desired
    data_desired = pd.read_csv('../records/DesiredTrajectory100.csv', header=None)  # It has no header
    traj_size = data_desired.shape[0]-1 - 10
    list_com_desired, list_lin_desired, list_foot_pos_desired = store_multicol_data(data_desired, step, traj_size)

    # results
    data = pd.read_csv('../results/solution.csv', header=None)  # It has no header
    list_com_solution, list_lin_solution, list_ang_solution, list_forces_solution, \
        list_foot_pos_solution = store_onecol_data(data, step, traj_size)
    # CoM plot solution desired
    plt.plot(list_com_solution[0], label="x solution")
    plt.plot(list_com_solution[1], label="y solution")
    plt.plot(list_com_solution[2], label="z solution")

    plt.plot(list_com_desired[0], label="x desired")
    plt.plot(list_com_desired[1], label="y desired")
    plt.plot(list_com_desired[2], label="z desired")
    plt.legend()
    plt.title("CoM")
    plt.show()

    # Lin plot solution desired
    plt.plot(list_lin_solution[0], label="x solution")
    plt.plot(list_lin_solution[1], label="y solution")
    plt.plot(list_lin_solution[2], label="z solution")

    plt.plot(list_lin_desired[0], label="x desired")
    plt.plot(list_lin_desired[1], label="y desired")
    plt.plot(list_lin_desired[2], label="z desired")
    plt.legend()
    plt.title("Linear Vel.")
    plt.show()

    # Foot pos plot solution desired
    for leg_item in range(2):
        plt.plot(list_foot_pos_solution[leg_item][0], label="x solution")
        plt.plot(list_foot_pos_solution[leg_item][1], label="y solution")
        plt.plot(list_foot_pos_solution[leg_item][2], label="z solution")

        plt.plot(list_foot_pos_desired[leg_item][0], label="x desired")
        plt.plot(list_foot_pos_desired[leg_item][1], label="y desired")
        plt.plot(list_foot_pos_desired[leg_item][2], label="z desired")
        plt.legend()
        plt.title("Foot positions leg")
        plt.show()