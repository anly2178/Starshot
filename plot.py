"""
Realised I was probably making something unnecessarily difficult. If I think
it's worthwhile, I will come back to this.
"""
# """
# General plotting functions.
# """
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator
# import numpy as np
#
# def plot_temp_vs_t(time, state, **kwargs):
#     """
#     Plots the temperature of the sail as a function of time. Also,
#     includes a label and title, shows the plot, and saves to a filepath
#     if given optional keyworded arguments. By default, will show plot.
#
#     Inputs:
#         time -- array of int or float (units: s)
#         state -- array of states including the temperature (units: K)
#         title -- string
#         label (optional) -- string
#         show -- bool (True or False)
#         filepath -- string
#
#         Example which does not use optional arguments:
#             plot_temp_vs_t(time, state)
#         Example which uses all optional arguments:
#             plot_temp_vs_t(time, state, title = 'Temperature vs Time',
#              label = 'absorption = 1', show = True, filepath = 'sample.png')
#     """
#     t = time
#     temp = state[2,:]
#     title = None
#     label = None
#     show = True
#     filepath = None
#
#     #Get optional arguments; raise error if wrong type.
#     for key,value in kwargs.items():
#         if key == 'title':
#             if type(value) != str:
#                 raise TypeError('title = <string>')
#             else:
#                 title = value
#         elif key == 'label':
#             if type(value) != str:
#                 raise TypeError('label = <string>')
#             else:
#                 label = value
#         elif key == 'show':
#             if type(value) != bool:
#                 raise TypeError('show = <bool>')
#             else:
#                 show = value
#         elif key == 'filepath':
#             if type(value) != str:
#                 raise TypeError('filepath = <string>')
#             else:
#                 filepath = value
#     #Plot
#
# #     plt.plot()
# # plt.plot(t, temp1, label = 'Kipping')
# # plt.xscale('log')
# # plt.xlabel('Time (s)')
# # plt.ylabel('Equilibrium temperature (K)')
# # plt.legend()
# # plt.savefig('Corrected/P0mp=1000/abs_coeff=1e-4/tempvst_diff.png')
# # state = np.array([[1],[2],[3]])
# # plot_temp_vs_t(1, state,title = 'fdas', label = 'fdsa', show = False, filepath = 'fdsa')
