# import matplotlib
# matplotlib.use("Qt4Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys
# from matplotlib import colors
# from matplotlib.ticker import PercentFormatter
def draw_bks(flname):
    with open(flname) as f:
        for l in f:
            lnlist = l.strip().split('\t')
            if lnlist[0] == "S":
                left_x = [] # left and right list
                right_x = []
                left_x_num = right_x_num = 0
                opt_elem_list = lnlist[4:]
                for e in opt_elem_list:
                    e_list = e.split(':')
                    if e_list[0] == "BK":
                        x_1 = int(e_list[2][1:])
                        x_2 = int(e_list[3])
                        # print x_1, x_2
                        if x_1 < x_2:
                            right_x.append(x_1)
                            right_x_num += 1
                        else:
                            left_x.append(x_1)
                            left_x_num += 1
                #plot here
                # print left_x, right_x
                plot_left_right_bks([left_x, right_x], int(lnlist[2]), left_x_num, right_x_num) 
                #wait to process another one
                # k = raw_input("Press q to quit else to continue...")
                # if k == "q":
                    # break
            else:
                continue
    f.close()

def plot_left_right_bks(x, s_len, l_n, r_n):
    # print x
    left_x = x[0]
    right_x = x[1]
    w = 250 
    
    left_counts, left_bins = np.histogram(left_x, bins=w, range=(0,s_len))
    right_counts, right_bins = np.histogram(right_x, bins=w, range = (0,s_len))
    left_counts_list = np.negative(left_counts).tolist()
    right_counts_list = right_counts.tolist()
    x = left_bins.tolist()
    y_min = min(left_counts_list)
    y_max = max(right_counts_list)
     
    xb = np.arange(w)
    fig = plt.figure(figsize=(20,12))
    ax  = plt.subplot(111)
    # plt.get_current_fig_manager().window.setGeometry(3,3,12,20)
    ax.bar(xb , left_counts_list, color='g')
    ax.bar(xb, right_counts_list, color='b')
    # ax.set_xlim(0, s_len)
    for i, v in enumerate(left_counts_list):
        if v == y_min:
            ax.text(i, v-1, str(v), color='b', fontweight='bold')
    for i, v in enumerate(right_counts_list):
        if v == y_max:
            ax.text(i, v+1, str(v), color='b', fontweight='bold')
    ax.set_ylim(y_min, y_max)
    # ax.set_ylim(-8000, 8000)
    ax.set_ylabel("count"+" MLB: "+str(-y_min)+" MRB: "+str(y_max))
    ax.set_xlabel("RL: "+str(s_len)+"\tLB: "+str(l_n)+" RB: "+str(r_n))
    # plt.figure(figsize=(1,1))
    # ax.set_xticks(np.linspace(0,s_len, num=50))
    # plt.bar(xb, left_counts_list, color='b', width=0.25)
    # plt.bar(xb, right_counts_list, color='g', width=0.25)
    # plt.bar(x[:-1], left_counts_list, color='b', width=0.25)
    # plt.bar(x[:-1], right_counts_list, color='g', width=0.25)
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    draw_bks(sys.argv[1])        


