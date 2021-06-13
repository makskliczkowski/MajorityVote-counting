import os
import numpy as np
import pandas as pd
from math import sqrt


directory = "C:\\Users\\Wiktor\\Desktop\\tak\\Studia\\Big Data Analytics\\1 semestr\\Introduction to Complex System\\" \
            "MaksRepo\\MajorityVote-counting\\cpp\\Project1\\Results\\Non_infty\\Triangular\\"
directory_save = "C:\\Users\\Wiktor\\Desktop\\tak\\Studia\\Big Data Analytics\\1 semestr\\Introduction to " \
                 "Complex System\\MaksRepo\\MajorityVote-counting\\cpp\\Project1\\Results\\Non_infty\\Triangular_format\\"

for file in os.listdir(directory):
    if file[-4:] == ".dat":
        data_m = []
        data_sus = []
        data_bin = []
        f = pd.read_csv(directory+"\\"+file, names=["p", 'm', 'susc', 'binder'], skiprows=1, sep="\t")
        f = f.astype(np.float32)
        splt = file.split(",")
        q = float(splt[1].split("=")[-1])
        N = int(splt[2].split("=")[-1])
        L = sqrt(N)
        ps = f['p']
        m = f['m']
        sus = f['susc']
        binder = f['binder']
        file_m = open(directory_save+"\\mmmm"+file, 'w')
        file_s = open(directory_save+"\\ssss"+file, 'w')
        file_b = open(directory_save+"\\bbbb"+file, 'w')

        for i in range(len(ps)):
            file_m.write(str(ps[i]) + "\t" + str(m[i]) + "\t" + str(m[i]/3000) + "\n")
            file_s.write(str(ps[i]) + "\t" + str(sus[i]) + "\t" + str(sus[i]/3000) + "\n")
            file_b.write(str(ps[i]) + "\t" + str(binder[i]) + "\t" + str(binder[i]/3000) + "\n")
            data_m.append([ps[i], m[i], m[i]/3000])
            data_m.append([ps[i], sus[i], sus[i]/3000])
            data_m.append([ps[i], binder[i], binder[i]/3000])

        file_m.close()
        file_s.close()
        file_b.close()
        print(splt)
        print("\n")
