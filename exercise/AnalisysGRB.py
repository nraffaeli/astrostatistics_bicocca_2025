import numpy as np

data = np.loadtxt("Summary_GRB.txt", dtype='str',unpack='True')

T90 = np.array(data[6],dtype=float) 
T90_err = np.array(data[7],dtype=float)
Red = np.array(data[11], dtype=float)
Red_err = np.array(data[12], dtype=float)
Fluence = np.array(data[9], dtype=float)
Fluence_err = np.array(data[10], dtype=float)

for i in range(len(T90)) :
    if T90[i] <= -998 :
        T90[i] = np.nan
        Red[i] = np.nan
        T90_err[i] = np.nan
        Red_err[i] = np.nan
    elif Red[i] <= -998 :
        T90[i] = np.nan
        Red[i] = np.nan
        T90_err[i] = np.nan
        Red_err[i] = np.nan
    
valid_indices = ~np.isnan(T90) & ~np.isnan(Red)

with open("GRB_Red_T90.txt", "w") as f:
    for i in np.where(valid_indices)[0]:
        f.write(f"{Red[i]}\t{T90[i]}\t{Red_err[i]}\t{T90_err[i]}\n")

for i in range(len(T90)) :
    if T90[i] <= -998 :
        T90[i] = np.nan
        Fluence[i] = np.nan
        T90_err[i] = np.nan
        Fluence_err[i] = np.nan
    elif Fluence[i] <= -998 :
        T90[i] = np.nan
        Fluence[i] = np.nan
        T90_err[i] = np.nan
        Fluence_err[i] = np.nan

valid_indices = ~np.isnan(T90) & ~np.isnan(Fluence)

with open("GRB_Fluence_T90.txt", "w") as f:
    for i in np.where(valid_indices)[0]:
        f.write(f"{Fluence[i]}\t{T90[i]}\t{Fluence_err[i]}\t{T90_err[i]}\n")
