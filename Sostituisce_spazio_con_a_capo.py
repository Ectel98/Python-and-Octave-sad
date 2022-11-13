f = open('intervalli risultanti a16\\a16_rr_interval.txt',"r+")

n_lines = sum(1 for line in f)

f.seek(0)

l = f.readlines()
f.seek(0)

for i in range (0,n_lines):
    f.write(l[i].replace(' ','\n'))

f.close()

