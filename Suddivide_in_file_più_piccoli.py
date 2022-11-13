from pathlib import Path


for i in range(2,10):
    s = "c0" + str(i)
    in_f = open(s + '.txt',"r")
    n_lines = sum(1 for line in in_f)
    in_f.seek(0)
    st = 0

    for e in range(0,round(n_lines/10000)):
        Path('divisi ' + s).mkdir(parents=True, exist_ok=True)
        out_f = open('divisi ' + s +'\\' + s +'-' + str(round((st+1)/10000)) + '.txt',"w")
        for i in range (st,st+10000):
            out_f.write(in_f.readline())
        st += 10000
        out_f.close()
    in_f.close()
